"""
CLI Progress Monitor using curses for real-time algorithm execution monitoring.
"""

import argparse
import curses
import signal
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import List

from src.infrastructure.persistence.work_state.queries import (
    ExecutionDetail,
    ProgressSummary,
    WorkStateQueries,
)


def show_loading_progress(message: str, progress: float, width: int = 50):
    """
    Mostra uma barra de progresso durante o carregamento.
    
    Args:
        message: Mensagem a ser exibida
        progress: Progresso entre 0.0 e 1.0
        width: Largura da barra de progresso
    """
    import sys
    
    if progress >= 1.0:
        # Quando completo, substitui a linha atual com a mensagem final
        sys.stdout.write(f'\r{message}' + ' ' * (width + 10))
        sys.stdout.flush()
    else:
        # Durante o progresso, mostra a barra
        filled_len = int(width * progress)
        bar = '█' * filled_len + '░' * (width - filled_len)
        percent = int(progress * 100)
        
        # Limpa a linha atual e mostra a barra
        sys.stdout.write(f'\r{message}: [{bar}] {percent}%')
        sys.stdout.flush()


def wait_for_work_running_status(work_id: str, queries: WorkStateQueries, timeout: int = 30) -> bool:
    """
    Aguarda especificamente que o work mude para status 'running'.

    Args:
        work_id: ID do work a ser monitorado
        queries: Instância do WorkStateQueries
        timeout: Tempo limite em segundos para aguardar

    Returns:
        True se o work estiver em execução, False se timeout
    """
    start_time = time.time()
    
    while time.time() - start_time < timeout:
        try:
            status = queries.get_work_status(work_id)
            
            if status == 'running':
                return True
            elif status in ['completed', 'failed', 'canceled', 'error']:
                return True  # Permite mostrar o status final
            
        except Exception as e:
            pass
            
        time.sleep(2)

    return False


def wait_for_work_service_initialization(timeout: int = 30) -> bool:
    """
    Aguarda a inicialização do WorkServicePersistence.

    Args:
        timeout: Tempo limite em segundos para aguardar a inicialização

    Returns:
        True se inicializado com sucesso, False se timeout
    """
    start_time = time.time()

    while time.time() - start_time < timeout:
        try:
            # Tenta obter o work service
            from src.application.services.work_service import get_work_service
            work_service = get_work_service()

            # Verifica se o repositório está pronto
            if hasattr(work_service, 'repository') and hasattr(work_service.repository, 'is_ready'):
                if work_service.repository.is_ready():
                    # Testa uma operação básica para confirmar funcionamento
                    work_service.list()
                    return True
            else:
                # Fallback: testa operação básica diretamente
                work_service.list()
                return True

        except Exception as e:
            # Se houver erro, aguarda um pouco e tenta novamente
            time.sleep(1)

    return False


def wait_for_work_and_database(work_id: str, state_db_path: Path, timeout: int = 60) -> bool:
    """
    Aguarda que o banco de estado seja criado, a tabela work seja escrita e o work esteja em execução.

    Args:
        work_id: ID do work a ser monitorado
        state_db_path: Caminho para o banco de dados de estado
        timeout: Tempo limite em segundos para aguardar

    Returns:
        True se tudo estiver pronto, False se timeout
    """
    start_time = time.time()
    
    # Primeiro aguarda a criação do arquivo do banco
    while time.time() - start_time < timeout and not state_db_path.exists():
        time.sleep(1)

    if not state_db_path.exists():
        return False

    # Agora aguarda a inicialização do banco e criação da tabela work
    while time.time() - start_time < timeout:
        try:            
            queries = WorkStateQueries(str(state_db_path))
            
            # Verifica se o work existe no banco (tabela work foi criada)
            if queries.work_exists(work_id):
                # Verifica se o work está em um estado válido (não mais 'unknown' inicial)
                status = queries.get_work_status(work_id)
                
                if status and status in ['running', 'queued']:
                    # Aguarda povoamento do banco com dados de progresso
                    progress_data = queries.get_work_progress_summary(work_id)
                    if progress_data and progress_data.global_execution.get('Total', 0) > 0:
                        queries.close()
                        return True
                
                elif status in ['completed', 'failed', 'canceled', 'error']:
                    queries.close()
                    return True  # Permite mostrar o status final
            
            queries.close()
            
        except Exception as e:
            # Se ainda não conseguiu conectar ou ler do banco, continua tentando
            pass
            
        time.sleep(1)

    return False


class ProgressMonitor:
    """Real-time progress monitor using curses, htop-style."""

    def __init__(self, work_id: str):
        """Initialize monitor."""
        print(f"🚀 Inicializando monitor para work: {work_id}")
        
        # Inicia barra de progresso única
        total_steps = 5
        current_step = 0
        
        # Etapa 1: Inicialização do WorkService
        current_step += 1
        show_loading_progress("⚙️ Inicializando monitor", current_step / total_steps)
        
        if not wait_for_work_service_initialization():
            raise RuntimeError("Timeout aguardando inicialização do WorkService")

        from src.application.services.work_service import get_work_service

        self.work_id = work_id
        self.work_service = get_work_service()
        self.running = True
        self.refresh_interval = 0.1  # seconds

        # Etapa 2: Buscar detalhes do work
        current_step += 1
        show_loading_progress("⚙️ Inicializando monitor", current_step / total_steps)
        
        max_retries = 10
        retry_count = 0
        work_details = None

        while retry_count < max_retries:
            work_details = self.work_service.get(work_id)
            if work_details and work_details.get("output_path"):
                break
            time.sleep(1)
            retry_count += 1

        if not work_details or not work_details.get("output_path"):
            raise ValueError(f"Work {work_id} not found or has no output path")

        # Etapa 3: Preparar banco de dados
        current_step += 1
        show_loading_progress("⚙️ Inicializando monitor", current_step / total_steps)
        
        state_db_path = Path(work_details["output_path"]) / "state.db"

        # Etapa 4: Aguardar dados de execução
        current_step += 1
        show_loading_progress("⚙️ Inicializando monitor", current_step / total_steps)
        
        if not wait_for_work_and_database(work_id, state_db_path, timeout=60):
            raise ValueError(f"Timeout aguardando criação do banco de estado e dados de execução")

        # Etapa 5: Finalizar inicialização
        current_step += 1
        show_loading_progress("⚙️ Inicializando monitor", current_step / total_steps)
        
        self.queries = WorkStateQueries(str(state_db_path))
        self.start_time = time.time()
        self.scroll_pos = 0
        self.executions: List[ExecutionDetail] = []
        self.logs: List[dict] = []

        # Captura o status inicial para detectar mudanças
        try:
            self.previous_status = self.queries.get_work_status(self.work_id)
        except Exception:
            self.previous_status = None

        # Finalizar barra de progresso
        show_loading_progress("✅ Monitor inicializado com sucesso", 1.0)
        print()  # Nova linha após finalizar
        print("🎯 Iniciando interface...")

    def signal_handler(self, signum, frame):
        """Handle Ctrl+C gracefully by canceling the work."""
        try:
            self.work_service.cancel(self.work_id)
        except Exception:
            pass
        self.running = False

    def format_progress_bar(self, progress: float, width: int, show_percent: bool = True) -> str:
        """Creates a text-based progress bar."""
        filled_len = int(width * progress)
        bar = '█' * filled_len + '─' * (width - filled_len)
        if show_percent:
            return f"[{bar}] {progress:.1%}"
        else:
            return bar

    def format_execution_time(self, elapsed_seconds: float) -> str:
        """Format execution time as mm:ss.mmm"""
        total_ms = int(elapsed_seconds * 1000)
        minutes = total_ms // 60000
        seconds = (total_ms % 60000) // 1000
        milliseconds = total_ms % 1000

        if minutes > 0:
            return f"{minutes:02d}:{seconds:02d}.{milliseconds:03d}"
        else:
            return f"{seconds:02d}.{milliseconds:03d}s"

    def draw_header(self, stdscr, progress: ProgressSummary, work_status: str):
        """Draws the header with summary information."""
        h, w = stdscr.getmaxyx()
        if h < 5:
            return 0

        # Line 1: Work Info
        uptime = timedelta(seconds=int(time.time() - self.start_time))
        work_info = f"Work: {self.work_id}  Status: {work_status.upper()}  Uptime: {uptime}"
        stdscr.addstr(0, 1, work_info.ljust(w - 2), curses.color_pair(4) | curses.A_BOLD)

        # Line 2: Global Progress
        glob_prog = progress.global_progress
        glob_label = f"Overall Progress: {progress.global_execution['Finished']}/{progress.global_execution['Total']}"

        # Calculate available width for progress bar
        available_width = w - len(glob_label) - 10  # Reserve space for label + padding
        bar_width = max(10, min(40, available_width))  # Between 10-40 chars, but fit screen

        glob_bar = self.format_progress_bar(glob_prog, bar_width)

        stdscr.addstr(1, 1, glob_label, curses.A_BOLD)
        if len(glob_label) + len(glob_bar) + 4 <= w - 2:
            stdscr.addstr(1, len(glob_label) + 2, glob_bar, curses.color_pair(3))
        else:
            # If still too long, show on next line
            stdscr.addstr(2, 3, glob_bar, curses.color_pair(3))
            # Adjust return value to account for extra line
            return 7

        # Line 3-4: Hierarchical Progress in aligned columns
        tasks = progress.tasks
        dsets = progress.datasets
        cfgs = progress.configs
        algs = progress.algorithms

        task_running_count = 1 if tasks.get('Running') else 0
        dset_running_count = 1 if dsets.get('Running') else 0
        cfg_running_count = 1 if cfgs.get('Running') else 0
        alg_running_count = 1 if algs.get('Running') else 0

        # Calculate column widths for alignment
        col1_width = w // 2 - 2  # Half screen minus padding
        col2_width = w - col1_width - 4  # Remaining width

        # Line 3: Task and Dataset
        task_str = f"Task: {tasks.get('Running', '-') or '-'} ({tasks['Finished'] + task_running_count}/{tasks['Total']})"
        dset_str = f"Dataset: {dsets.get('Running', '-') or '-'} ({dsets['Finished'] + dset_running_count}/{dsets['Total']})"

        stdscr.addstr(2, 1, task_str[:col1_width].ljust(col1_width))
        stdscr.addstr(2, col1_width + 2, dset_str[:col2_width])

        # Line 4: Config and Algorithm
        cfg_str = f"Config: {cfgs.get('Running', '-') or '-'} ({cfgs['Finished'] + cfg_running_count}/{cfgs['Total']})"
        alg_str = f"Algorithm: {algs.get('Running', '-') or '-'} ({algs['Finished'] + alg_running_count}/{algs['Total']})"

        stdscr.addstr(3, 1, cfg_str[:col1_width].ljust(col1_width))
        stdscr.addstr(3, col1_width + 2, alg_str[:col2_width])

        # Line 5: Current Combination Description
        combo = progress.current_combination_details
        if combo:
            exec_info = progress.execution
            current_seq = exec_info['Finished'] + exec_info['Running']
            total_seq = exec_info['Total']
            combo_desc = f"Running combination: {combo['task_id']} > {combo['dataset_id']} > {combo['preset_id']}. Sequences: {current_seq}/{total_seq}"
            stdscr.addstr(4, 1, combo_desc.ljust(w-2)[:w-2], curses.A_BOLD)

        # Line 6: Execution Table Header
        header = f"{'ST SEQ':<6} {'ALGORITHM':<19} {'PROGRESS':<17} {'TIME':<10} {'DETAILS'}"
        stdscr.addstr(5, 0, header.ljust(w), curses.color_pair(5) | curses.A_REVERSE)

        return 6

    def draw_executions_list(self, stdscr, y_pos: int):
        """Draws the list of executions with visual status indicators."""
        h, w = stdscr.getmaxyx()
        max_rows = h - y_pos - 1

        for i in range(max_rows):
            exec_idx = self.scroll_pos + i
            if exec_idx >= len(self.executions):
                break

            exec_detail = self.executions[exec_idx]

            seq = exec_detail.sequencia
            algo = exec_detail.algorithm_id[:19]

            # Visual status indicator
            if exec_detail.status == 'completed':
                status_icon = "✓"
                color = curses.color_pair(4)  # Cyan for completed
                progress_bar = f"[{'█' * 15}]"  # Full bar for completed
            elif exec_detail.status == 'running':
                status_icon = "▶"
                color = curses.color_pair(3) if exec_detail.progress > 0 else curses.color_pair(2)  # Green/Yellow
                progress_bar = f"[{self.format_progress_bar(exec_detail.progress, 15, show_percent=False)}]"
            elif exec_detail.status in ('failed', 'error'):
                status_icon = "✗"
                color = curses.color_pair(1)  # Red for failed
                progress_bar = f"[{'▓' * 15}]"  # Error pattern
            elif exec_detail.status == 'queued':
                status_icon = "⋯"
                color = curses.color_pair(6) if hasattr(curses, 'COLOR_MAGENTA') else curses.color_pair(2)  # Magenta or fallback
                progress_bar = f"[{'─' * 15}]"  # Empty bar for queued
            else:
                status_icon = "?"
                color = curses.color_pair(2)  # Yellow for unknown
                progress_bar = f"[{self.format_progress_bar(exec_detail.progress, 15, show_percent=False)}]"

            # Calculate elapsed time only if execution has started
            if exec_detail.finished_at and exec_detail.started_at:
                elapsed_time = exec_detail.finished_at - exec_detail.started_at
                time_str = self.format_execution_time(elapsed_time)
            elif exec_detail.started_at:
                elapsed_time = time.time() - exec_detail.started_at
                time_str = self.format_execution_time(elapsed_time)
            else:
                # No time display for queued executions
                time_str = "-"

            details = (exec_detail.progress_message or '')[:w - 70]

            line = f"{status_icon} {seq:<4} {algo:<19} {progress_bar:<17} {time_str:<10} {details}"

            stdscr.addstr(y_pos + i, 0, line.ljust(w)[:w], color)

    def draw_footer(self, stdscr):
        """Draws the footer with controls, legend and status."""
        h, w = stdscr.getmaxyx()
        if h <= 1:
            return

        controls = "Q:Quit P:Pause C:Cancel ↑↓:Scroll"
        legend = "✓:Done ▶:Running ✗:Failed ⋯:Queued"
        timestamp = datetime.now().strftime("%H:%M:%S")

        # If screen is wide enough, show legend and controls on same line
        if len(controls) + len(legend) + len(timestamp) + 6 <= w:
            footer_text = f"{controls}  {legend}".ljust(w - len(timestamp) - 1) + timestamp
        else:
            # Otherwise just show controls and timestamp
            footer_text = f"{controls}".ljust(w - len(timestamp) - 1) + timestamp

        stdscr.addstr(h - 1, 0, footer_text, curses.color_pair(5) | curses.A_REVERSE)

    def draw_log_panel(self, stdscr, y_pos: int, max_height: int):
        """Draws a small panel with recent error and warning logs."""
        h, w = stdscr.getmaxyx()
        if max_height <= 1 or not self.logs:
            return

        # Small panel header
        title = " Recent Events "
        stdscr.addstr(y_pos, 0, "┌" + "─" * (len(title)) + "┐" + "─" * (w - len(title) - 3), curses.color_pair(5))
        stdscr.addstr(y_pos, 1, title, curses.color_pair(5) | curses.A_REVERSE)
        y_pos += 1

        # Show only the most recent 2 entries
        for i, log_entry in enumerate(self.logs[:min(max_height - 1, 2)]):
            if y_pos >= h - 1:
                break

            log_time = datetime.fromtimestamp(log_entry['timestamp']).strftime('%H:%M:%S')

            if log_entry['type'] == 'error':
                color = curses.color_pair(1)
                icon = "✗"
            else:
                color = curses.color_pair(2)
                icon = "!"

            # Truncate message to fit width
            max_msg_len = w - 15  # Reserve space for time and icon
            message = log_entry['message']
            if len(message) > max_msg_len:
                message = message[:max_msg_len-3] + "..."

            line = f"│ {icon} [{log_time}] {message}"
            stdscr.addstr(y_pos + i, 0, line.ljust(w)[:w], color)

    def run(self, stdscr):
        """Main monitoring loop."""
        curses.curs_set(0)
        stdscr.nodelay(1)
        stdscr.timeout(50) # 50ms timeout for getch()

        curses.start_color()
        curses.init_pair(1, curses.COLOR_RED, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)
        curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)
        curses.init_pair(4, curses.COLOR_CYAN, curses.COLOR_BLACK)
        curses.init_pair(5, curses.COLOR_WHITE, curses.COLOR_BLUE)
        curses.init_pair(6, curses.COLOR_MAGENTA, curses.COLOR_BLACK)

        last_refresh = 0

        while self.running:
            h, w = stdscr.getmaxyx()
            current_time = time.time()

            key = stdscr.getch()
            if key in [ord('q'), ord('Q'), 27]:
                break
            elif key in [ord('p'), ord('P')]:
                self.handle_pause(stdscr)
            elif key in [ord('c'), ord('C')]:
                self.handle_cancel(stdscr)
            elif key == curses.KEY_UP:
                self.scroll_pos = max(0, self.scroll_pos - 1)
            elif key == curses.KEY_DOWN:
                # Adjust scroll calculation for new layout
                log_panel_height = 3 if self.logs else 0  # Fixed small size
                max_scroll = max(0, len(self.executions) - (h - 6 - log_panel_height - 1))
                self.scroll_pos = min(max_scroll, self.scroll_pos + 1)

            # Verifica se foi pausado ou cancelado e deve sair
            if not self.running:
                break

            if current_time - last_refresh >= self.refresh_interval:
                try:
                    stdscr.clear()

                    progress = self.queries.get_work_progress_summary(self.work_id)
                    if not progress:
                        stdscr.addstr(0, 0, f"Work '{self.work_id}' not found!", curses.color_pair(1))
                        stdscr.refresh()
                        time.sleep(1)
                        continue

                    work_status = self.queries.get_work_status(self.work_id) or 'unknown'

                    # Detecta mudança de status de 'queued' para algo que não seja 'running'
                    if (self.previous_status == 'queued' and
                        work_status not in ['queued', 'running'] and
                        work_status != 'unknown'):
                        self.show_status_change_screen(stdscr, self.previous_status, work_status, progress)
                        break

                    # Atualiza o status anterior
                    if work_status != 'unknown':
                        self.previous_status = work_status

                    if work_status == 'paused':
                        self.show_paused_screen(stdscr, work_status, progress)
                        # Exit monitor when paused
                        break
                    elif work_status not in ['running', 'queued']:
                        self.show_final_screen(stdscr, work_status, progress)
                        break

                    # Fetch executions for the current combination
                    current_combo_id = progress.current_combination_details.get('combination_id') if progress.current_combination_details else None
                    if current_combo_id:
                        old_executions_count = len(self.executions)
                        # Get ALL executions for the combination, not just running ones
                        self.executions = self.queries.get_combination_executions_detail(current_combo_id)
                        # Sort by sequence number for better visualization
                        self.executions.sort(key=lambda x: x.sequencia)
                        new_executions_count = len(self.executions)
                    else:
                        old_executions_count = 0
                        new_executions_count = 0
                        self.executions = []

                    # Fetch and combine logs
                    errors = self.queries.get_error_summary(self.work_id, limit=2)
                    warnings = self.queries.get_execution_warnings(self.work_id, limit=2)

                    logs = []
                    for e in errors:
                        logs.append({'type': 'error', 'timestamp': e.timestamp, 'message': f"Unit {e.unit_id[:8]}: {e.error_message}"})
                    for w in warnings:
                        logs.append({'type': 'warning', 'timestamp': w['timestamp'], 'message': f"Unit {(w.get('unit_id') or 'N/A')[:8]}: {w.get('message', 'Unknown warning')}"})
                    self.logs = sorted(logs, key=lambda x: x.get('timestamp', 0), reverse=True)[:2]  # Keep only 2 most recent

                    # Calculate layout - fixed small log panel
                    log_panel_height = 3 if self.logs else 0 # 1 header + max 2 log lines
                    exec_list_height = h - 6 - log_panel_height - 1  # header - logs - footer

                    # Smart auto-scroll logic - focus on running/next executions
                    if exec_list_height > 0 and new_executions_count > 0:
                        max_scroll = max(0, new_executions_count - exec_list_height)

                        # Find first running or queued execution
                        first_active_idx = None
                        for i, exec_detail in enumerate(self.executions):
                            if exec_detail.status in ('running', 'queued'):
                                first_active_idx = i
                                break

                        # Check if we were showing the last line before update
                        was_at_bottom = (old_executions_count <= exec_list_height) or \
                                       (self.scroll_pos >= max(0, old_executions_count - exec_list_height))

                        # If new executions were added and we were at bottom, or if no manual scrolling occurred
                        if new_executions_count > old_executions_count and was_at_bottom:
                            # Try to center on first active execution, but don't go past bottom
                            if first_active_idx is not None:
                                target_scroll = max(0, min(max_scroll, first_active_idx - exec_list_height // 2))
                                self.scroll_pos = target_scroll
                            else:
                                self.scroll_pos = max_scroll
                        elif self.scroll_pos > max_scroll:
                            # Ensure scroll position is valid
                            self.scroll_pos = max_scroll

                    # Draw interface
                    y_pos = self.draw_header(stdscr, progress, work_status)

                    if exec_list_height > 0:
                        self.draw_executions_list(stdscr, y_pos)

                    # Draw small log panel just above footer
                    if self.logs:
                        log_y_pos = h - 1 - log_panel_height
                        self.draw_log_panel(stdscr, log_y_pos, log_panel_height)

                    self.draw_footer(stdscr)

                    stdscr.refresh()
                    last_refresh = current_time

                except Exception as e:
                    stdscr.clear()
                    stdscr.addstr(0, 0, f"Error: {str(e)}", curses.color_pair(1))
                    stdscr.addstr(1, 0, "Press 'q' to quit")
                    stdscr.refresh()
                    time.sleep(1)

            time.sleep(0.02)

    def handle_pause(self, stdscr):
        """Shows a confirmation dialog for pausing the work."""
        h, w = stdscr.getmaxyx()
        
        # Criar caixa de diálogo com fundo
        dialog_width = 50
        dialog_height = 7
        start_y = (h - dialog_height) // 2
        start_x = (w - dialog_width) // 2
        
        # Criar janela de diálogo
        dialog_win = curses.newwin(dialog_height, dialog_width, start_y, start_x)
        dialog_win.bkgd(' ', curses.color_pair(5))  # Fundo azul
        dialog_win.box()
        
        # Título
        title = "⏸️ PAUSAR TRABALHO"
        dialog_win.addstr(1, (dialog_width - len(title)) // 2, title, curses.A_BOLD)
        
        # Mensagem
        msg1 = f"Pausar work '{self.work_id[:20]}...'?"
        msg2 = "Pressione Y para confirmar ou N para cancelar"
        dialog_win.addstr(3, (dialog_width - len(msg1)) // 2, msg1)
        dialog_win.addstr(4, (dialog_width - len(msg2)) // 2, msg2)
        
        # Opções
        options = "[Y] Sim    [N] Não"
        dialog_win.addstr(5, (dialog_width - len(options)) // 2, options, curses.A_BOLD)
        
        dialog_win.refresh()

        confirm_key = -1
        while confirm_key not in [ord('y'), ord('Y'), ord('n'), ord('N'), 27]:  # 27 = ESC
            confirm_key = stdscr.getch()
            time.sleep(0.1)

        if confirm_key in [ord('y'), ord('Y')]:
            try:
                self.work_service.pause(self.work_id)
                # Mostrar tela de confirmação de pausa
                self.show_paused_confirmation_screen(stdscr)
                # Após pausar com sucesso, sair do monitor
                self.running = False
            except Exception:
                pass # Ignore errors, will be reflected in status

    def show_paused_confirmation_screen(self, stdscr):
        """Displays a screen when the work is successfully paused."""
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = "⏸️ TRABALHO PAUSADO"
        msg1 = f"O work '{self.work_id}' foi pausado com sucesso."
        msg2 = "Monitor será finalizado em 3 segundos..."

        stdscr.addstr(h//2 - 1, (w - len(title)) // 2, title, curses.color_pair(2) | curses.A_BOLD)
        stdscr.addstr(h//2, (w - len(msg1)) // 2, msg1)
        stdscr.addstr(h//2 + 1, (w - len(msg2)) // 2, msg2)
        stdscr.refresh()

        for _ in range(30):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def show_paused_screen(self, stdscr, status: str, progress):
        """Displays a paused screen when the work is paused."""
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = f"⏸️ Work Pausado: {status.upper()}"
        finished = progress.global_execution["Finished"]
        total = progress.global_execution["Total"]
        pause_progress = f"Progresso no momento da pausa: {finished}/{total} execuções completadas."

        stdscr.addstr(h//2 - 1, (w - len(title)) // 2, title, curses.color_pair(2) | curses.A_BOLD)
        stdscr.addstr(h//2, (w - len(pause_progress)) // 2, pause_progress)
        stdscr.addstr(h//2 + 1, (w - 60) // 2, "O trabalho foi pausado com sucesso.")
        stdscr.addstr(h//2 + 2, (w - 40) // 2, "Monitor será finalizado em 3 segundos...")
        stdscr.refresh()

        for _ in range(30):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def handle_cancel(self, stdscr):
        """Shows a confirmation dialog for cancelling the work."""
        h, w = stdscr.getmaxyx()
        
        # Criar caixa de diálogo com fundo
        dialog_width = 50
        dialog_height = 7
        start_y = (h - dialog_height) // 2
        start_x = (w - dialog_width) // 2
        
        # Criar janela de diálogo
        dialog_win = curses.newwin(dialog_height, dialog_width, start_y, start_x)
        dialog_win.bkgd(' ', curses.color_pair(1))  # Fundo vermelho
        dialog_win.box()
        
        # Título
        title = "🛑 CANCELAR TRABALHO"
        dialog_win.addstr(1, (dialog_width - len(title)) // 2, title, curses.A_BOLD)
        
        # Mensagem
        msg1 = f"Cancelar work '{self.work_id[:20]}...'?"
        msg2 = "Pressione Y para confirmar ou N para cancelar"
        dialog_win.addstr(3, (dialog_width - len(msg1)) // 2, msg1)
        dialog_win.addstr(4, (dialog_width - len(msg2)) // 2, msg2)
        
        # Opções
        options = "[Y] Sim    [N] Não"
        dialog_win.addstr(5, (dialog_width - len(options)) // 2, options, curses.A_BOLD)
        
        dialog_win.refresh()

        confirm_key = -1
        while confirm_key not in [ord('y'), ord('Y'), ord('n'), ord('N'), 27]:  # 27 = ESC
            confirm_key = stdscr.getch()
            time.sleep(0.1)

        if confirm_key in [ord('y'), ord('Y')]:
            try:
                self.work_service.cancel(self.work_id)
                # Mostrar tela de confirmação de cancelamento
                self.show_canceled_screen(stdscr)
                # Após cancelar com sucesso, sair do monitor
                self.running = False
            except Exception:
                pass # Ignore errors, will be reflected in status

    def show_canceled_screen(self, stdscr):
        """Displays a screen when the work is successfully canceled."""
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = "🛑 TRABALHO CANCELADO"
        msg1 = f"O work '{self.work_id}' foi cancelado com sucesso."
        msg2 = "Monitor será finalizado em 3 segundos..."

        stdscr.addstr(h//2 - 1, (w - len(title)) // 2, title, curses.color_pair(1) | curses.A_BOLD)
        stdscr.addstr(h//2, (w - len(msg1)) // 2, msg1)
        stdscr.addstr(h//2 + 1, (w - len(msg2)) // 2, msg2)
        stdscr.refresh()

        for _ in range(30):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def show_status_change_screen(self, stdscr, old_status: str, new_status: str, progress: ProgressSummary):
        """Displays a screen when work status changes from queued to non-running status."""
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        # Determina a cor e mensagem baseada no novo status
        if new_status in ['failed', 'error']:
            color = curses.color_pair(1)  # Vermelho
            title = f"Work Failed: {new_status.upper()}"
            message = "O trabalho falhou antes de começar a executar."
        elif new_status == 'canceled':
            color = curses.color_pair(2)  # Amarelo
            title = f"Work Canceled: {new_status.upper()}"
            message = "O trabalho foi cancelado antes de começar a executar."
        elif new_status == 'completed':
            color = curses.color_pair(3)  # Verde
            title = f"Work Completed: {new_status.upper()}"
            message = "O trabalho foi concluído rapidamente."
        else:
            color = curses.color_pair(2)  # Amarelo
            title = f"Status Changed: {old_status.upper()} → {new_status.upper()}"
            message = f"O trabalho mudou de '{old_status}' para '{new_status}' sem executar."

        finished = progress.global_execution["Finished"]
        total = progress.global_execution["Total"]
        progress_info = f"Progresso: {finished}/{total} execuções completadas."

        # Calcula posições centralizadas
        title_y = h//2 - 2
        message_y = h//2 - 1
        progress_y = h//2
        exit_y = h//2 + 2

        stdscr.addstr(title_y, (w - len(title)) // 2, title, color | curses.A_BOLD)
        stdscr.addstr(message_y, (w - len(message)) // 2, message)
        stdscr.addstr(progress_y, (w - len(progress_info)) // 2, progress_info)
        stdscr.addstr(exit_y, (w - 40) // 2, "Monitor will exit in 5 seconds...")
        stdscr.refresh()

        # Aguarda 5 segundos ou tecla pressionada
        for _ in range(50):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def show_final_screen(self, stdscr, status: str, progress: ProgressSummary):
        """Displays a final summary screen when the work is done."""
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        # Escolhe emoji baseado no status
        if status == 'completed':
            emoji = "🎉"
        elif status in ['failed', 'error']:
            emoji = "❌"
        elif status == 'canceled':
            emoji = "🛑"
        else:
            emoji = "ℹ️"

        title = f"{emoji} Work Finished with status: {status.upper()}"
        finished = progress.global_execution["Finished"]
        total = progress.global_execution["Total"]
        final_progress = f"Final Progress: {finished}/{total} executions completed."

        stdscr.addstr(h//2 - 1, (w - len(title)) // 2, title, curses.color_pair(3) | curses.A_BOLD)
        stdscr.addstr(h//2, (w - len(final_progress)) // 2, final_progress)
        stdscr.addstr(h//2 + 2, (w - 30) // 2, "Monitor will exit in 3 seconds...")
        stdscr.refresh()

        for _ in range(30):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def start(self):
        """Start the monitor."""
        signal.signal(signal.SIGINT, self.signal_handler)
        try:
            if not self.queries.work_exists(self.work_id):
                print(f"Error: Work '{self.work_id}' not found!")
                return 1
            curses.wrapper(self.run)
            return 0
        except KeyboardInterrupt:
            return 0
        except Exception as e:
            print(f"Error starting monitor: {e}")
            return 1
        finally:
            self.queries.close()


def main():
    """CLI entry point."""

    parser = argparse.ArgumentParser(description="CSPBench Progress Monitor (htop-style)")
    parser.add_argument("work_id", help="Work ID to monitor")

    args = parser.parse_args()

    try:
        monitor = ProgressMonitor(args.work_id)
        return monitor.start()
    except RuntimeError as e:
        if "Timeout" in str(e):
            if "WorkService" in str(e):
                print("\n❌ Erro: Sistema WorkService não foi inicializado a tempo.", file=sys.stderr)
            else:
                print(f"\n❌ Erro de execução: {e}", file=sys.stderr)
        else:
            print(f"\n❌ Erro de execução: {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        if "not found" in str(e):
            print(f"\n🔍 Erro: Work '{args.work_id}' não encontrado.", file=sys.stderr)
        elif "Timeout" in str(e):
            print(f"\n⏰ Erro: Timeout aguardando inicialização do work '{args.work_id}'.", file=sys.stderr)
        else:
            print(f"\n⚠️ Erro de configuração: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\n💥 Erro inesperado: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
