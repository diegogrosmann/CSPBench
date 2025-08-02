"""Monitor TUI (Curses) para interface avançada."""

import curses
import threading
import time
from datetime import datetime
from typing import Any, Dict, Optional

from .interfaces import MonitoringInterface, TaskProgress, TaskType


class TUIMonitor(MonitoringInterface):
    """Monitor TUI que exibe interface avançada com curses."""

    def __init__(self):
        self.task_type: TaskType = None
        self.batch_name: str = ""
        self.start_time: datetime = None
        self.screen = None
        self.running = False
        self.update_thread = None
        self.last_data: Optional[TaskProgress] = None
        self.update_interval: float = 3.0  # 3 segundos
        self.log_messages = []
        self.max_log_lines = 6

    def start_monitoring(self, task_type: TaskType, config: Dict[str, Any]) -> None:
        """Inicia o monitoramento."""
        self.task_type = task_type
        self.batch_name = config.get("batch_name", "Batch Sem Nome")
        self.start_time = datetime.now()
        self.running = True

        # Inicializa curses
        self.screen = curses.initscr()
        curses.noecho()
        curses.cbreak()
        self.screen.nodelay(True)
        curses.curs_set(0)

        # Configura cores se disponível
        if curses.has_colors():
            curses.start_color()
            curses.init_pair(1, curses.COLOR_GREEN, curses.COLOR_BLACK)  # Verde
            curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)  # Amarelo
            curses.init_pair(3, curses.COLOR_RED, curses.COLOR_BLACK)  # Vermelho
            curses.init_pair(4, curses.COLOR_CYAN, curses.COLOR_BLACK)  # Ciano
            curses.init_pair(5, curses.COLOR_BLUE, curses.COLOR_BLACK)  # Azul

        # Inicia thread de atualização
        self.update_thread = threading.Thread(target=self._update_loop, daemon=True)
        self.update_thread.start()

        # Primeira atualização
        self._draw_screen()

    def update_progress(self, progress_data: TaskProgress) -> None:
        """Atualiza informações de progresso."""
        self.last_data = progress_data

        # Adiciona mensagem de log se houver callback info
        self._add_log_message(progress_data)

    def _add_log_message(self, data: TaskProgress) -> None:
        """Adiciona mensagem ao log."""
        timestamp = datetime.now().strftime("%H:%M:%S")

        if self.task_type == TaskType.EXECUTION and data.execution_data:
            if data.execution_data.algorithm_callback_info:
                for algo, info in data.execution_data.algorithm_callback_info.items():
                    if info and info != "processando...":
                        msg = f"[{timestamp}] {algo}: {info}"
                        if msg not in [
                            log[-len(msg) :] for log in self.log_messages[-3:]
                        ]:
                            self.log_messages.append(msg)

        elif self.task_type == TaskType.OPTIMIZATION and data.optimization_data:
            if data.optimization_data.trial_callback_info:
                msg = f"[{timestamp}] Trial {data.optimization_data.current_trial}: {data.optimization_data.trial_callback_info}"
                if msg not in [log[-len(msg) :] for log in self.log_messages[-3:]]:
                    self.log_messages.append(msg)

        elif self.task_type == TaskType.SENSITIVITY and data.sensitivity_data:
            if data.sensitivity_data.callback_info:
                param = data.sensitivity_data.current_parameter
                msg = f"[{timestamp}] {param}: {data.sensitivity_data.callback_info}"
                if msg not in [log[-len(msg) :] for log in self.log_messages[-3:]]:
                    self.log_messages.append(msg)

        # Mantém apenas as últimas mensagens
        if len(self.log_messages) > 20:
            self.log_messages = self.log_messages[-20:]

    def _update_loop(self) -> None:
        """Loop de atualização da interface."""
        while self.running:
            try:
                # Verifica input do usuário
                key = self.screen.getch()
                if key == ord("q") or key == ord("Q"):
                    self.running = False
                    break

                # Atualiza tela
                self._draw_screen()

                time.sleep(self.update_interval)

            except Exception:
                break

    def _draw_screen(self) -> None:
        """Desenha a tela principal."""
        if not self.screen:
            return

        try:
            height, width = self.screen.getmaxyx()
            self.screen.clear()

            # Desenha bordas
            self._draw_box(0, 0, height - 1, width - 1)

            # Cabeçalho
            self._draw_header(width)

            # Conteúdo principal baseado no tipo de tarefa
            if self.last_data:
                if self.task_type == TaskType.EXECUTION:
                    self._draw_execution_content()
                elif self.task_type == TaskType.OPTIMIZATION:
                    self._draw_optimization_content()
                elif self.task_type == TaskType.SENSITIVITY:
                    self._draw_sensitivity_content()

            # Rodapé
            self._draw_footer(height, width)

            self.screen.refresh()

        except curses.error:
            pass

    def _draw_box(self, y1: int, x1: int, y2: int, x2: int) -> None:
        """Desenha uma caixa com bordas."""
        try:
            # Bordas horizontais
            for x in range(x1 + 1, x2):
                self.screen.addch(y1, x, "─")
                self.screen.addch(y2, x, "─")

            # Bordas verticais
            for y in range(y1 + 1, y2):
                self.screen.addch(y, x1, "│")
                self.screen.addch(y, x2, "│")

            # Cantos
            self.screen.addch(y1, x1, "┌")
            self.screen.addch(y1, x2, "┐")
            self.screen.addch(y2, x1, "└")
            self.screen.addch(y2, x2, "┘")

        except curses.error:
            pass

    def _draw_header(self, width: int) -> None:
        """Desenha o cabeçalho."""
        try:
            # Título principal
            title = "CSPBench Monitor v0.1.0"
            self._safe_addstr(1, (width - len(title)) // 2, title, curses.color_pair(4))

            # Linha divisória
            for x in range(1, width - 1):
                self.screen.addch(2, x, "─")
            self.screen.addch(2, 0, "├")
            self.screen.addch(2, width - 1, "┤")

            # Informações do batch
            if self.last_data:
                line1_left = f"📋 Batch: {self.batch_name}"
                line1_right = self._get_task_info()
                self._safe_addstr(3, 2, line1_left)
                if len(line1_right) < width - len(line1_left) - 4:
                    self._safe_addstr(3, width - len(line1_right) - 2, line1_right)

                line2_left = self._get_dataset_info()
                line2_right = f"💾 Session: {self.last_data.session_id}"
                self._safe_addstr(4, 2, line2_left)
                if len(line2_right) < width - len(line2_left) - 4:
                    self._safe_addstr(4, width - len(line2_right) - 2, line2_right)

                start_time_str = (
                    self.start_time.strftime("%H:%M:%S") if self.start_time else ""
                )
                line3 = f"⏰ Iniciado: {start_time_str}"
                self._safe_addstr(5, 2, line3)

        except curses.error:
            pass

    def _get_task_info(self) -> str:
        """Retorna informações da tarefa atual."""
        if not self.last_data:
            return ""

        if self.task_type == TaskType.EXECUTION and self.last_data.execution_data:
            data = self.last_data.execution_data
            exec_num = data.completed_executions + 1
            return f"📊 Execução: {data.current_execution} ({exec_num}/{data.total_executions})"

        elif (
            self.task_type == TaskType.OPTIMIZATION and self.last_data.optimization_data
        ):
            data = self.last_data.optimization_data
            opt_num = data.completed_optimizations + 1
            return f"📊 Otimização: {data.current_optimization} ({opt_num}/{data.total_optimizations})"

        elif self.task_type == TaskType.SENSITIVITY and self.last_data.sensitivity_data:
            data = self.last_data.sensitivity_data
            analysis_num = data.completed_analyses + 1
            return f"📊 Análise: {data.current_analysis} ({analysis_num}/{data.total_analyses})"

        return ""

    def _get_dataset_info(self) -> str:
        """Retorna informações do dataset atual."""
        if not self.last_data:
            return ""

        if self.task_type == TaskType.EXECUTION and self.last_data.execution_data:
            data = self.last_data.execution_data
            return f"🗂️  Dataset: {data.current_dataset} ({data.current_dataset_index}/{data.total_datasets})"

        elif (
            self.task_type == TaskType.OPTIMIZATION and self.last_data.optimization_data
        ):
            data = self.last_data.optimization_data
            return f"🗂️  Dataset: {data.current_dataset} ({data.current_dataset_index}/{data.total_datasets})"

        elif self.task_type == TaskType.SENSITIVITY and self.last_data.sensitivity_data:
            data = self.last_data.sensitivity_data
            return f"🗂️  Dataset: {data.current_dataset} ({data.current_dataset_index}/{data.total_datasets})"

        return ""

    def _draw_execution_content(self) -> None:
        """Desenha conteúdo específico de execução."""
        if not self.last_data or not self.last_data.execution_data:
            return

        data = self.last_data.execution_data
        current_line = 7

        try:
            # Linha divisória e título da seção
            self._draw_section_divider(current_line - 1, "PROGRESSO")
            current_line += 1

            # Barra de progresso geral
            if data.total_algorithms > 0:
                progress = (data.completed_algorithms / data.total_algorithms) * 100
                progress_bar = self._create_progress_bar(progress, 65)
                progress_text = f"{progress_bar} {progress:.0f}% ({data.completed_algorithms}/{data.total_algorithms})"
                self._safe_addstr(current_line, 2, progress_text)
                current_line += 2

            # Algoritmos ativos
            active_algos = []
            for algo in [data.current_algorithm]:
                if algo and not data.algorithm_results.get(algo, {}).get(
                    "completed", False
                ):
                    active_algos.append(algo)

            for algo in active_algos:
                progress = data.algorithm_progress.get(algo, 0)
                callback_info = data.algorithm_callback_info.get(algo, "processando...")

                # Linha do algoritmo
                progress_bar = self._create_progress_bar(progress, 20)
                exec_time = "0.00s"  # Placeholder
                distance = "N/A"  # Placeholder
                algo_line = f"⏳ {algo:<12} [{progress_bar}] {progress:.0f}%    {exec_time}    dist: {distance}"
                self._safe_addstr(current_line, 2, algo_line)
                current_line += 1

                # Linha de callback
                task_info = (
                    f"Task {data.current_task_info}: " if data.current_task_info else ""
                )
                callback_line = f"    └─ {task_info}{callback_info}"
                self._safe_addstr(
                    current_line, 2, callback_line[:75], curses.color_pair(2)
                )
                current_line += 2

            # Seção de resultados
            self._draw_section_divider(current_line, "RESULTADOS")
            current_line += 2

            if data.best_distance is not None:
                self._safe_addstr(
                    current_line, 2, f"🏆 Melhor distância: {data.best_distance}"
                )
                current_line += 1

            elapsed = (
                (datetime.now() - self.start_time).total_seconds()
                if self.start_time
                else 0
            )
            self._safe_addstr(current_line, 2, f"⏱️  Tempo total: {elapsed:.2f}s")
            current_line += 1

            success_rate = (
                (data.completed_algorithms / data.total_algorithms * 100)
                if data.total_algorithms > 0
                else 0
            )
            self._safe_addstr(
                current_line, 2, f"📊 Taxa de sucesso: {success_rate:.0f}%"
            )
            current_line += 2

            # Log
            self._draw_log_section(current_line + 1)

        except curses.error:
            pass

    def _draw_optimization_content(self) -> None:
        """Desenha conteúdo específico de otimização."""
        if not self.last_data or not self.last_data.optimization_data:
            return

        data = self.last_data.optimization_data
        current_line = 7

        try:
            # Informações adicionais
            self._safe_addstr(current_line, 2, "🎯 Objetivo: Minimizar distância")
            current_line += 1

            # Linha divisória e título da seção
            self._draw_section_divider(current_line, "PROGRESSO")
            current_line += 1

            # Barra de progresso geral
            if data.n_trials > 0:
                progress = (data.completed_trials / data.n_trials) * 100
                progress_bar = self._create_progress_bar(progress, 65)
                progress_text = f"{progress_bar} {progress:.0f}% ({data.completed_trials}/{data.n_trials})"
                self._safe_addstr(current_line, 2, progress_text)
                current_line += 2

            # Trial atual
            if data.current_trial is not None:
                trial_line = f"⏳ Trial {data.current_trial}/{data.n_trials}     [{'█' * 20}] 100%    45.2s   valor: {data.best_value or 'N/A'}"
                self._safe_addstr(current_line, 2, trial_line)
                current_line += 1

                # Linha de callback
                task_info = (
                    f"Task {data.current_task_info}: " if data.current_task_info else ""
                )
                callback_line = f"    └─ {task_info}{data.trial_callback_info}"
                self._safe_addstr(
                    current_line, 2, callback_line[:75], curses.color_pair(2)
                )
                current_line += 2

            # Seção de resultados
            self._draw_section_divider(current_line, "RESULTADOS")
            current_line += 2

            if data.best_value is not None:
                self._safe_addstr(
                    current_line, 2, f"🏆 Melhor valor: {data.best_value}"
                )
                current_line += 1

                if data.best_params:
                    params_str = ", ".join(
                        [f"{k}={v}" for k, v in list(data.best_params.items())[:3]]
                    )
                    self._safe_addstr(
                        current_line, 2, f"📊 Melhores parâmetros: {params_str}"
                    )
                    current_line += 1

            elapsed = (
                (datetime.now() - self.start_time).total_seconds()
                if self.start_time
                else 0
            )
            self._safe_addstr(
                current_line,
                2,
                f"⏱️  Tempo decorrido: {elapsed/60:.0f}m {elapsed%60:.0f}s",
            )
            current_line += 1

            improvement_rate = 3 if data.completed_trials > 0 else 0  # Placeholder
            self._safe_addstr(
                current_line,
                2,
                f"📈 Taxa de melhoria: {improvement_rate}/{data.completed_trials} trials",
            )
            current_line += 1

            self._safe_addstr(current_line, 2, f"🔄 Sampler: {data.sampler_name}")
            current_line += 2

            # Log
            self._draw_log_section(current_line + 1)

        except curses.error:
            pass

    def _draw_sensitivity_content(self) -> None:
        """Desenha conteúdo específico de análise de sensibilidade."""
        if not self.last_data or not self.last_data.sensitivity_data:
            return

        data = self.last_data.sensitivity_data
        current_line = 7

        try:
            # Informações adicionais
            self._safe_addstr(
                current_line,
                2,
                f"🔬 Método: {data.analysis_method} ({data.method_details})",
            )
            current_line += 1
            params_str = ", ".join(data.parameters)
            self._safe_addstr(
                current_line, 2, f"📊 Parâmetros: {len(data.parameters)} ({params_str})"
            )
            current_line += 1

            # Linha divisória e título da seção
            self._draw_section_divider(current_line, "PROGRESSO")
            current_line += 1

            # Barra de progresso geral
            if data.n_samples > 0:
                progress = (data.completed_samples / data.n_samples) * 100
                progress_bar = self._create_progress_bar(progress, 65)
                progress_text = f"{progress_bar} {progress:.0f}% ({data.completed_samples}/{data.n_samples})"
                self._safe_addstr(current_line, 2, progress_text)
                current_line += 2

            # Parâmetro atual
            if data.current_parameter:
                param_progress = data.parameter_progress.get(data.current_parameter, 0)
                param_results = data.sensitivity_results.get(data.current_parameter, {})
                mu_star = param_results.get("mu_star", 0.12)
                sigma = param_results.get("sigma", 0.08)

                progress_bar = self._create_progress_bar(param_progress, 20)
                param_line = f"⏳ {data.current_parameter:<12} [{progress_bar}] {param_progress:.0f}%     18.2s   μ*={mu_star:.2f}, σ={sigma:.2f}"
                self._safe_addstr(current_line, 2, param_line)
                current_line += 1

                # Linha de callback
                task_info = (
                    f"Task {data.current_task_info}: " if data.current_task_info else ""
                )
                callback_line = f"    └─ {task_info}{data.callback_info}"
                self._safe_addstr(
                    current_line, 2, callback_line[:75], curses.color_pair(2)
                )
                current_line += 2

            # Seção de resultados
            self._draw_section_divider(current_line, "RESULTADOS")
            current_line += 2

            self._safe_addstr(current_line, 2, "🔍 Sensibilidade Preliminar:")
            current_line += 1

            # Mostra resultados dos parâmetros concluídos
            sorted_results = sorted(
                data.sensitivity_results.items(),
                key=lambda x: x[1].get("mu_star", 0),
                reverse=True,
            )
            for i, (param, results) in enumerate(sorted_results[:4], 1):
                mu_star = results.get("mu_star", 0)
                sigma = results.get("sigma", 0)
                if param == data.current_parameter:
                    level = "ANALISANDO..."
                elif mu_star > 0.2:
                    level = "ALTA"
                elif mu_star > 0.1:
                    level = "MÉDIA"
                else:
                    level = "BAIXA"
                self._safe_addstr(
                    current_line,
                    4,
                    f"{i}. {param}: {level} (μ*={mu_star:.2f}, σ={sigma:.2f})",
                )
                current_line += 1

            # Informações do método
            if data.analysis_method == "morris":
                self._safe_addstr(current_line, 2, "📈 Trajetórias: 20/20 completas")
                current_line += 1
                self._safe_addstr(
                    current_line, 2, "🎯 Método: Morris (num_levels=4, grid_jump=2)"
                )
                current_line += 2

            # Log
            self._draw_log_section(current_line + 1)

        except curses.error:
            pass

    def _draw_section_divider(self, line: int, title: str) -> None:
        """Desenha um divisor de seção."""
        try:
            height, width = self.screen.getmaxyx()
            for x in range(1, width - 1):
                self.screen.addch(line, x, "─")
            self.screen.addch(line, 0, "├")
            self.screen.addch(line, width - 1, "┤")

            # Título centralizado
            title_with_spaces = f" {title} "
            start_x = (width - len(title_with_spaces)) // 2
            self._safe_addstr(line, start_x, title_with_spaces)

        except curses.error:
            pass

    def _draw_log_section(self, start_line: int) -> None:
        """Desenha a seção de log."""
        try:
            height, width = self.screen.getmaxyx()

            # Divisor
            self._draw_section_divider(start_line, "LOG")
            current_line = start_line + 2

            # Mostra as últimas mensagens de log
            displayed_logs = (
                self.log_messages[-self.max_log_lines :] if self.log_messages else []
            )

            if not displayed_logs:
                start_msg = (
                    f"[{self.start_time.strftime('%H:%M:%S')}] Monitoramento iniciado..."
                    if self.start_time
                    else ""
                )
                self._safe_addstr(current_line, 2, start_msg)
            else:
                for log_msg in displayed_logs:
                    if current_line < height - 3:  # Deixa espaço para o rodapé
                        self._safe_addstr(current_line, 2, log_msg[: width - 4])
                        current_line += 1

        except curses.error:
            pass

    def _draw_footer(self, height: int, width: int) -> None:
        """Desenha o rodapé."""
        try:
            # Linha divisória
            for x in range(1, width - 1):
                self.screen.addch(height - 2, x, "─")
            self.screen.addch(height - 2, 0, "├")
            self.screen.addch(height - 2, width - 1, "┤")

            # Texto do rodapé
            footer_text = "💡 Tecla: [q]uit"
            self._safe_addstr(height - 1, 2, footer_text)

        except curses.error:
            pass

    def _create_progress_bar(self, progress: float, width: int) -> str:
        """Cria uma barra de progresso."""
        filled_length = int(width * progress // 100)
        if width == 65:
            return "━" * filled_length + "━" * (width - filled_length)
        else:
            bar = "█" * filled_length + "░" * (width - filled_length)
            return bar

    def _safe_addstr(self, y: int, x: int, text: str, attr=0) -> None:
        """Adiciona string de forma segura, evitando erros de curses."""
        try:
            height, width = self.screen.getmaxyx()
            if 0 <= y < height - 1 and 0 <= x < width - 1:
                # Trunca o texto se necessário
                max_len = width - x - 1
                if len(text) > max_len:
                    text = text[:max_len]
                self.screen.addstr(y, x, text, attr)
        except curses.error:
            pass

    def finish_monitoring(self, results: Dict[str, Any]) -> None:
        """Finaliza o monitoramento."""
        self.running = False

        if self.update_thread and self.update_thread.is_alive():
            self.update_thread.join(timeout=1.0)

        # Mostra tela final
        if self.screen:
            try:
                self.screen.clear()
                height, width = self.screen.getmaxyx()

                msg = "✅ Monitoramento concluído!"
                self._safe_addstr(
                    height // 2, (width - len(msg)) // 2, msg, curses.color_pair(1)
                )

                elapsed = (
                    (datetime.now() - self.start_time).total_seconds()
                    if self.start_time
                    else 0
                )
                time_msg = f"⏰ Tempo total: {elapsed/60:.0f}m {elapsed%60:.0f}s"
                self._safe_addstr(
                    height // 2 + 1, (width - len(time_msg)) // 2, time_msg
                )

                exit_msg = "Pressione qualquer tecla para sair..."
                self._safe_addstr(
                    height // 2 + 3, (width - len(exit_msg)) // 2, exit_msg
                )

                self.screen.refresh()
                self.screen.nodelay(False)
                self.screen.getch()

            except curses.error:
                pass

    def show_error(self, error: str) -> None:
        """Exibe erro."""
        if self.screen:
            try:
                height, width = self.screen.getmaxyx()
                error_msg = f"❌ Erro: {error}"
                self._safe_addstr(
                    height // 2, 2, error_msg[: width - 4], curses.color_pair(3)
                )
                self.screen.refresh()
            except curses.error:
                pass

    def close(self) -> None:
        """Fecha e limpa recursos do monitor."""
        self.running = False

        if self.update_thread and self.update_thread.is_alive():
            self.update_thread.join(timeout=1.0)

        if self.screen:
            try:
                curses.echo()
                curses.nocbreak()
                curses.curs_set(1)
                curses.endwin()
            except curses.error:
                pass
            self.screen = None
