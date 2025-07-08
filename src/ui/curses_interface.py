"""
Interface curses para acompanhar execução de algoritmos CSP.

Classes:
    CursesInterface: Interface principal para acompanhar execução
    ExecutionStatus: Status de execução de algoritmos
    ProgressTracker: Rastreador de progresso visual

Funções:
    run_with_curses_interface: Executa algoritmos com interface curses
"""

import curses
import logging
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import TYPE_CHECKING, Any, Dict, Optional

from src.core.interfaces.executor import TaskStatus

if TYPE_CHECKING:
    from .curses_integration import AlgorithmExecutionTracker

logger = logging.getLogger(__name__)


class WorkerStatus(Enum):
    """Status dos workers de execução."""

    IDLE = "idle"
    RUNNING = "running"
    FINISHED = "finished"
    ERROR = "error"
    FAILED = "failed"
    COMPLETED = "completed"
    TIMEOUT = "timeout"


class ExecutionStatus(Enum):
    """Status de execução de algoritmos."""

    PENDING = "pendente"
    RUNNING = "executando"
    COMPLETED = "concluído"
    ERROR = "erro"
    TIMEOUT = "timeout"


class AlgorithmProgress:
    """Progresso de execução de um algoritmo."""

    def __init__(self, name: str, total_runs: int):
        self.name = name
        self.total_runs = total_runs
        self.current_run = 0
        self.status = ExecutionStatus.PENDING
        self.start_time: float | None = None
        self.end_time: float | None = None
        self.best_distance = float("inf")
        self.current_distance = float("inf")
        self.error_message = ""
        self.iterations = 0
        self.progress_message = ""

    def start(self):
        """Inicia a execução."""
        self.status = ExecutionStatus.RUNNING
        self.start_time = time.time()

    def update(
        self,
        run_index: int,
        distance: float | None = None,
        iterations: int | None = None,
        message: str = "",
    ):
        """Atualiza progresso."""
        self.current_run = run_index + 1
        if distance is not None:
            self.current_distance = distance
            if distance < self.best_distance:
                self.best_distance = distance
        if iterations is not None:
            self.iterations = iterations
        self.progress_message = message

    def complete(self, success: bool = True, error_message: str = ""):
        """Completa a execução."""
        self.status = ExecutionStatus.COMPLETED if success else ExecutionStatus.ERROR
        self.end_time = time.time()
        self.error_message = error_message

    def get_elapsed_time(self) -> float:
        """Retorna tempo decorrido."""
        if self.start_time is None:
            return 0.0
        end_time = self.end_time or time.time()
        return end_time - self.start_time

    def get_progress_percent(self) -> float:
        """Retorna porcentagem de progresso."""
        if self.total_runs == 0:
            return 0.0
        return (self.current_run / self.total_runs) * 100


class BatchProgress:
    """Progresso de execução de um batch."""

    def __init__(self, name: str):
        self.name = name
        self.total_configs = 0
        self.current_config = 0
        self.config_name = ""
        self.current_base = 0
        self.total_bases = 0
        self.distancia_string_base: int | None = None
        self.dataset_info = {}
        self.algorithms: dict[str, AlgorithmProgress] = {}
        self.start_time: float | None = None
        self.end_time: float | None = None
        self.status = ExecutionStatus.PENDING

    def start(self):
        """Inicia o batch."""
        self.status = ExecutionStatus.RUNNING
        self.start_time = time.time()

    def set_current_config(
        self,
        config_index: int,
        config_name: str | None = None,
        total_configs: int | None = None,
    ):
        """Define configuração atual."""
        self.current_config = config_index
        if config_name:
            self.config_name = config_name
        if total_configs:
            self.total_configs = total_configs
        elif self.total_configs is None:
            self.total_configs = 1

    def set_current_base(
        self,
        base_index: int,
        total_bases: int | None = None,
        distancia_string_base: Optional[int] = None,
        dataset_info: Optional[Dict[str, Any]] = None,
    ):
        """Define base atual."""
        self.current_base = base_index
        if total_bases:
            self.total_bases = total_bases
        elif self.total_bases is None:
            self.total_bases = 1

        if distancia_string_base is not None:
            self.distancia_string_base = distancia_string_base

        if dataset_info is not None:
            self.dataset_info = dataset_info

    def add_algorithm(self, name: str, total_runs: int):
        """Adiciona algoritmo ao rastreamento."""
        self.algorithms[name] = AlgorithmProgress(name, total_runs)

    def get_algorithm(self, name: str) -> AlgorithmProgress | None:
        """Obtém progresso de algoritmo."""
        return self.algorithms.get(name)

    def complete(self, success: bool = True):
        """Completa o batch."""
        self.status = ExecutionStatus.COMPLETED if success else ExecutionStatus.ERROR
        self.end_time = time.time()

    def get_elapsed_time(self) -> float:
        """Retorna tempo decorrido."""
        if self.start_time is None:
            return 0.0
        end_time = self.end_time or time.time()
        return end_time - self.start_time

    def get_overall_progress(self) -> float:
        """Retorna progresso geral."""
        if self.total_configs == 0:
            return 0.0
        progress = 0.0
        current = 0.0
        if self.total_configs > 0:
            current = (self.current_config - 1) / self.total_configs
        elif self.total_bases > 0:
            current = (self.current_base - 1) / self.total_bases

        progress = current * 100 if current > 0 else 0

        return progress


class CursesInterface:
    """Interface curses para acompanhar execução."""

    def __init__(self, stdscr=None):
        self.stdscr = stdscr
        self.height = 24  # Valor padrão
        self.width = 80  # Valor padrão
        if stdscr:
            self.height, self.width = stdscr.getmaxyx()

        self.batch_progress: BatchProgress | None = None
        self.lock = threading.Lock()
        self.running = True

        # Cores serão configuradas em start()
        self.colors = {}

        # Sistema de slots para compatibilidade com CursesExecutionMonitor
        self.max_workers = 8
        self.algorithm_slots: dict[int, AlgorithmSlot] = {}
        self.algorithm_trackers: Optional[Dict[str, "AlgorithmExecutionTracker"]] = (
            None  # Para armazenar referência aos trackers do CursesExecutionMonitor
        )

        # Estado do countdown automático
        self.countdown_active = False
        self.countdown_message = ""

        # Inicializar slots
        for i in range(self.max_workers):
            self.algorithm_slots[i] = AlgorithmSlot(slot_id=i)

        # Thread para loop principal (compatibilidade)
        self.main_thread: Optional[threading.Thread] = None

    def set_batch_progress(self, batch_progress: BatchProgress):
        """Define progresso do batch."""
        with self.lock:
            self.batch_progress = batch_progress

    def draw_header(self, y: int) -> int:
        """Desenha cabeçalho."""
        if not self.stdscr or not self.colors:
            return y + 2

        title = "CSP-BLFGA - Execução de Algoritmos"
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Centralizar título
        title_x = max(0, (self.width - len(title)) // 2)
        cyan_attr = self.colors.get("cyan", curses.A_BOLD)
        self.stdscr.addstr(y, title_x, title, cyan_attr | curses.A_BOLD)

        # Timestamp no canto direito
        timestamp_x = max(0, self.width - len(timestamp) - 1)
        blue_attr = self.colors.get("blue", curses.A_NORMAL)
        self.stdscr.addstr(y, timestamp_x, timestamp, blue_attr)

        return y + 2

    def draw_batch_info(self, y: int) -> int:
        """Desenha informações do batch."""
        if not self.batch_progress or not self.stdscr:
            return y

        batch = self.batch_progress

        # Título do batch
        self.stdscr.addstr(
            y,
            0,
            f"Batch: {batch.name}",
            self.colors.get("yellow", curses.A_BOLD) | curses.A_BOLD,
        )
        y += 1

        # Progresso geral
        progress = batch.get_overall_progress()
        elapsed = batch.get_elapsed_time()

        self.stdscr.addstr(y, 0, f"Progresso Geral: {progress:.1f}%")
        self.stdscr.addstr(y, 25, f"Tempo Decorrido: {elapsed:.1f}s")
        y += 1

        # Configuração atual
        self.stdscr.addstr(
            y, 0, f"Config: {batch.current_config}/{batch.total_configs}"
        )
        if batch.total_bases > 0:
            self.stdscr.addstr(y, 20, f"Base: {batch.current_base}/{batch.total_bases}")
        y += 2

        # Informações do dataset
        if batch.dataset_info and self.stdscr:
            info = batch.dataset_info
            dataset_str = f"Dataset: {info.get('type', 'N/A')}"
            if "n" in info and "L" in info:
                dataset_str += f" | Seqs: {info['n']} | Tamanho: {info['L']}"
            if "alphabet" in info:
                dataset_str += f" | Alfabeto: {info['alphabet']}"

            green_attr = self.colors.get("green", curses.A_NORMAL)
            self.stdscr.addstr(y, 0, dataset_str[: self.width - 1], green_attr)
            y += 1

        # Distância da string base (se disponível)
        if self.batch_progress.distancia_string_base is not None:
            dist_base = self.batch_progress.distancia_string_base
            if dist_base is not None and dist_base != "-":
                base_str = f"Distância String Base: {dist_base}"
                cyan_attr = self.colors.get("cyan", curses.A_NORMAL)
                self.stdscr.addstr(y, 0, base_str[: self.width - 1], cyan_attr)
                y += 1

        return y + 1

    def draw_algorithms(self, y: int) -> int:
        """Desenha progresso dos algoritmos."""
        if not self.batch_progress or not self.stdscr:
            return y

        batch = self.batch_progress

        # Cabeçalho
        yellow_attr = self.colors.get("yellow", curses.A_BOLD)
        self.stdscr.addstr(y, 0, "Algoritmos:", yellow_attr | curses.A_BOLD)
        y += 1

        # Linha de separação
        self.stdscr.addstr(y, 0, "-" * min(80, self.width - 1))
        y += 1

        # Mostrar algoritmos do batch se existirem
        if batch.algorithms:
            for alg_name, alg_progress in batch.algorithms.items():
                if y >= self.height - 2:
                    break

                # Nome do algoritmo
                color = self.colors.get("green", curses.A_NORMAL)
                if alg_progress.status == ExecutionStatus.RUNNING:
                    color = self.colors.get("yellow", curses.A_BOLD)
                elif alg_progress.status == ExecutionStatus.ERROR:
                    color = self.colors.get("red", curses.A_REVERSE)
                elif alg_progress.status == ExecutionStatus.COMPLETED:
                    color = self.colors.get("green", curses.A_NORMAL)

                self.stdscr.addstr(y, 0, f"{alg_name:<15}", color | curses.A_BOLD)

                # Status
                status_str = alg_progress.status.value.upper()
                self.stdscr.addstr(y, 16, f"{status_str:<12}")

                # Progresso
                progress = alg_progress.get_progress_percent()
                self.stdscr.addstr(y, 29, f"{progress:5.1f}%")

                # Execução atual
                if alg_progress.total_runs > 0:
                    self.stdscr.addstr(
                        y, 36, f"({alg_progress.current_run}/{alg_progress.total_runs})"
                    )

                # Melhor distância
                if alg_progress.best_distance != float("inf"):
                    dist_str = f"{alg_progress.best_distance:.4f}"
                    self.stdscr.addstr(y, 45, f"Dist: {dist_str}")

                # Tempo decorrido
                elapsed = alg_progress.get_elapsed_time()
                if elapsed > 0:
                    time_x = min(60, self.width - 10)
                    self.stdscr.addstr(y, time_x, f"{elapsed:.1f}s")

                y += 1

                # Exibir execuções ativas agrupadas (múltiplas por linha)
                if (
                    alg_progress.status == ExecutionStatus.RUNNING
                    and self.algorithm_trackers
                    and alg_name in self.algorithm_trackers
                    and y < self.height - 2
                ):

                    tracker = self.algorithm_trackers[alg_name]
                    active_executions = tracker.get_active_executions()

                    if active_executions:
                        cyan_attr = self.colors.get("cyan", curses.A_BOLD)
                        for exec_info in active_executions:
                            if y >= self.height - 2:
                                break
                            exec_msg = f"→ Executando ({exec_info.run_index + 1}/{tracker.total_runs})"

                            # Adicionar mensagem de progresso se disponível
                            if exec_info.progress_message:
                                exec_msg += f" - {exec_info.progress_message}"

                            # Truncar mensagem se muito longa
                            max_msg_len = self.width - 5
                            if len(exec_msg) > max_msg_len:
                                exec_msg = exec_msg[: max_msg_len - 3] + "..."

                            self.stdscr.addstr(y, 2, exec_msg, cyan_attr)
                            y += 1

                # Exibir mensagem de progresso apenas se não houver algorithm_trackers (fallback)
                elif (
                    alg_progress.progress_message
                    and alg_progress.status == ExecutionStatus.RUNNING
                    and not self.algorithm_trackers
                    and y < self.height - 2
                ):
                    msg = alg_progress.progress_message[: self.width - 5]
                    cyan_attr = self.colors.get("cyan", curses.A_BOLD)
                    self.stdscr.addstr(y, 2, f"→ {msg}", cyan_attr)
                    y += 1

        # Mostrar informações dos slots ativos se não houver algoritmos no batch
        elif any(
            slot.status != WorkerStatus.IDLE for slot in self.algorithm_slots.values()
        ):
            y += 1
            cyan_attr = self.colors.get("cyan", curses.A_BOLD)
            self.stdscr.addstr(y, 0, "Workers Ativos:", cyan_attr | curses.A_BOLD)
            y += 1

            for slot in self.algorithm_slots.values():
                if slot.status == WorkerStatus.IDLE or y >= self.height - 2:
                    continue

                # Mapear status para cor
                color = self.colors.get("green", curses.A_NORMAL)
                if slot.status == WorkerStatus.RUNNING:
                    color = self.colors.get("yellow", curses.A_BOLD)
                elif slot.status == WorkerStatus.FAILED:
                    color = self.colors.get("red", curses.A_REVERSE)
                elif slot.status == WorkerStatus.COMPLETED:
                    color = self.colors.get("green", curses.A_NORMAL)

                # Formato: [slot] algoritmo status tempo
                slot_info = f"[{slot.slot_id:2d}] {slot.algorithm_name:<15} {slot.status.value:<10}"
                if slot.elapsed_time > 0:
                    slot_info += f" {slot.elapsed_time:6.1f}s"

                self.stdscr.addstr(y, 2, slot_info, color)
                y += 1

        return y

    def draw_progress_bar(
        self,
        y: int,
        x: int,
        width: int,
        progress: float,
        filled_char: str = "█",
        empty_char: str = "░",
    ) -> None:
        """Desenha barra de progresso."""
        if width <= 0 or not self.stdscr:
            return

        filled_width = int((progress / 100.0) * width)
        filled_width = max(0, min(filled_width, width))

        # Barra preenchida
        if filled_width > 0:
            green_attr = self.colors.get("green", curses.A_NORMAL)
            self.stdscr.addstr(y, x, filled_char * filled_width, green_attr)

        # Barra vazia
        if filled_width < width:
            self.stdscr.addstr(y, x + filled_width, empty_char * (width - filled_width))

    def draw_footer(self, y: int) -> int:
        """Desenha rodapé."""
        if not self.stdscr:
            return y + 1

        # Usar mensagem de countdown se estiver ativa, caso contrário usar rodapé normal
        if self.countdown_active and self.countdown_message:
            footer = self.countdown_message
            footer_attr = self.colors.get("yellow", curses.A_BOLD) | curses.A_BOLD
        else:
            footer = "Pressione 'q' para sair | 'r' para atualizar"
            footer_attr = self.colors.get("blue", curses.A_NORMAL)

        footer_x = max(0, (self.width - len(footer)) // 2)

        if y < self.height - 1:
            self.stdscr.addstr(y, footer_x, footer, footer_attr)

        return y + 1

    def refresh_display(self):
        """Atualiza a tela."""
        if not self.stdscr:
            return

        with self.lock:
            try:
                self.stdscr.clear()

                y = 0
                y = self.draw_header(y)
                y = self.draw_batch_info(y)
                y = self.draw_algorithms(y)

                # Rodapé
                if y < self.height - 2:
                    y = self.height - 2
                self.draw_footer(y)

                self.stdscr.refresh()

            except curses.error:
                # Ignorar erros de desenho (terminal muito pequeno, etc.)
                pass

    def run(self):
        """Loop principal da interface."""
        if not self.stdscr:
            return

        # Configurar input não-bloqueante
        self.stdscr.nodelay(True)
        self.stdscr.timeout(100)  # 100ms timeout

        while self.running:
            try:
                # Atualizar tela
                self.refresh_display()

                # Verificar input do usuário
                key = self.stdscr.getch()
                if key == ord("q") or key == ord("Q"):
                    self.running = False
                elif key == ord("r") or key == ord("R"):
                    # Forçar atualização
                    pass

            except KeyboardInterrupt:
                self.running = False
                break
            except Exception as e:
                logger.error("Erro na interface curses: %s", e)
                break

        # Limpar curses
        curses.curs_set(1)

    def start(self, stdscr):
        """
        Inicia a interface curses (compatibilidade com CursesExecutionMonitor).

        Args:
            stdscr: Tela curses padrão
        """
        self.stdscr = stdscr
        self.height, self.width = stdscr.getmaxyx()
        self.running = True

        # Configurar curses
        curses.curs_set(0)  # Esconder cursor
        stdscr.nodelay(True)  # Não bloquear em getch()
        stdscr.timeout(100)  # Timeout de 100ms

        # Inicializar cores
        self.colors = {}
        if curses.has_colors():
            curses.start_color()
            curses.use_default_colors()

            # Definir pares de cores
            curses.init_pair(1, curses.COLOR_GREEN, -1)  # Verde
            curses.init_pair(2, curses.COLOR_YELLOW, -1)  # Amarelo
            curses.init_pair(3, curses.COLOR_RED, -1)  # Vermelho
            curses.init_pair(4, curses.COLOR_BLUE, -1)  # Azul
            curses.init_pair(5, curses.COLOR_CYAN, -1)  # Ciano
            curses.init_pair(6, curses.COLOR_MAGENTA, -1)  # Magenta

            self.colors = {
                "green": curses.color_pair(1),
                "yellow": curses.color_pair(2),
                "red": curses.color_pair(3),
                "blue": curses.color_pair(4),
                "cyan": curses.color_pair(5),
                "magenta": curses.color_pair(6),
            }
        else:
            # Terminal não suporta cores, usar atributos padrão
            self.colors = {
                "green": curses.A_NORMAL,
                "yellow": curses.A_BOLD,
                "red": curses.A_REVERSE,
                "blue": curses.A_UNDERLINE,
                "cyan": curses.A_BOLD,
                "magenta": curses.A_STANDOUT,
            }

        # Criar um batch progress básico apenas se não houver um já configurado
        # Não sobrescrever informações reais do dataset
        if not self.batch_progress:
            self.batch_progress = BatchProgress("Execução de Algoritmos")
            self.batch_progress.start()
            self.batch_progress.set_current_config(0, "Execução Principal", 1)

        # Iniciar thread principal
        self.main_thread = threading.Thread(target=self._main_loop, daemon=True)
        self.main_thread.start()

    def stop(self):
        """Para a interface curses (compatibilidade com CursesExecutionMonitor)."""
        self.running = False
        if self.main_thread:
            self.main_thread.join(timeout=1.0)

    def update_worker(
        self,
        slot_id: int,
        algorithm_name: str,
        status: WorkerStatus,
        task_id: Optional[str] = None,
        result: Optional[Any] = None,
    ):
        """
        Atualiza um worker específico (compatibilidade com CursesExecutionMonitor).

        Args:
            slot_id: ID do slot do worker
            algorithm_name: Nome do algoritmo
            status: Status atual
            task_id: ID da tarefa (opcional)
            result: Resultado da tarefa (opcional)
        """
        with self.lock:
            if slot_id not in self.algorithm_slots:
                return

            slot = self.algorithm_slots[slot_id]

            # Atualizar dados do slot
            slot.algorithm_name = algorithm_name
            slot.status = status
            slot.task_id = task_id
            slot.result = result

            # Gerenciar tempo de início
            if status == WorkerStatus.RUNNING and slot.start_time is None:
                slot.start_time = time.time()
            elif status in [
                WorkerStatus.COMPLETED,
                WorkerStatus.FAILED,
                WorkerStatus.TIMEOUT,
            ]:
                if slot.start_time:
                    slot.elapsed_time = time.time() - slot.start_time
                slot.start_time = None
            elif status == WorkerStatus.IDLE:
                slot.start_time = None
                slot.elapsed_time = 0.0

            # Sincronizar com sistema de algoritmos do batch
            if self.batch_progress:
                alg_progress = self.batch_progress.get_algorithm(algorithm_name)
                if not alg_progress:
                    # Criar novo progresso de algoritmo
                    self.batch_progress.add_algorithm(algorithm_name, 1)
                    alg_progress = self.batch_progress.get_algorithm(algorithm_name)

                if alg_progress:
                    # Mapear status
                    if status == WorkerStatus.RUNNING:
                        alg_progress.start()
                        alg_progress.progress_message = f"Executando no slot {slot_id}"
                    elif status == WorkerStatus.COMPLETED:
                        if result and hasattr(result, "result_data"):
                            # Extrair distância do resultado se disponível
                            result_data = result.result_data
                            if isinstance(result_data, dict):
                                distance = result_data.get(
                                    "distance",
                                    result_data.get("distancia", float("inf")),
                                )
                                alg_progress.best_distance = distance
                        alg_progress.complete(success=True)
                    elif status in [WorkerStatus.FAILED, WorkerStatus.TIMEOUT]:
                        error_msg = (
                            "Timeout"
                            if status == WorkerStatus.TIMEOUT
                            else "Falha na execução"
                        )
                        if result and hasattr(result, "error_message"):
                            error_msg = result.error_message
                        alg_progress.complete(success=False, error_message=error_msg)

    def _main_loop(self):
        """Loop principal da interface (compatibilidade)."""
        while self.running:
            try:
                # Atualizar tempo decorrido dos slots em execução
                with self.lock:
                    current_time = time.time()
                    for slot in self.algorithm_slots.values():
                        if slot.status == WorkerStatus.RUNNING and slot.start_time:
                            slot.elapsed_time = current_time - slot.start_time

                # Atualizar tela
                self.refresh_display()

                # Verificar input do usuário
                if self.stdscr:
                    key = self.stdscr.getch()
                    if key == ord("q") or key == ord("Q"):
                        self.running = False
                        break
                    elif key == ord("r") or key == ord("R"):
                        # Forçar atualização
                        pass

                time.sleep(0.1)

            except KeyboardInterrupt:
                self.running = False
                break
            except Exception as e:
                logger.error("Erro no loop principal da interface curses: %s", e)
                time.sleep(0.1)

    def _safe_addstr(self, y: int, x: int, text: str, attr=0):
        """Adiciona string de forma segura."""
        if not self.stdscr or y >= self.height or x >= self.width:
            return
        try:
            # Truncar texto se necessário
            max_len = self.width - x - 1
            if len(text) > max_len:
                text = text[:max_len]
            self.stdscr.addstr(y, x, text, attr)
        except curses.error:
            pass

    def set_countdown_message(self, message: str):
        """Define mensagem de countdown para exibir no rodapé."""
        with self.lock:
            self.countdown_active = bool(message)
            self.countdown_message = message

    def clear_countdown_message(self):
        """Limpa mensagem de countdown."""
        with self.lock:
            self.countdown_active = False
            self.countdown_message = ""


def create_curses_console_adapter(
    interface: CursesInterface, batch_progress: BatchProgress
):
    """Cria adaptador para console que atualiza interface curses."""

    class CursesConsoleAdapter:
        def __init__(self, interface: CursesInterface, batch_progress: BatchProgress):
            self.interface = interface
            self.batch_progress = batch_progress

        def print(self, message: str, end: str = "\n", flush: bool = False):
            """Processa mensagens do console."""
            # Extrair informações da mensagem para atualizar interface
            if "Executando" in message and "execução" in message:
                # Exemplo: "Executando BLF-GA (execução 2/5)..."
                parts = message.split()
                if len(parts) >= 2:
                    alg_name = parts[1]
                    alg_progress = self.batch_progress.get_algorithm(alg_name)
                    if alg_progress:
                        alg_progress.progress_message = message.strip()

        def print_inline(self, message: str, flush: bool = True):
            """Processa mensagens inline."""
            pass

        def print_warning(self, prefix: str, message: str):
            """Processa warnings."""
            pass

    return CursesConsoleAdapter(interface, batch_progress)


def run_batch_with_curses(batch_executor: Any, execucoes: list[Any]):
    """Executa batch com interface curses."""

    def run_curses_interface(stdscr):
        # Criar interface
        interface = CursesInterface(stdscr)

        # Criar progresso do batch
        batch_progress = BatchProgress("Execução Batch")
        batch_progress.start()
        interface.set_batch_progress(batch_progress)

        # Thread para execução do batch
        batch_thread = threading.Thread(
            target=execute_batch_with_progress,
            args=(batch_executor, execucoes, batch_progress, interface),
        )
        batch_thread.daemon = True
        batch_thread.start()

        # Executar interface
        interface.run()

        # Aguardar conclusão
        batch_thread.join(timeout=1.0)

    return curses.wrapper(run_curses_interface)


def execute_batch_with_progress(
    batch_executor: Any,
    execucoes: list[Any],
    batch_progress: BatchProgress,
    interface: CursesInterface,
):
    """Executa batch atualizando progresso."""

    try:
        # Executar configurações
        for config_index, config in enumerate(execucoes):
            batch_progress.set_current_config(
                config_index, config.nome, config.num_bases
            )

            # Adicionar algoritmos ao rastreamento
            for alg_name in config.algoritmos:
                batch_progress.add_algorithm(
                    alg_name, config.execucoes_por_algoritmo_por_base
                )

            # Executar bases
            for base_index in range(config.num_bases):
                # Gerar dataset
                seqs, dataset_params = batch_executor._generate_dataset(
                    config.dataset_config
                )

                # Atualizar informações do dataset
                dataset_info = {
                    "type": config.dataset_config.get("tipo", "unknown"),
                    "n": len(seqs),
                    "L": len(seqs[0]) if seqs else 0,
                    "alphabet": "".join(sorted(set("".join(seqs)))) if seqs else "",
                }
                batch_progress.set_current_base(base_index, dataset_info=dataset_info)

                # Executar algoritmos
                for alg_name in config.algoritmos:
                    alg_progress = batch_progress.get_algorithm(alg_name)
                    if alg_progress:
                        alg_progress.start()

                        try:
                            # Executar algoritmo com callback de progresso
                            executions = execute_algorithm_with_progress(
                                alg_name, seqs, config, alg_progress
                            )

                            # Atualizar com melhor resultado
                            if executions:
                                valid_results = [
                                    e
                                    for e in executions
                                    if e.get("distancia", float("inf")) != float("inf")
                                ]
                                if valid_results:
                                    best_exec = min(
                                        valid_results, key=lambda e: e["distancia"]
                                    )
                                    alg_progress.best_distance = best_exec["distancia"]
                                    alg_progress.complete(success=True)
                                else:
                                    alg_progress.complete(
                                        success=False,
                                        error_message="Nenhum resultado válido",
                                    )
                            else:
                                alg_progress.complete(
                                    success=False, error_message="Execução falhou"
                                )

                        except Exception as e:
                            alg_progress.complete(success=False, error_message=str(e))

        batch_progress.complete(success=True)

    except Exception as e:
        batch_progress.complete(success=False)
        logger.error("Erro na execução do batch: %s", e)


def execute_algorithm_with_progress(
    alg_name: str, seqs: list[str], config: Any, alg_progress: AlgorithmProgress
):
    """Executa algoritmo com atualização de progresso."""

    from algorithms.base import global_registry

    if alg_name not in global_registry:
        raise ValueError(f"Algoritmo não encontrado: {alg_name}")

    AlgClass = global_registry[alg_name]
    alphabet = "".join(sorted(set("".join(seqs))))

    # Criar callbacks de progresso
    def progress_callback(msg: str):
        alg_progress.progress_message = msg

    def iteration_callback(iteration: int, distance: float):
        alg_progress.update(alg_progress.current_run, distance, iteration)

    # Executar runs
    executions = []
    total_runs = getattr(config, "execucoes_por_algoritmo_por_base", 1)
    timeout = getattr(config, "timeout", 300)

    for run_index in range(total_runs):
        alg_progress.update(run_index)

        # Executar run individual
        try:
            from src.core.interfaces import create_executor

            # Criar executor com timeout
            executor = create_executor(
                executor_type="scheduler", timeout_seconds=timeout
            )

            try:
                instance = AlgClass(seqs, alphabet)

                # Configurar callbacks se suportados
                if hasattr(instance, "set_progress_callback"):
                    instance.set_progress_callback(progress_callback)
                if hasattr(instance, "set_iteration_callback"):
                    instance.set_iteration_callback(iteration_callback)

                # Submeter tarefa e aguardar resultado
                handle = executor.submit(instance)

                # Aguardar conclusão
                while executor.poll(handle) == TaskStatus.RUNNING:
                    import time

                    time.sleep(0.1)

                # Obter resultado
                result = executor.result(handle)

                if isinstance(result, Exception):
                    raise result

                if result and hasattr(result, "center") and hasattr(result, "distance"):
                    # Resultado do novo sistema
                    executions.append(
                        {
                            "tempo": result["metadata"].get("execution_time", 0.0),
                            "distancia": result["distance"],
                            "melhor_string": (
                                result["center"] if result["center"] else ""
                            ),
                            "iteracoes": result["metadata"].get("iteracoes", 0),
                        }
                    )
                else:
                    executions.append(
                        {
                            "tempo": 0.0,
                            "distancia": float("inf"),
                            "erro": "Resultado vazio",
                        }
                    )
            finally:
                # Garantir que o executor seja encerrado
                if hasattr(executor, "shutdown"):
                    executor.shutdown(wait=True)

        except Exception as e:
            logger.error("Erro executando %s run %s: %s", alg_name, run_index, e)
            executions.append({"tempo": 0.0, "distancia": float("inf"), "erro": str(e)})

    return executions


# Mapeamento de slots para algoritmos
@dataclass
class AlgorithmSlot:
    """Representa um slot de algoritmo para execução."""

    slot_id: int
    algorithm_name: str = "---"
    status: WorkerStatus = WorkerStatus.IDLE
    task_id: Optional[str] = None
    start_time: Optional[float] = None
    elapsed_time: float = 0.0
    result: Optional[Any] = None
