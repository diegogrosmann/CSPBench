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
from datetime import datetime
from enum import Enum
from typing import Any

from src.core.exec.batch_executor import BatchExecutor

logger = logging.getLogger(__name__)


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

    def update(self, run_index: int, distance: float | None = None, iterations: int | None = None, message: str = ""):
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

    def __init__(self, name: str, total_configs: int):
        self.name = name
        self.total_configs = total_configs
        self.current_config = 0
        self.current_base = 0
        self.total_bases = 0
        self.dataset_info = {}
        self.algorithms: dict[str, AlgorithmProgress] = {}
        self.start_time: float | None = None
        self.end_time: float | None = None
        self.status = ExecutionStatus.PENDING

    def start(self):
        """Inicia o batch."""
        self.status = ExecutionStatus.RUNNING
        self.start_time = time.time()

    def set_current_config(self, config_index: int, config_name: str, total_bases: int):
        """Define configuração atual."""
        self.current_config = config_index + 1
        self.total_bases = total_bases
        self.current_base = 0

    def set_current_base(self, base_index: int, dataset_info: dict):
        """Define base atual."""
        self.current_base = base_index + 1
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
        base_progress = 0.0
        if self.total_bases > 0:
            base_progress = (self.current_base / self.total_bases) * (1.0 / self.total_configs)
        config_progress = (self.current_config - 1) / self.total_configs
        return (config_progress + base_progress) * 100


class CursesInterface:
    """Interface curses para acompanhar execução."""

    def __init__(self, stdscr):
        self.stdscr = stdscr
        self.height, self.width = stdscr.getmaxyx()
        self.batch_progress: BatchProgress | None = None
        self.lock = threading.Lock()
        self.running = True

        # Configurar curses
        curses.curs_set(0)  # Esconder cursor
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

    def set_batch_progress(self, batch_progress: BatchProgress):
        """Define progresso do batch."""
        with self.lock:
            self.batch_progress = batch_progress

    def draw_header(self, y: int) -> int:
        """Desenha cabeçalho."""
        title = "CSP-BLFGA - Execução de Algoritmos"
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Centralizar título
        title_x = max(0, (self.width - len(title)) // 2)
        self.stdscr.addstr(y, title_x, title, self.colors["cyan"] | curses.A_BOLD)

        # Timestamp no canto direito
        timestamp_x = max(0, self.width - len(timestamp) - 1)
        self.stdscr.addstr(y, timestamp_x, timestamp, self.colors["blue"])

        return y + 2

    def draw_batch_info(self, y: int) -> int:
        """Desenha informações do batch."""
        if not self.batch_progress:
            return y

        batch = self.batch_progress

        # Título do batch
        self.stdscr.addstr(y, 0, f"Batch: {batch.name}", self.colors["yellow"] | curses.A_BOLD)
        y += 1

        # Progresso geral
        progress = batch.get_overall_progress()
        elapsed = batch.get_elapsed_time()

        self.stdscr.addstr(y, 0, f"Progresso Geral: {progress:.1f}%")
        self.stdscr.addstr(y, 25, f"Tempo Decorrido: {elapsed:.1f}s")
        y += 1

        # Configuração atual
        self.stdscr.addstr(y, 0, f"Config: {batch.current_config}/{batch.total_configs}")
        if batch.total_bases > 0:
            self.stdscr.addstr(y, 20, f"Base: {batch.current_base}/{batch.total_bases}")
        y += 2

        # Informações do dataset
        if batch.dataset_info:
            info = batch.dataset_info
            dataset_str = f"Dataset: {info.get('type', 'N/A')}"
            if "n" in info and "L" in info:
                dataset_str += f" | Seqs: {info['n']} | Tamanho: {info['L']}"
            if "alphabet" in info:
                dataset_str += f" | Alfabeto: {info['alphabet']}"

            self.stdscr.addstr(y, 0, dataset_str[: self.width - 1], self.colors["green"])
            y += 2

        return y

    def draw_algorithms(self, y: int) -> int:
        """Desenha progresso dos algoritmos."""
        if not self.batch_progress:
            return y

        batch = self.batch_progress

        if not batch.algorithms:
            return y

        # Cabeçalho
        self.stdscr.addstr(y, 0, "Algoritmos:", self.colors["yellow"] | curses.A_BOLD)
        y += 1

        # Linha de separação
        self.stdscr.addstr(y, 0, "-" * min(80, self.width - 1))
        y += 1

        # Algoritmos
        for alg_name, alg_progress in batch.algorithms.items():
            if y >= self.height - 2:
                break

            # Nome do algoritmo
            color = self.colors["green"]
            if alg_progress.status == ExecutionStatus.RUNNING:
                color = self.colors["yellow"]
            elif alg_progress.status == ExecutionStatus.ERROR:
                color = self.colors["red"]
            elif alg_progress.status == ExecutionStatus.COMPLETED:
                color = self.colors["green"]

            self.stdscr.addstr(y, 0, f"{alg_name:<15}", color | curses.A_BOLD)

            # Status
            status_str = alg_progress.status.value.upper()
            self.stdscr.addstr(y, 16, f"{status_str:<12}")

            # Progresso
            progress = alg_progress.get_progress_percent()
            self.stdscr.addstr(y, 29, f"{progress:5.1f}%")

            # Execução atual
            if alg_progress.total_runs > 0:
                self.stdscr.addstr(y, 36, f"({alg_progress.current_run}/{alg_progress.total_runs})")

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

            # Mensagem de progresso se houver
            if alg_progress.progress_message and y < self.height - 2:
                msg = alg_progress.progress_message[: self.width - 5]
                self.stdscr.addstr(y, 2, f"→ {msg}", self.colors["cyan"])
                y += 1

        return y

    def draw_progress_bar(
        self, y: int, x: int, width: int, progress: float, filled_char: str = "█", empty_char: str = "░"
    ) -> None:
        """Desenha barra de progresso."""
        if width <= 0:
            return

        filled_width = int((progress / 100.0) * width)
        filled_width = max(0, min(filled_width, width))

        # Barra preenchida
        if filled_width > 0:
            self.stdscr.addstr(y, x, filled_char * filled_width, self.colors["green"])

        # Barra vazia
        if filled_width < width:
            self.stdscr.addstr(y, x + filled_width, empty_char * (width - filled_width))

    def draw_footer(self, y: int) -> int:
        """Desenha rodapé."""
        footer = "Pressione 'q' para sair | 'r' para atualizar"
        footer_x = max(0, (self.width - len(footer)) // 2)

        if y < self.height - 1:
            self.stdscr.addstr(y, footer_x, footer, self.colors["blue"])

        return y + 1

    def refresh_display(self):
        """Atualiza a tela."""
        with self.lock:
            self.stdscr.clear()

            try:
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
                logger.error(f"Erro na interface curses: {e}")
                break

        # Limpar curses
        curses.curs_set(1)


def create_curses_console_adapter(interface: CursesInterface, batch_progress: BatchProgress):
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


def run_batch_with_curses(batch_executor: BatchExecutor, execucoes: list[Any]):
    """Executa batch com interface curses."""

    def run_curses_interface(stdscr):
        # Criar interface
        interface = CursesInterface(stdscr)

        # Criar progresso do batch
        batch_progress = BatchProgress("Execução Batch", len(execucoes))
        batch_progress.start()
        interface.set_batch_progress(batch_progress)

        # Thread para execução do batch
        batch_thread = threading.Thread(
            target=execute_batch_with_progress, args=(batch_executor, execucoes, batch_progress, interface)
        )
        batch_thread.daemon = True
        batch_thread.start()

        # Executar interface
        interface.run()

        # Aguardar conclusão
        batch_thread.join(timeout=1.0)

    return curses.wrapper(run_curses_interface)


def execute_batch_with_progress(
    batch_executor: BatchExecutor, execucoes: list[Any], batch_progress: BatchProgress, interface: CursesInterface
):
    """Executa batch atualizando progresso."""

    try:
        # Executar configurações
        for config_index, config in enumerate(execucoes):
            batch_progress.set_current_config(config_index, config.nome, config.num_bases)

            # Adicionar algoritmos ao rastreamento
            for alg_name in config.algoritmos:
                batch_progress.add_algorithm(alg_name, config.execucoes_por_algoritmo_por_base)

            # Executar bases
            for base_index in range(config.num_bases):
                # Gerar dataset
                seqs, dataset_params = batch_executor._generate_dataset(config.dataset_config)

                # Atualizar informações do dataset
                dataset_info = {
                    "type": config.dataset_config.get("tipo", "unknown"),
                    "n": len(seqs),
                    "L": len(seqs[0]) if seqs else 0,
                    "alphabet": "".join(sorted(set("".join(seqs)))) if seqs else "",
                }
                batch_progress.set_current_base(base_index, dataset_info)

                # Executar algoritmos
                for alg_name in config.algoritmos:
                    alg_progress = batch_progress.get_algorithm(alg_name)
                    if alg_progress:
                        alg_progress.start()

                        try:
                            # Executar algoritmo com callback de progresso
                            executions = execute_algorithm_with_progress(alg_name, seqs, config, alg_progress)

                            # Atualizar com melhor resultado
                            if executions:
                                valid_results = [
                                    e for e in executions if e.get("distancia", float("inf")) != float("inf")
                                ]
                                if valid_results:
                                    best_exec = min(valid_results, key=lambda e: e["distancia"])
                                    alg_progress.best_distance = best_exec["distancia"]
                                    alg_progress.complete(success=True)
                                else:
                                    alg_progress.complete(success=False, error_message="Nenhum resultado válido")
                            else:
                                alg_progress.complete(success=False, error_message="Execução falhou")

                        except Exception as e:
                            alg_progress.complete(success=False, error_message=str(e))

        batch_progress.complete(success=True)

    except Exception as e:
        batch_progress.complete(success=False)
        logger.error(f"Erro na execução do batch: {e}")


def execute_algorithm_with_progress(alg_name: str, seqs: list[str], config: Any, alg_progress: AlgorithmProgress):
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
    for run_index in range(config.execucoes_por_algoritmo_por_base):
        alg_progress.update(run_index)

        # Executar run individual
        try:
            from src.core.exec.algorithm_executor import AlgorithmExecutor

            executor = AlgorithmExecutor(config.timeout)
            instance = AlgClass(seqs, alphabet)

            # Configurar callbacks se suportados
            if hasattr(instance, "set_progress_callback"):
                instance.set_progress_callback(progress_callback)
            if hasattr(instance, "set_iteration_callback"):
                instance.set_iteration_callback(iteration_callback)

            result = executor.execute_with_timeout(instance)[0]
            executions.append(
                {
                    "tempo": 0.0,  # Será preenchido pela execução
                    "distancia": result.get("distancia", float("inf")) if result else float("inf"),
                    "melhor_string": result.get("melhor_string", "") if result else "",
                    "iteracoes": result.get("iteracoes", 0) if result else 0,
                }
            )

        except Exception as e:
            executions.append({"tempo": 0.0, "distancia": float("inf"), "erro": str(e)})

    return executions
