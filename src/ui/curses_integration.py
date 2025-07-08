"""
Integração da Interface Curses com Sistema de Execução

Este módulo integra a interface curses com o sistema de execução paralela,
fornecendo monitoramento em tempo real de algoritmos CSP agrupados por algoritmo.
"""

import curses
import time
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

from algorithms.base import global_registry
from src.core.interfaces import TaskStatus, create_executor
from src.core.interfaces.task_result import TaskResult

from .curses_interface import (
    AlgorithmProgress,
    BatchProgress,
    CursesInterface,
    ExecutionStatus,
    WorkerStatus,
)


@dataclass
class ExecutionInfo:
    """Informações sobre uma execução específica."""

    run_index: int
    status: str
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    distance: Optional[float] = None
    error: Optional[str] = None
    progress_message: Optional[str] = None


@dataclass
class AlgorithmExecutionTracker:
    """Rastreador de execuções para um algoritmo específico."""

    name: str
    total_runs: int
    completed_runs: int = 0
    failed_runs: int = 0
    running_runs: int = 0
    pending_runs: int = 0
    executions: Optional[Dict[int, ExecutionInfo]] = None
    best_distance: Optional[float] = None
    status: str = "PENDENTE"

    def __post_init__(self):
        if self.executions is None:
            self.executions = {}
            for i in range(self.total_runs):
                self.executions[i] = ExecutionInfo(run_index=i, status="PENDENTE")
            self.pending_runs = self.total_runs

    def submit_execution(self, run_index: int):
        """Submete uma execução específica para o executor."""
        if self.executions and run_index in self.executions:
            exec_info = self.executions[run_index]

            # Atualizar contadores baseado no status anterior
            if exec_info.status == "PENDENTE":
                self.pending_runs -= 1

            # Atualizar status para submetido (ainda não executando)
            exec_info.status = "SUBMETIDO"
            self.status = "SUBMETIDO"

    def start_execution(self, run_index: int):
        """Inicia uma execução específica (chamado quando realmente começa)."""
        if self.executions and run_index in self.executions:
            exec_info = self.executions[run_index]

            # Só atualizar se estava submetido
            if exec_info.status == "SUBMETIDO":
                exec_info.status = "EXECUTANDO"
                exec_info.start_time = time.time()
                self.running_runs += 1
                self.status = "EXECUTANDO"

    def complete_execution(
        self,
        run_index: int,
        distance: Optional[float] = None,
        error: Optional[str] = None,
    ):
        """Completa uma execução específica."""
        if self.executions and run_index in self.executions:
            exec_info = self.executions[run_index]
            exec_info.end_time = time.time()
            exec_info.distance = distance
            exec_info.error = error

            # Atualizar contadores
            if exec_info.status == "EXECUTANDO":
                self.running_runs -= 1

            if error:
                exec_info.status = "ERRO"
                self.failed_runs += 1
            else:
                exec_info.status = "CONCLUÍDO"
                self.completed_runs += 1
                if distance is not None and (
                    self.best_distance is None or distance < self.best_distance
                ):
                    self.best_distance = distance

            # Atualizar status geral
            if self.completed_runs + self.failed_runs == self.total_runs:
                self.status = "CONCLUÍDO" if self.completed_runs > 0 else "ERRO"

    def get_progress_percent(self) -> float:
        """Retorna a porcentagem de progresso."""
        if self.total_runs == 0:
            return 0.0
        return (self.completed_runs / self.total_runs) * 100

    def get_active_executions(self) -> List[ExecutionInfo]:
        """Retorna execuções ativas (realmente executando)."""
        if not self.executions:
            return []
        return [
            exec_info
            for exec_info in self.executions.values()
            if exec_info.status == "EXECUTANDO"  # Só as que realmente iniciaram
        ]

    def get_submitted_executions(self) -> List[ExecutionInfo]:
        """Retorna execuções submetidas (aguardando início)."""
        if not self.executions:
            return []
        return [
            exec_info
            for exec_info in self.executions.values()
            if exec_info.status == "SUBMETIDO"
        ]


class CursesExecutionMonitor:
    """
    Monitor de execução com interface curses integrada.

    Combina o sistema de execução com a CursesInterface para
    fornecer monitoramento visual em tempo real agrupado por algoritmo.
    """

    def __init__(self, max_workers: int = 8, timeout: int = 300):
        """
        Inicializa o monitor de execução.

        Args:
            max_workers: Número máximo de workers
            timeout: Timeout por tarefa em segundos
        """
        self.max_workers = max_workers
        self.timeout = timeout
        self.interface = CursesInterface(None)  # Será configurado em start()
        self.executor = None

        # Rastreamento de algoritmos
        self.algorithm_trackers: Dict[str, AlgorithmExecutionTracker] = {}

        # Informações do dataset
        self.dataset_info: Dict[str, Any] = {}

    def _safe_getch(self):
        """
        Captura input de forma segura, ignorando eventos especiais que podem causar travamento.

        Returns:
            int: Código da tecla pressionada ou -1 se timeout/erro
        """
        if not self.interface.stdscr:
            return -1

        try:
            key = self.interface.stdscr.getch()

            # Filtrar apenas teclas que nos interessam
            if key == -1:  # Timeout
                return -1
            elif key == ord("q") or key == ord("Q"):
                return key
            elif key == ord("r") or key == ord("R"):  # Refresh
                return key
            elif key == 27:  # ESC
                return key
            elif key >= 32 and key <= 126:  # Caracteres ASCII imprimíveis
                return key
            else:
                # Ignorar todos os outros eventos especiais (mouse, teclas especiais, etc.)
                return -1

        except Exception:
            # Se houver qualquer erro na captura, retornar timeout
            return -1

    def update_algorithm_display(
        self, algorithm_trackers: Dict[str, AlgorithmExecutionTracker]
    ):
        """
        Atualiza a exibição dos algoritmos agrupados usando o sistema da CursesInterface.

        Args:
            algorithm_trackers: Dicionário com rastreadores de algoritmos
        """
        if not self.interface:
            return  # Criar ou atualizar BatchProgress
        if not self.interface.batch_progress:
            # Verificar se é execução em lote
            batch_name = self.dataset_info.get("batch_name", "Execução de Algoritmos")
            total_executions = self.dataset_info.get("total_executions", 1)

            self.interface.batch_progress = BatchProgress(name=batch_name)
            self.interface.batch_progress.start()

        batch = self.interface.batch_progress

        # Atualizar configuração e base atuais
        execution_index = self.dataset_info.get("execution_index", 1)
        execution_name = self.dataset_info.get("execution_name", "runtime")
        distancia_string_base = self.dataset_info.get("distancia_string_base", None)
        total_executions = self.dataset_info.get("total_executions", 1)
        total_bases = self.dataset_info.get("total_bases", 1)
        base_index = self.dataset_info.get("base_index", 1)

        # Configurar execução atual
        batch.set_current_config(execution_index, execution_name, total_executions)

        # Usar informações reais do dataset
        dataset_info = self.dataset_info.copy()
        if not dataset_info:
            dataset_info = {"type": "runtime", "n": 0, "L": 0, "alphabet": "N/A"}

        # Configurar base atual
        batch.set_current_base(
            base_index, total_bases, distancia_string_base, dataset_info
        )

        # Atualizar algoritmos no BatchProgress
        for alg_name, tracker in algorithm_trackers.items():
            if alg_name not in batch.algorithms:
                batch.add_algorithm(alg_name, tracker.total_runs)

            alg_progress = batch.get_algorithm(alg_name)
            if alg_progress:
                # Mapear status do tracker para ExecutionStatus
                if tracker.status == "PENDENTE":
                    alg_progress.status = ExecutionStatus.PENDING
                elif tracker.status == "EXECUTANDO":
                    alg_progress.status = ExecutionStatus.RUNNING
                elif tracker.status == "CONCLUÍDO":
                    alg_progress.status = ExecutionStatus.COMPLETED
                elif tracker.status == "ERRO":
                    alg_progress.status = ExecutionStatus.ERROR

                # Atualizar progresso
                alg_progress.current_run = tracker.completed_runs
                alg_progress.total_runs = tracker.total_runs

                # Atualizar melhor distância
                if tracker.best_distance is not None:
                    alg_progress.best_distance = tracker.best_distance

                # Armazenar informações de execuções ativas no progress_message
                active_executions = tracker.get_active_executions()
                if active_executions:
                    # Criar lista de execuções ativas para exibição
                    exec_list = [
                        f"({exec_info.run_index + 1}/{tracker.total_runs})"
                        for exec_info in active_executions
                    ]
                    alg_progress.progress_message = "|".join(exec_list)
                else:
                    alg_progress.progress_message = ""

        # Armazenar referência aos trackers na interface para acesso posterior
        self.interface.algorithm_trackers = algorithm_trackers

        # Usar o sistema de refresh da CursesInterface
        self.interface.refresh_display()

    def execute_algorithms(
        self,
        algorithm_names: List[str],
        seqs: List[str],
        alphabet: str,
        console=None,
        num_execs: int = 1,
        dataset_params: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, List[TaskResult]]:
        """
        Executa algoritmos com monitoramento curses.

        Args:
            algorithm_names: Lista de nomes dos algoritmos
            seqs: Sequências de entrada
            alphabet: Alfabeto utilizado
            console: Console para fallback (não usado com curses)
            num_execs: Número de execuções por algoritmo
            dataset_params: Parâmetros do dataset para exibição

        Returns:
            Dicionário com listas de resultados por algoritmo
        """
        # Extrair informações do dataset
        self.dataset_info = self._extract_dataset_info(seqs, alphabet, dataset_params)

        def curses_main(stdscr):
            return self._execute_with_curses(
                stdscr, algorithm_names, seqs, alphabet, num_execs
            )

        return curses.wrapper(curses_main)

    def _execute_with_curses(
        self,
        stdscr,
        algorithm_names: List[str],
        seqs: List[str],
        alphabet: str,
        num_execs: int = 1,
    ) -> Dict[str, List[TaskResult]]:
        """
        Execução interna com interface curses.

        Args:
            stdscr: Tela curses
            algorithm_names: Algoritmos a executar
            seqs: Sequências
            alphabet: Alfabeto

        Returns:
            Resultados da execução
        """
        # Inicializar interface curses
        # Primeiro, configurar o BatchProgress com as informações corretas do dataset
        if not self.interface.batch_progress:
            self.interface.batch_progress = BatchProgress(name="Execução de Algoritmos")
            self.interface.batch_progress.start()

            # Usar informações reais do dataset
            dataset_info = self.dataset_info.copy()
            if not dataset_info:
                dataset_info = {"type": "runtime", "n": 0, "L": 0, "alphabet": "N/A"}

            # Definir configuração e base com as informações do dataset
            self.interface.batch_progress.set_current_config(
                dataset_info.get("execution_index", 0),
                dataset_info.get("execution_name", "runtime"),
                dataset_info.get("total_executions", 1),
            )
            self.interface.batch_progress.set_current_base(
                dataset_info.get("base_index", 0),
                dataset_info.get("total_bases", "runtime"),
                dataset_info.get("distancia_string_base", None),
                dataset_info,
            )

        # Agora iniciar a interface (que não vai sobrescrever o BatchProgress)
        self.interface.start(stdscr)

        try:
            # Preparar tarefas
            tasks = self._prepare_tasks(algorithm_names, seqs, alphabet)

            if not tasks:
                return {}

            # Executar com monitoramento (múltiplas execuções)
            results = self._execute_tasks_with_monitoring(tasks, num_execs)

            # Mostrar resultados finais
            self._show_final_results(results)

            return results

        finally:
            self.interface.stop()

    def _prepare_tasks(
        self, algorithm_names: List[str], seqs: List[str], alphabet: str
    ) -> List[Dict[str, Any]]:
        """
        Prepara tarefas para execução.

        Args:
            algorithm_names: Lista de algoritmos
            seqs: Sequências
            alphabet: Alfabeto

        Returns:
            Lista de tarefas preparadas
        """
        tasks = []

        for algorithm_name in algorithm_names:
            try:
                # Obter classe do algoritmo
                if algorithm_name not in global_registry:
                    print(f"❌ Algoritmo {algorithm_name} não encontrado no registry")
                    continue

                algorithm_class = global_registry[algorithm_name]
                print(f"✓ Preparando {algorithm_name}...")

                # Criar instância
                instance = algorithm_class(seqs, alphabet)
                print(f"✓ {algorithm_name} instanciado com sucesso")

                # Adicionar tarefa
                tasks.append(
                    {
                        "name": algorithm_name,
                        "instance": instance,
                        "strings": seqs,
                        "alphabet": alphabet,
                    }
                )

            except Exception as e:
                # Log do erro
                print(f"❌ Erro ao preparar {algorithm_name}: {e}")
                import traceback

                traceback.print_exc()
                continue

        return tasks

    def _execute_tasks_with_monitoring(
        self, tasks: List[Dict[str, Any]], num_execs: int = 1
    ) -> Dict[str, List[TaskResult]]:
        """
        Executa tarefas com monitoramento visual agrupado por algoritmo.

        Para cada algoritmo, submete todas as suas execuções em paralelo,
        aguarda todas terminarem, e só então passa para o próximo algoritmo.

        Args:
            tasks: Lista de tarefas
            num_execs: Número de execuções por algoritmo

        Returns:
            Resultados da execução (lista de resultados por algoritmo)
        """
        results = {}

        # Configurar executor
        self.executor = create_executor(
            max_workers=min(self.max_workers, num_execs),
            timeout_seconds=self.timeout,
        )

        try:
            # Inicializar rastreadores de algoritmos
            for task in tasks:
                algorithm_name = task["name"]
                self.algorithm_trackers[algorithm_name] = AlgorithmExecutionTracker(
                    name=algorithm_name, total_runs=num_execs
                )
                results[algorithm_name] = []

            # Exibição inicial
            self.update_algorithm_display(self.algorithm_trackers)

            # Para cada algoritmo, executar múltiplas vezes EM PARALELO
            for task in tasks:
                algorithm_name = task["name"]
                tracker = self.algorithm_trackers[algorithm_name]
                algorithm_results = []

                # Submeter TODAS as execuções do algoritmo atual em paralelo
                handles = []
                for exec_idx in range(num_execs):
                    # Atualizar tracker - submetendo execução
                    tracker.submit_execution(exec_idx)

                    # Criar nova instância do algoritmo para cada execução
                    algorithm_class = task["instance"].__class__
                    algorithm_instance = algorithm_class(
                        task["strings"], task["alphabet"]
                    )

                    # Configurar callback de progresso para esta execução
                    def make_progress_callback(alg_name, exec_idx):
                        def progress_callback(message):
                            if alg_name in self.algorithm_trackers:
                                tracker = self.algorithm_trackers[alg_name]
                                if (
                                    tracker.executions
                                    and exec_idx in tracker.executions
                                ):
                                    exec_info = tracker.executions[exec_idx]
                                    exec_info.progress_message = message
                                    # Atualizar display em tempo real
                                    self.update_algorithm_display(
                                        self.algorithm_trackers
                                    )

                        return progress_callback

                    # Configurar callback se o algoritmo suporta
                    if hasattr(algorithm_instance, "set_progress_callback"):
                        algorithm_instance.set_progress_callback(
                            make_progress_callback(algorithm_name, exec_idx)
                        )

                    # Submeter para execução (não aguardar)
                    handle = self.executor.submit(algorithm_instance)
                    handles.append((handle, exec_idx))

                # Atualizar display após submeter todas as execuções
                self.update_algorithm_display(self.algorithm_trackers)

                # Aguardar TODAS as execuções do algoritmo atual terminarem
                while handles:
                    # Verificar handles concluídos
                    completed_handles = []
                    for handle, exec_idx in handles:
                        status = self.executor.poll(handle)

                        # Verificar se a execução realmente iniciou
                        if status == TaskStatus.RUNNING:
                            # Marcar como executando se ainda não marcou
                            if tracker.executions and exec_idx in tracker.executions:
                                exec_info = tracker.executions[exec_idx]
                                if exec_info.status == "SUBMETIDO":
                                    tracker.start_execution(exec_idx)

                        # Verificar se concluiu (não está mais em fila ou executando)
                        if status not in (TaskStatus.QUEUED, TaskStatus.RUNNING):
                            completed_handles.append((handle, exec_idx))

                    # Processar handles concluídos
                    for handle, exec_idx in completed_handles:
                        handles.remove((handle, exec_idx))

                        # Obter resultado
                        try:
                            result = self.executor.result(handle)

                            if isinstance(result, Exception):
                                # Erro na execução
                                tracker.complete_execution(exec_idx, error=str(result))

                                task_result = TaskResult(
                                    success=False,
                                    distance=None,
                                    center=None,
                                    time=0.0,
                                    error=str(result),
                                    traceback=None,
                                    metadata={"run_index": exec_idx},
                                )
                            else:
                                # Extrair dados do resultado
                                distance = result.get(
                                    "distance", result.get("distancia", float("inf"))
                                )
                                center = result.get(
                                    "center", result.get("melhor_string", "")
                                )

                                # Extrair tempo do metadata se disponível
                                execution_time = result.get("metadata", {}).get(
                                    "execution_time", 0.0
                                )
                                if execution_time == 0.0:
                                    execution_time = result.get(
                                        "tempo", result.get("time", 0.0)
                                    )

                                # Atualizar tracker - execução concluída
                                tracker.complete_execution(exec_idx, distance=distance)

                                task_result = TaskResult(
                                    success=True,
                                    distance=distance,
                                    center=center,
                                    time=execution_time,
                                    error=None,
                                    traceback=None,
                                    metadata={
                                        "original_result": result,
                                        "run_index": exec_idx,
                                    },
                                )

                            algorithm_results.append(task_result)

                        except Exception as e:
                            # Erro ao obter resultado
                            tracker.complete_execution(exec_idx, error=str(e))

                            task_result = TaskResult(
                                success=False,
                                distance=None,
                                center=None,
                                time=0.0,
                                error=str(e),
                                traceback=None,
                                metadata={"run_index": exec_idx},
                            )
                            algorithm_results.append(task_result)

                    # Atualizar display periodicamente
                    self.update_algorithm_display(self.algorithm_trackers)

                    # Breve pausa para evitar uso excessivo de CPU
                    if handles:  # Só pausar se ainda há handles pendentes
                        time.sleep(0.1)

                # Armazenar resultados do algoritmo (ordenar por run_index)
                algorithm_results.sort(key=lambda x: x.metadata.get("run_index", 0))
                results[algorithm_name] = algorithm_results

            return results

        finally:
            # Garantir que o executor seja encerrado
            if self.executor and hasattr(self.executor, "shutdown"):
                self.executor.shutdown(wait=True)

    def _show_final_results(self, results: Dict[str, List[TaskResult]]):
        """
        Exibe resultados finais na interface e fecha automaticamente em 5 segundos.

        Args:
            results: Resultados da execução (listas de resultados por algoritmo)
        """
        # Marcar batch como concluído
        if self.interface.batch_progress:
            self.interface.batch_progress.complete(success=True)

        # Exibir resultados finais
        self.update_algorithm_display(self.algorithm_trackers)

        # Aguardar 5 segundos ou tecla do usuário para encerrar
        if self.interface.stdscr:
            # Configurar timeout para leitura de tecla
            self.interface.stdscr.nodelay(True)

            # Tempo limite de 5 segundos
            start_time = time.time()
            timeout_seconds = 5.0

            while True:
                current_time = time.time()
                elapsed = current_time - start_time
                remaining = timeout_seconds - elapsed

                # Verificar se o tempo limite foi atingido
                if elapsed >= timeout_seconds:
                    break

                # Verificar entrada do usuário
                key = self._safe_getch()
                if key == ord("q") or key == ord("Q") or key == 27:  # ESC
                    break

                # Atualizar a tela com contador regressivo
                self._update_display_with_countdown(remaining)

                # Pequena pausa
                time.sleep(0.1)

    def _update_display_with_countdown(self, remaining_seconds: float):
        """
        Atualiza a tela mostrando contador regressivo.

        Args:
            remaining_seconds: Segundos restantes até fechamento automático
        """
        # Atualizar display normal
        self.update_algorithm_display(self.algorithm_trackers)

        # Adicionar mensagem de countdown na interface
        if remaining_seconds > 0:
            countdown_msg = f"Fechando automaticamente em {remaining_seconds:.1f}s (pressione 'q' para fechar agora)"
            self.interface.set_countdown_message(countdown_msg)
        else:
            self.interface.clear_countdown_message()

    def _extract_dataset_info(
        self,
        seqs: List[str],
        alphabet: str,
        dataset_params: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Extrai informações do dataset para exibição na interface.

        Args:
            seqs: Sequências do dataset
            alphabet: Alfabeto utilizado
            dataset_params: Parâmetros retornados pelo gerador/carregador de dataset

        Returns:
            Dicionário com informações formatadas do dataset
        """
        # Informações básicas
        info = {"n": len(seqs), "L": len(seqs[0]) if seqs else 0, "alphabet": alphabet}

        # Se temos parâmetros do dataset, copiar todas as informações
        if dataset_params:
            # Copiar todos os parâmetros do dataset
            info.update(dataset_params)

            # Substituir apenas as informações básicas que calculamos
            info["n"] = len(seqs)
            info["L"] = len(seqs[0]) if seqs else 0
            info["alphabet"] = alphabet
            # Determinar tipo de dataset baseado na fonte
            dataset_source = dataset_params.get("dataset_source", "unknown")

            if dataset_source == "1":  # Sintético
                info["type"] = "synthetic"
                # Adicionar parâmetros específicos do sintético
                if "noise" in dataset_params:
                    info["noise"] = dataset_params["noise"]
                if "seed" in dataset_params:
                    info["seed"] = dataset_params["seed"]
                if "fully_random" in dataset_params:
                    info["fully_random"] = dataset_params["fully_random"]

            elif dataset_source == "2":  # Arquivo
                info["type"] = "file"
                if "filename" in dataset_params:
                    info["filename"] = dataset_params["filename"]

            elif dataset_source == "3":  # Entrez
                info["type"] = "entrez"
                if "query" in dataset_params:
                    info["query"] = dataset_params["query"]
                if "db" in dataset_params:
                    info["db"] = dataset_params["db"]
            else:
                info["type"] = "unknown"
        else:
            # Sem parâmetros, usar tipo genérico
            info["type"] = "sequences"

        return info


def execute_algorithms_with_curses_monitoring(
    algorithm_names: List[str],
    seqs: List[str],
    alphabet: str,
    max_workers: int = 8,
    timeout: int = 300,
    dataset_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, TaskResult]:
    """
    Executa algoritmos com monitoramento curses (compatibilidade - apenas 1 execução).

    Args:
        algorithm_names: Lista de algoritmos
        seqs: Sequências de entrada
        alphabet: Alfabeto
        max_workers: Número máximo de workers
        timeout: Timeout por tarefa
        dataset_params: Parâmetros do dataset para exibição

    Returns:
        Dicionário com resultados por algoritmo (apenas primeiro resultado)
    """
    monitor = CursesExecutionMonitor(max_workers=max_workers, timeout=timeout)
    multi_results = monitor.execute_algorithms(
        algorithm_names, seqs, alphabet, num_execs=1, dataset_params=dataset_params
    )

    # Converter de lista para resultado único para compatibilidade
    single_results = {}
    for alg_name, results_list in multi_results.items():
        if results_list:
            single_results[alg_name] = results_list[0]
        else:
            single_results[alg_name] = TaskResult(
                success=False,
                distance=None,
                center=None,
                time=0.0,
                error="Nenhum resultado encontrado",
                traceback=None,
                metadata={},
            )

    return single_results


if __name__ == "__main__":
    # Exemplo de uso
    results = execute_algorithms_with_curses_monitoring(
        algorithm_names=["BLF-GA", "Baseline"],
        seqs=["ATCG", "ATCC", "ATCA"],
        alphabet="ATCG",
        max_workers=4,
        timeout=60,
    )

    print("Resultados:", results)
