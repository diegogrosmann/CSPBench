"""
Executor de algoritmos com timeout e monitoramento de recursos.

Classes:
    TimeoutException: Exceção para timeout de execução.
    ResourceLimitException: Exceção para violação de recursos.
    AlgorithmExecutor: Executor de algoritmos com timeout e monitoramento.
    ParallelAlgorithmExecutor: Executor paralelo usando ProcessPoolExecutor.
"""

import logging
import multiprocessing
import threading
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, TimeoutError, as_completed
from typing import Any

from csp_blfga.utils.resource_monitor import (
    ResourceLimits,
    ResourceMonitor,
    force_garbage_collection,
)

logger = logging.getLogger(__name__)


class TimeoutException(Exception):
    """
    Exceção lançada quando um algoritmo excede o tempo limite de execução.
    """

    pass


class ResourceLimitException(Exception):
    """
    Exceção lançada quando há violação dos limites de recursos do sistema.
    """

    pass


class AlgorithmExecutor:
    """
    Executa algoritmos em threads isoladas com timeout e monitoramento de recursos.

    Args:
        timeout_seconds (int): Tempo limite de execução em segundos.

    Métodos:
        execute_with_timeout(...): Executa algoritmo com timeout e monitoração.
    """

    def __init__(self, timeout_seconds: int):
        """
        Inicializa o executor com timeout e limites de recursos.

        Args:
            timeout_seconds (int): Tempo limite em segundos.
        """
        self.timeout = timeout_seconds
        self.stop_event = threading.Event()

        # Configurar limites usando o método from_config
        limits = ResourceLimits.from_config()
        # Ajustar limites baseados no timeout
        limits.max_iterations = min(100000, timeout_seconds * 1000)
        limits.check_interval = min(2.0, timeout_seconds / 10)

        logger.info(f"AlgorithmExecutor criado: timeout={timeout_seconds}s")
        self.resource_monitor = ResourceMonitor(limits)
        self.resource_violation = False

    def execute_with_timeout(
        self,
        algorithm_instance: Any,
        progress_callback: Callable[[str], None] | None = None,
        warning_callback: Callable[[str], None] | None = None,
    ) -> tuple[Any | None, int | float, dict]:
        """
        Executa um algoritmo com limite de tempo e monitoramento de recursos.

        Args:
            algorithm_instance (Any): Instância do algoritmo a ser executado.
            progress_callback (Optional[Callable[[str], None]]): Callback para atualizações de progresso.
            warning_callback (Optional[Callable[[str], None]]): Callback para avisos.

        Returns:
            Tuple[Optional[Any], Union[int, float], dict]: Resultado da execução, incluindo centro, distância e informações adicionais.
        """
        logger.info(f"Executando algoritmo {algorithm_instance.__class__.__name__} com timeout={self.timeout}s")
        result_queue = multiprocessing.Queue()
        progress_queue = multiprocessing.Queue()
        self.stop_event.clear()
        self.resource_violation = False

        # Força garbage collection antes de começar
        force_garbage_collection()

        # Configurar monitor de recursos
        def resource_violation_handler(violation_msg: str):
            logger.warning(f"Violação de recursos detectada: {violation_msg}")
            self.resource_violation = True
            self.stop_event.set()
            if warning_callback:
                warning_callback(f"Recursos limitados: {violation_msg}")

        self.resource_monitor.set_violation_callback(resource_violation_handler)

        def run_algorithm_mp(result_queue, progress_queue):
            """Função executada no processo do algoritmo."""
            try:

                def progress(msg):
                    progress_queue.put(msg)

                if hasattr(algorithm_instance, "set_progress_callback"):
                    algorithm_instance.set_progress_callback(progress)
                logger.info(f"Iniciando processo do algoritmo {algorithm_instance.__class__.__name__}")

                # Executar algoritmo e processar resultado
                result = algorithm_instance.run()

                if len(result) == 3:
                    # Nova interface CSPAlgorithm
                    center, distance, metadata = result
                    info = metadata
                elif len(result) == 2:
                    # Interface legacy
                    center, distance = result
                    info = {"melhor_string": center}

                    # Tentar obter número de iterações de diferentes atributos
                    if hasattr(algorithm_instance, "geracao"):
                        info["iteracoes"] = algorithm_instance.geracao
                    elif hasattr(algorithm_instance, "iterations"):
                        info["iteracoes"] = algorithm_instance.iterations
                    elif hasattr(algorithm_instance, "num_iteracoes"):
                        info["iteracoes"] = algorithm_instance.num_iteracoes
                    else:
                        info["iteracoes"] = 0
                else:
                    # Formato inesperado
                    center, distance = result[0], result[1]
                    info = {"erro": "Formato de resultado inesperado"}

                logger.info(f"Algoritmo finalizado: dist={distance}, info={info}")
                result_queue.put(("success", center, distance, info))

            except Exception as e:
                error_msg = str(e)
                logger.error(f"Erro na execução do algoritmo: {error_msg}", exc_info=True)
                result_queue.put(("error", None, float("inf"), {"erro": error_msg}))

        # Iniciar monitoramento de recursos
        self.resource_monitor.start_monitoring()

        # Iniciar processo do algoritmo
        algorithm_process = multiprocessing.Process(target=run_algorithm_mp, args=(result_queue, progress_queue))
        algorithm_process.daemon = True
        algorithm_process.start()

        # Aguardar resultado ou timeout com verificações mais frequentes
        start_time = time.time()
        check_interval = 0.5  # Verificar a cada 0.5 segundos

        while algorithm_process.is_alive():
            elapsed = time.time() - start_time

            # Exibir progresso recebido do processo filho
            try:
                while True:
                    msg = progress_queue.get_nowait()
                    if progress_callback:
                        progress_callback(msg)
            except Exception:
                pass

            if elapsed >= self.timeout:
                logger.warning("Timeout atingido, sinalizando parada")
                # Timeout atingido - sinalizar parada
                self.stop_event.set()
                algorithm_process.terminate()
                algorithm_process.join(timeout=2.0)

                if algorithm_process.is_alive():
                    logger.warning("Processo do algoritmo não terminou graciosamente após timeout")

                self.resource_monitor.stop_monitoring()
                raise TimeoutException(f"Algoritmo excedeu tempo limite de {self.timeout}s")

            # Verificar se há resultado disponível
            try:
                result = result_queue.get_nowait()
                status, center, distance, info = result

                self.resource_monitor.stop_monitoring()
                logger.info(f"Resultado recebido do processo: status={status}")

                if status == "timeout":
                    raise TimeoutException("Algoritmo cancelado por timeout durante execução")
                elif status == "resource_limit":
                    raise ResourceLimitException("Algoritmo cancelado por limite de recursos")
                elif status == "error":
                    return center, distance, info
                else:
                    return center, distance, info

            except Exception:
                time.sleep(check_interval)

        # Processo terminou - obter resultado final
        self.resource_monitor.stop_monitoring()

        try:
            result = result_queue.get_nowait()
            status, center, distance, info = result
            logger.info(f"Processo terminou: status={status}")
            return center, distance, info
        except Exception:
            logger.error("Processo terminou sem resultado")
            return None, float("inf"), {"erro": "Processo terminou sem resultado"}

    def cancel(self):
        """
        Cancela a execução atual.
        """
        logger.info("Cancelando execução atual")
        self.stop_event.set()
        self.resource_monitor.stop_monitoring()
        self.resource_monitor.stop_monitoring()


def _execute_algorithm_in_process(alg_class, strings, alphabet, params, timeout):
    """
    Função auxiliar para executar algoritmo em processo separado.
    Esta função é importada e executada no processo worker.
    """
    try:
        # Criar instância do algoritmo
        instance = alg_class(strings, alphabet, **params)

        # Executar algoritmo
        result = instance.run()

        # Retornar resultado padronizado
        if len(result) == 2:
            # Algoritmo legacy - adicionar metadata vazio
            center, dist = result
            metadata = {"iteracoes": 1, "legacy": True}
            return center, dist, metadata
        else:
            # Algoritmo moderno
            return result

    except Exception as e:
        # Retornar erro estruturado
        return None, float("inf"), {"erro": str(e)}


class ParallelAlgorithmExecutor:
    """
    Executor paralelo de algoritmos usando ProcessPoolExecutor.

    Permite execução de múltiplos algoritmos em paralelo com timeout
    e controle de recursos avançado.
    """

    def __init__(self, max_workers: int | None = None, timeout_seconds: int = 300):
        """
        Inicializa o executor paralelo.

        Args:
            max_workers: Número máximo de workers (padrão: CPU count)
            timeout_seconds: Timeout por execução
        """
        self.max_workers = max_workers or multiprocessing.cpu_count()
        self.timeout_seconds = timeout_seconds
        self.executor = None

    def __enter__(self):
        """Context manager entry."""
        self.executor = ProcessPoolExecutor(max_workers=self.max_workers)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self.executor:
            self.executor.shutdown(wait=True)

    def execute_parallel(self, tasks: list[dict]) -> list[dict]:
        """
        Executa múltiplas tarefas em paralelo.

        Args:
            tasks: Lista de dicionários com:
                - alg_class: Classe do algoritmo
                - strings: Strings de entrada
                - alphabet: Alfabeto
                - params: Parâmetros
                - name: Nome da tarefa (opcional)

        Returns:
            Lista de resultados com status de cada execução
        """
        if not self.executor:
            raise RuntimeError("Executor não inicializado. Use como context manager.")

        # Submeter todas as tarefas
        future_to_task = {}
        for i, task in enumerate(tasks):
            future = self.executor.submit(
                _execute_algorithm_in_process,
                task["alg_class"],
                task["strings"],
                task["alphabet"],
                task.get("params", {}),
                self.timeout_seconds,
            )
            future_to_task[future] = {
                "index": i,
                "name": task.get("name", f"Task-{i}"),
                "alg_class": task["alg_class"],
            }

        # Coletar resultados
        results: list[dict] = [{} for _ in range(len(tasks))]

        for future in as_completed(future_to_task, timeout=self.timeout_seconds + 10):
            task_info = future_to_task[future]
            index = task_info["index"]

            try:
                center, dist, metadata = future.result(timeout=1)
                results[index] = {
                    "name": task_info["name"],
                    "success": True,
                    "center": center,
                    "distance": dist,
                    "metadata": metadata,
                    "error": None,
                }
            except TimeoutError:
                results[index] = {
                    "name": task_info["name"],
                    "success": False,
                    "center": None,
                    "distance": float("inf"),
                    "metadata": {},
                    "error": f"Timeout ({self.timeout_seconds}s)",
                }
            except Exception as e:
                results[index] = {
                    "name": task_info["name"],
                    "success": False,
                    "center": None,
                    "distance": float("inf"),
                    "metadata": {},
                    "error": str(e),
                }

        return results

    def execute_single_parallel(self, alg_class, strings, alphabet, params=None, name=None):
        """
        Executa um único algoritmo em processo separado.

        Args:
            alg_class: Classe do algoritmo
            strings: Strings de entrada
            alphabet: Alfabeto
            params: Parâmetros
            name: Nome da tarefa

        Returns:
            Resultado da execução
        """
        task = {
            "alg_class": alg_class,
            "strings": strings,
            "alphabet": alphabet,
            "params": params or {},
            "name": name or alg_class.__name__,
        }

        results = self.execute_parallel([task])
        return results[0] if results else None


class ModernParallelExecutor:
    """
    Executor paralelo moderno usando ProcessPoolExecutor.

    Compatível com CSPAlgorithm e otimizado para performance.
    """

    def __init__(self, max_workers: int | None = None, timeout: int = 300):
        """
        Inicializa o executor paralelo.

        Args:
            max_workers: Número máximo de processos (padrão: CPU count)
            timeout: Timeout por execução em segundos
        """
        self.max_workers = max_workers or multiprocessing.cpu_count()
        self.timeout = timeout

    def execute_algorithm_batch(self, tasks: list[dict]) -> list[dict]:
        """
        Executa múltiplos algoritmos em paralelo.

        Args:
            tasks: Lista de tarefas, cada uma contendo:
                - alg_class: Classe do algoritmo
                - strings: Strings de entrada
                - alphabet: Alfabeto
                - params: Parâmetros (opcional)
                - name: Nome da tarefa (opcional)

        Returns:
            Lista de resultados
        """
        logger.info(f"Executando {len(tasks)} tarefas em paralelo com {self.max_workers} workers")

        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # Submeter todas as tarefas
            future_to_task = {executor.submit(self._execute_single_task, task): task for task in tasks}

            results = []

            # Coletar resultados conforme completam
            for future in as_completed(future_to_task, timeout=self.timeout * len(tasks)):
                task = future_to_task[future]
                try:
                    result = future.result(timeout=self.timeout)
                    results.append(result)
                except TimeoutError:
                    logger.warning(f"Timeout na tarefa {task.get('name', 'unknown')}")
                    results.append(
                        {
                            "name": task.get("name", "unknown"),
                            "success": False,
                            "center": "",
                            "distance": float("inf"),
                            "metadata": {"erro": "timeout"},
                            "execution_time": self.timeout,
                            "error": "timeout",
                        }
                    )
                except Exception as e:
                    logger.error(f"Erro na tarefa {task.get('name', 'unknown')}: {e}")
                    results.append(
                        {
                            "name": task.get("name", "unknown"),
                            "success": False,
                            "center": "",
                            "distance": float("inf"),
                            "metadata": {"erro": str(e)},
                            "execution_time": 0,
                            "error": str(e),
                        }
                    )

            return results

    @staticmethod
    def _execute_single_task(task: dict) -> dict:
        """
        Executa uma única tarefa de algoritmo.

        Args:
            task: Dicionário com parâmetros da tarefa

        Returns:
            Resultado da execução
        """
        start_time = time.time()

        try:
            # Extrair parâmetros
            alg_class = task["alg_class"]
            strings = task["strings"]
            alphabet = task["alphabet"]
            params = task.get("params", {})
            name = task.get("name", alg_class.__name__)

            # Criar instância do algoritmo
            algorithm = alg_class(strings, alphabet, **params)

            # Executar algoritmo
            result = algorithm.run()

            # Processar resultado baseado na interface
            if len(result) == 3:
                # Nova interface CSPAlgorithm
                center, distance, metadata = result
            elif len(result) == 2:
                # Interface legacy
                center, distance = result
                metadata = {"iteracoes": getattr(algorithm, "geracao", 0)}
            else:
                raise ValueError(f"Formato de resultado inesperado: {len(result)} elementos")

            execution_time = time.time() - start_time

            return {
                "name": name,
                "success": True,
                "center": center,
                "distance": distance,
                "metadata": metadata,
                "execution_time": execution_time,
                "error": None,
            }

        except Exception as e:
            execution_time = time.time() - start_time
            return {
                "name": task.get("name", "unknown"),
                "success": False,
                "center": "",
                "distance": float("inf"),
                "metadata": {"erro": str(e)},
                "execution_time": execution_time,
                "error": str(e),
            }
