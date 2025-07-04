"""
Executor de algoritmos em threads isoladas com timeout e monitoramento de recursos.

Classes:
    TimeoutException: Exceção para timeout de execução.
    ResourceLimitException: Exceção para violação de recursos.
    AlgorithmExecutor: Executor de algoritmos com timeout e monitoramento.
"""

import logging
import multiprocessing
import threading
import time
from collections.abc import Callable
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

        # Configurar limites baseados no timeout e dataset
        limits = ResourceLimits(
            max_memory_mb=2048,  # 2GB máximo
            max_iterations=min(100000, timeout_seconds * 1000),  # Baseado no timeout
            check_interval=min(
                2.0, timeout_seconds / 10
            ),  # Verificar mais frequentemente em timeouts curtos
        )
        logger.info(
            f"AlgorithmExecutor criado: timeout={timeout_seconds}s, limites={limits}"
        )
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
        logger.info(
            f"Executando algoritmo {algorithm_instance.__class__.__name__} com timeout={self.timeout}s"
        )
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
                logger.info(
                    f"Iniciando processo do algoritmo {algorithm_instance.__class__.__name__}"
                )

                center, distance = algorithm_instance.run()

                # Coletar informações adicionais
                info = {"melhor_string": center}

                # Tentar obter número de iterações de diferentes atributos
                if hasattr(algorithm_instance, "geracao"):
                    info["iteracoes"] = algorithm_instance.geracao
                elif hasattr(algorithm_instance, "iterations"):
                    info["iteracoes"] = algorithm_instance.iterations
                elif hasattr(algorithm_instance, "num_iteracoes"):
                    info["iteracoes"] = algorithm_instance.num_iteracoes
                else:
                    info["iteracoes"] = 1

                logger.info(f"Algoritmo finalizado: dist={distance}, info={info}")
                result_queue.put(("success", center, distance, info))

            except Exception as e:
                error_msg = str(e)
                logger.error(
                    f"Erro na execução do algoritmo: {error_msg}", exc_info=True
                )
                result_queue.put(("error", None, float("inf"), {"erro": error_msg}))

        # Iniciar monitoramento de recursos
        self.resource_monitor.start_monitoring()

        # Iniciar processo do algoritmo
        algorithm_process = multiprocessing.Process(
            target=run_algorithm_mp, args=(result_queue, progress_queue)
        )
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
                    logger.warning(
                        "Processo do algoritmo não terminou graciosamente após timeout"
                    )

                self.resource_monitor.stop_monitoring()
                raise TimeoutException(
                    f"Algoritmo excedeu tempo limite de {self.timeout}s"
                )

            # Verificar se há resultado disponível
            try:
                result = result_queue.get_nowait()
                status, center, distance, info = result

                self.resource_monitor.stop_monitoring()
                logger.info(f"Resultado recebido do processo: status={status}")

                if status == "timeout":
                    raise TimeoutException(
                        "Algoritmo cancelado por timeout durante execução"
                    )
                elif status == "resource_limit":
                    raise ResourceLimitException(
                        "Algoritmo cancelado por limite de recursos"
                    )
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
