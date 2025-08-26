"""AlgorithmRunner: executa uma unidade (repetition/trial/sample).
Assume retorno (center, distance, metadata) para CSPAlgorithm.
"""

from __future__ import annotations
from typing import Any, Dict, Optional
import time
import traceback
import threading
from src.domain.algorithms import global_registry, CSPAlgorithm, AlgorithmResult
from src.domain.distance import DistanceCalculator
from src.infrastructure.logging_config import get_logger
from src.infrastructure.monitoring.persistence_monitor import PersistenceMonitor
from src.infrastructure.execution_control import ExecutionController
from src.domain.status import BaseStatus

# Create module logger
logger = get_logger("CSPBench.AlgorithmRunner")


class AlgorithmExecutionError(Exception):
    pass


def run_algorithm(
    algorithm_name: str,
    strings: list[str],
    alphabet: str,
    distance_calculator: DistanceCalculator,
    execution_controller: ExecutionController,
    monitor: PersistenceMonitor,
    seed: Optional[int] = None,
    params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Executa um algoritmo CSP com controle de timeout e monitoração persistente.

    Args:
        algorithm_name: Nome do algoritmo registrado.
        strings: Lista de strings de entrada.
        alphabet: Alfabeto utilizado.
        distance_calculator: Calculadora de distância injetada.
        execution_controller: Controller com limites e timeouts.
        exec_persistence: Persistência escopo execução para registrar eventos.
        seed: Semente opcional.
        params: Parâmetros adicionais do algoritmo.

    Returns:
        dict com chaves:
            status: BaseStatus (RUNNING/COMPLETED/CANCELED/ERROR/SKIPPED)
            algorithm_result: AlgorithmResult | None
            error: str | None
            duration_s: float
            timeout_triggered: bool
            actual_params: dict
    """
    if params is None:
        params = {}

    logger.info(f"Iniciando execução do algoritmo: {algorithm_name}")
    logger.debug(
        f"Parâmetros de entrada: strings={len(strings)}, alphabet='{alphabet}', params={params}"
    )

    t0 = time.time()
    # Usado para propagação em bloco de exceção externo
    error_message: str | None = None
    try:
        # Registry lookup with detailed logging
        logger.debug(f"Procurando algoritmo '{algorithm_name}' no registro global")
        cls = global_registry.get(algorithm_name)
        if cls is None:
            error_msg = f"Algorithm not registered: {algorithm_name}"
            logger.error(error_msg)
            logger.debug(f"Algoritmos disponíveis: {list(global_registry.keys())}")
            raise AlgorithmExecutionError(error_msg)

        logger.info(f"Algoritmo '{algorithm_name}' encontrado: {cls.__name__}")

        # Algorithm instantiation with detailed logging
        logger.debug(f"Criando instância do algoritmo com {len(strings)} strings")

        # Remove 'seed' from params to avoid duplicate keyword argument
        filtered_params = {k: v for k, v in params.items() if k != 'seed'}

        alg: CSPAlgorithm = cls(
            strings=strings,
            alphabet=alphabet,
            distance_calculator=distance_calculator,
            monitor=monitor,
            seed=seed,
            internal_jobs=execution_controller.internal_jobs,
            **filtered_params,
        )
        logger.info(f"Instância do algoritmo '{algorithm_name}' criada com sucesso")

        # Algorithm execution with timing and timeout control
        logger.info(f"Iniciando execução do algoritmo '{algorithm_name}'")
        start_time = time.time()

        # Determinar timeout efetivo (prioridade params.max_time > controller.timeout_per_item)
        max_time = float(params.get("max_time", execution_controller.timeout_per_item))
        
        # Check if we're in main thread (signals only work in main thread)
        use_signal_timeout = threading.current_thread() is threading.main_thread()
        
        algorithm_result: AlgorithmResult | None = None
        timeout_triggered = False
        error_message: str | None = None

        def _target():
            nonlocal algorithm_result, error_message
            try:
                algorithm_result = alg.run()
            except Exception as run_exc:  # noqa: BLE001
                error_message = f"Erro em execução de algoritmo: {run_exc}"
                logger.error(error_message, exc_info=True)
                if monitor:
                    monitor.on_error(str(run_exc), run_exc)

        # Apply timeout using ExecutionController only if in main thread
        try:
            if use_signal_timeout and "max_time" not in params:
                # Use signal-based timeout only in main thread
                with execution_controller.item_timeout():
                    thread = threading.Thread(
                        target=_target, name=f"Alg-{algorithm_name}", daemon=True
                    )
                    thread.start()
                    thread.join()
            else:
                # Use polling-based timeout for subprocesses or custom max_time
                thread = threading.Thread(
                    target=_target, name=f"Alg-{algorithm_name}", daemon=True
                )
                thread.start()

                # Loop de espera com checagem de timeout e cancelamento
                while thread.is_alive():
                    thread.join(timeout=1)
                    elapsed = time.time() - start_time
                    # Check external cancellation
                    status = execution_controller.check_status()
                    if status == BaseStatus.CANCELED or status == BaseStatus.PAUSED:
                        if monitor:
                            monitor.on_warning("Execução cancelada externamente")
                        timeout_triggered = False
                        break
                    if elapsed > max_time:
                        timeout_triggered = True
                        if monitor:
                            monitor.on_warning(
                                f"Timeout atingido ({max_time}s) - solicitando cancelamento"
                            )
                            monitor.cancel()
                        break
        except TimeoutError as timeout_exc:
            timeout_triggered = True
            error_message = str(timeout_exc)
            if monitor:
                monitor.on_warning(f"Timeout do ExecutionController: {timeout_exc}")
            logger.warning(f"Timeout aplicado pelo ExecutionController: {timeout_exc}")

        # Se ainda vivo após timeout/cancel, não há forma segura de matar thread Python puro.
        # Registramos estado como timeout; algoritmo deve checar monitor.is_cancelled periodicamente para finalizar.
        if thread.is_alive():
            logger.warning(
                f"Thread de algoritmo '{algorithm_name}' ainda ativa após timeout/cancel - marcando estado como não finalizado"
            )

        execution_time = time.time() - start_time
        total_duration = time.time() - t0

        # Determinar status final
        if timeout_triggered:
            final_status = BaseStatus.CANCELED
            if algorithm_result is None:
                error_message = error_message or "Timeout atingido"
        else:
            controller_status = execution_controller.check_status()
            if not controller_status == BaseStatus.RUNNING:
                final_status = controller_status
            elif algorithm_result and algorithm_result.get("success"):
                final_status = BaseStatus.COMPLETED
            elif algorithm_result and not algorithm_result.get("success"):
                final_status = BaseStatus.ERROR
            else:
                final_status = BaseStatus.FAILED

        # Construir retorno simplificado
        result_dict: Dict[str, Any] = {
            "status": final_status,
            "algorithm_result": algorithm_result,
            "error": error_message,
            "duration_s": total_duration,
            "execution_time_s": execution_time,
            "timeout_triggered": timeout_triggered,
            "actual_params": algorithm_result.get('parameters', getattr(alg, "params", params)) if algorithm_result else getattr(alg, "params", params),
        }

        # Log final amigável
        if algorithm_result and algorithm_result.get("success"):
            logger.info(
                f"Algoritmo '{algorithm_name}' concluído: dist={algorithm_result['max_distance']}, tempo_exec={execution_time:.3f}s"
            )
        elif timeout_triggered:
            logger.warning(
                f"Algoritmo '{algorithm_name}' finalizado por timeout ({max_time}s) em {execution_time:.3f}s"
            )
        elif final_status == BaseStatus.CANCELED:
            logger.warning(
                f"Algoritmo '{algorithm_name}' cancelado externamente após {execution_time:.3f}s"
            )
        else:
            logger.error(
                f"Algoritmo '{algorithm_name}' finalizado com erro em {execution_time:.3f}s: {error_message}"
            )

        return result_dict

    except Exception as e:  # noqa: BLE001
        duration = time.time() - t0
        try:
            monitor.on_error(f"Erro inesperado na execução do algoritmo: {e}", e)
        except Exception:
            pass

        logger.error(
            f"Erro na execução do algoritmo '{algorithm_name}': {e} ({e.__class__.__name__})"
        )
        logger.debug("Traceback completo:", exc_info=True)

        return {
            "status": BaseStatus.FAILED,
            "algorithm_result": None,
            "error": error_message or str(e),
            "duration_s": duration,
            "timeout_triggered": False,
            "actual_params": params,
            "traceback": traceback.format_exc(limit=5),
        }
