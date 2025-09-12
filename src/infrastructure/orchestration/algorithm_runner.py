"""AlgorithmRunner: executes a single unit (repetition/trial/sample).
Assumes return (center, distance, metadata) for CSPAlgorithm.
"""

from __future__ import annotations

import threading
import time
import traceback
from typing import Any, Dict, Optional

from src.domain.algorithms import AlgorithmResult, CSPAlgorithm, global_registry
from src.domain.distance import DistanceCalculator
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.infrastructure.monitoring.persistence_monitor import PersistenceMonitor

# Create module logger
logger = get_logger("CSPBench.AlgorithmRunner")


class AlgorithmExecutionError(Exception):
    """Exception raised during algorithm execution."""
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
    """Execute a CSP algorithm with timeout control and persistent monitoring.

    Args:
        algorithm_name: Name of the registered algorithm.
        strings: List of input strings.
        alphabet: Alphabet used.
        distance_calculator: Injected distance calculator.
        execution_controller: Controller with limits and timeouts.
        monitor: Execution scope persistence for recording events.
        seed: Optional seed.
        params: Additional algorithm parameters.

    Returns:
        Dict with keys:
            status: BaseStatus (RUNNING/COMPLETED/CANCELED/ERROR/SKIPPED)
            algorithm_result: AlgorithmResult | None
            error: str | None
            duration_s: float
            timeout_triggered: bool
            actual_params: dict
    """
    if params is None:
        params = {}

    # Check first if not paused or canceled
    initial_status = execution_controller.check_status()
    if initial_status in [BaseStatus.PAUSED, BaseStatus.CANCELED]:
        logger.info(f"Algorithm '{algorithm_name}' will not be executed - current status: {initial_status.value}")
        return {
            "status": initial_status,
            "algorithm_result": None,
            "error": f"Execution not started - work is {initial_status.value.lower()}",
            "duration_s": 0.0,
            "execution_time_s": 0.0,
            "timeout_triggered": False,
            "actual_params": params,
        }

    logger.info(f"Starting algorithm execution: {algorithm_name}")
    logger.debug(
        f"Input parameters: strings={len(strings)}, alphabet='{alphabet}', params={params}"
    )

    t0 = time.time()
    # Used for propagation in external exception block
    error_message: str | None = None
    try:
        # Registry lookup with detailed logging
        logger.debug(f"Looking for algorithm '{algorithm_name}' in global registry")
        cls = global_registry.get(algorithm_name)
        if cls is None:
            error_msg = f"Algorithm not registered: {algorithm_name}"
            logger.error(error_msg)
            logger.debug(f"Available algorithms: {list(global_registry.keys())}")
            raise AlgorithmExecutionError(error_msg)

        logger.info(f"Algorithm '{algorithm_name}' found: {cls.__name__}")

        # Algorithm instantiation with detailed logging
        logger.debug(f"Creating algorithm instance with {len(strings)} strings")

        # Remove 'seed' from params to avoid duplicate keyword argument
        filtered_params = {k: v for k, v in params.items() if k != "seed"}

        alg: CSPAlgorithm = cls(
            strings=strings,
            alphabet=alphabet,
            distance_calculator=distance_calculator,
            monitor=monitor,
            seed=seed,
            internal_jobs=execution_controller.internal_jobs,
            **filtered_params,
        )
        logger.info(f"Algorithm instance '{algorithm_name}' created successfully")

        # Check status again before starting execution
        pre_execution_status = execution_controller.check_status()
        if pre_execution_status in [BaseStatus.PAUSED, BaseStatus.CANCELED]:
            duration = time.time() - t0
            logger.info(f"Algorithm '{algorithm_name}' will not be executed - status changed to: {pre_execution_status.value}")
            return {
                "status": pre_execution_status,
                "algorithm_result": None,
                "error": f"Execution interrupted - work was {pre_execution_status.value.lower()}",
                "duration_s": duration,
                "execution_time_s": 0.0,
                "timeout_triggered": False,
                "actual_params": params,
            }

        # Algorithm execution with timing and timeout control
        logger.info(f"Starting algorithm execution '{algorithm_name}'")
        start_time = time.time()

        # Determine effective timeout (priority params.max_time > controller.timeout_per_item)
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
                error_message = f"Error in algorithm execution: {run_exc}"
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

                # Wait loop with timeout and cancellation checking
                while thread.is_alive():
                    thread.join(timeout=1)
                    elapsed = time.time() - start_time
                    # Check external cancellation
                    status = execution_controller.check_status()
                    if status == BaseStatus.CANCELED or status == BaseStatus.PAUSED:
                        if monitor:
                            monitor.on_warning("Execution canceled externally")
                        timeout_triggered = False
                        break
                    if elapsed > max_time:
                        timeout_triggered = True
                        if monitor:
                            monitor.on_warning(
                                f"Timeout reached ({max_time}s) - requesting cancellation"
                            )
                            monitor.cancel()
                        break
        except TimeoutError as timeout_exc:
            timeout_triggered = True
            error_message = str(timeout_exc)
            if monitor:
                monitor.on_warning(f"ExecutionController timeout: {timeout_exc}")
            logger.warning(f"Timeout applied by ExecutionController: {timeout_exc}")

        # If still alive after timeout/cancel, there's no safe way to kill pure Python thread.
        # We log state as timeout; algorithm should check monitor.is_cancelled periodically to finish.
        if thread.is_alive():
            logger.warning(
                f"Algorithm thread '{algorithm_name}' still active after timeout/cancel - marking state as unfinished"
            )

        execution_time = time.time() - start_time
        total_duration = time.time() - t0

        # Determine final status
        if timeout_triggered:
            final_status = BaseStatus.CANCELED
            if algorithm_result is None:
                error_message = error_message or "Timeout reached"
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

        # Build simplified return
        result_dict: Dict[str, Any] = {
            "status": final_status,
            "algorithm_result": algorithm_result,
            "error": error_message,
            "duration_s": total_duration,
            "execution_time_s": execution_time,
            "timeout_triggered": timeout_triggered,
            "actual_params": (
                algorithm_result.get("parameters", getattr(alg, "params", params))
                if algorithm_result
                else getattr(alg, "params", params)
            ),
        }

        # Friendly final log
        if algorithm_result and algorithm_result.get("success"):
            logger.info(
                f"Algorithm '{algorithm_name}' completed: dist={algorithm_result['max_distance']}, exec_time={execution_time:.3f}s"
            )
        elif timeout_triggered:
            logger.warning(
                f"Algorithm '{algorithm_name}' finished by timeout ({max_time}s) in {execution_time:.3f}s"
            )
        elif final_status == BaseStatus.CANCELED:
            logger.warning(
                f"Algorithm '{algorithm_name}' canceled externally after {execution_time:.3f}s"
            )
        else:
            logger.error(
                f"Algorithm '{algorithm_name}' finished with error in {execution_time:.3f}s: {error_message}"
            )

        return result_dict

    except Exception as e:  # noqa: BLE001
        duration = time.time() - t0
        try:
            monitor.on_error(f"Unexpected error in algorithm execution: {e}", e)
        except Exception:
            pass

        logger.error(
            f"Error executing algorithm '{algorithm_name}': {e} ({e.__class__.__name__})"
        )
        logger.debug("Full traceback:", exc_info=True)

        return {
            "status": BaseStatus.FAILED,
            "algorithm_result": None,
            "error": error_message or str(e),
            "duration_s": duration,
            "timeout_triggered": False,
            "actual_params": params,
            "traceback": traceback.format_exc(limit=5),
        }
