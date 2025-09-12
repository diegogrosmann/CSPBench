"""ExperimentExecutor: executes ExperimentTask.

Parallelizes repetitions via ProcessPool if ``max_workers > 1`` using ``ProcessPoolExecutor``.
Implements resource controls, resume functionality and complete data saving.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from src.application.ports.repositories import AbstractExecutionEngine
from src.domain.config import AlgParams, CSPBenchConfig, ExperimentTask
from src.domain.dataset import Dataset
from src.domain.distance import create_distance_calculator
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.infrastructure.monitoring.persistence_monitor import PersistenceMonitor
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
    ExecutionScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
    WorkScopedPersistence,
)

from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.ExperimentExecutor")


import random
import time

def _worker_exec(
    experiment_unit_id: str,
    seed: int | None,
    strings: list[str],
    alphabet: str,
    distance_method: str,
    use_cache: bool,
    params: dict[str, Any],
    work_id: str,  # Changed from db_path to work_id
    internal_jobs: int,
    algorithm_name: str,
    cpu_config: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], str]:
    """Function executed in subprocess (needs to be picklable).
    
    Args:
        experiment_unit_id: Unique identifier for this experiment unit
        seed: Random seed for this repetition
        strings: Input strings for the algorithm
        alphabet: Alphabet used
        distance_method: Distance calculation method
        use_cache: Whether to use distance cache
        params: Algorithm parameters
        work_id: Work identifier for persistence
        internal_jobs: Number of internal jobs for algorithms
        algorithm_name: Name of algorithm to execute
        cpu_config: CPU configuration for worker process
        
    Returns:
        Tuple of (execution_result_dict, experiment_unit_id)
    """
    import psutil

    from src.domain.status import BaseStatus as _BaseStatus
    from src.infrastructure.execution_control import ExecutionController
    from src.infrastructure.logging_config import get_logger
    from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
        WorkScopedPersistence,
    )

    logger = get_logger("CSPBench.WorkerProcess")

    logger.debug(f"[WORKER][{experiment_unit_id}] Starting worker_exec.")

    # Apply CPU configuration to this worker process
    if cpu_config:
        try:
            current_process = psutil.Process()

            # Apply CPU affinity if specified
            if cpu_config.get("affinity"):
                current_process.cpu_affinity(cpu_config["affinity"])
                logger.info(f"[WORKER][{experiment_unit_id}] CPU affinity set to: {cpu_config['affinity']}")

            # Apply CPU priority if max_workers is limited
            if cpu_config.get("exclusive_cores", False):
                current_process.nice(5)  # Slightly lower priority than main process
                logger.info(f"[WORKER][{experiment_unit_id}] CPU priority adjusted for exclusive cores")

        except Exception as e:
            logger.warning(f"[WORKER][{experiment_unit_id}] Cannot apply CPU configuration: {e}")

    try:
        # Use factory method to create work_scoped directly
        work_scoped = WorkScopedPersistence.submit(work_id)
        execution_store = work_scoped.for_execution(experiment_unit_id)
    except Exception as e:
        logger.error(f"[WORKER][{experiment_unit_id}] Error recreating store: {e}")
        raise

    # Controller with work_id for status checks
    dummy_controller = ExecutionController(work_id=work_id)
    dummy_controller._internal_jobs = internal_jobs  # type: ignore[attr-defined]

    # Apply CPU configuration manually if provided (since we can't pass ResourcesConfig to subprocess)
    if cpu_config:
        dummy_controller._exclusive_cores = cpu_config.get("exclusive_cores", False)
        if "max_workers" in cpu_config:
            dummy_controller._max_workers = cpu_config["max_workers"]
            dummy_controller._current_workers = cpu_config["max_workers"]

        # Try to apply CPU affinity in worker process (best effort)
        try:
            if cpu_config.get("exclusive_cores", False) and hasattr(
                dummy_controller, "apply_memory_limits"
            ):
                dummy_controller.apply_memory_limits()
        except Exception as cpu_exc:
            logger.warning(
                f"[WORKER][{experiment_unit_id}] CPU configuration partially applied: {cpu_exc}"
            )

    # Create monitor with controller for cancellation checks
    monitor = PersistenceMonitor(execution_store, execution_controller=dummy_controller)

    distance_calc = create_distance_calculator(
        distance_method=distance_method, strings=strings, use_cache=use_cache
    )
    try:
        execution_store.update_execution_status(_BaseStatus.RUNNING)
        monitor.on_progress(0.0, f"Starting algorithm instance")

        result = run_algorithm(
            algorithm_name=algorithm_name,
            strings=strings,
            alphabet=alphabet,
            distance_calculator=distance_calc,
            execution_controller=dummy_controller,
            monitor=monitor,
            seed=seed,
            params=params,
        )

        monitor.on_progress(1.0, f"[{experiment_unit_id}] Algorithm instance finished")
        status_value = result.get("status")
        if hasattr(status_value, "value"):
            status_value = status_value.value
        objective = None
        alg_result = result.get("algorithm_result", None)
        if alg_result is not None:
            max_distance = alg_result.get("max_distance", None)
            if max_distance is not None:
                objective = max_distance

        execution_store.update_execution_status(
            status=status_value,
            result=alg_result,
            params=result.get("actual_params", None),
            objective=objective,
        )
        logger.debug(f"[WORKER][{experiment_unit_id}] worker_exec finished with status: {status_value}")
        return result, experiment_unit_id
    except Exception as exc:  # pragma: no cover - defensive
        execution_store.update_execution_status("failed", result={"error": str(exc)})
        logger.error(f"[WORKER][{experiment_unit_id}] worker_exec failed: {exc}")
        return {
            "status": _BaseStatus.FAILED,
            "error": str(exc),
        }, experiment_unit_id


class ExperimentExecutor(AbstractExecutionEngine):
    """Executor for experiment tasks with parallel repetition support."""
    
    _combination_store: CombinationScopedPersistence = None
    _execution_controller: ExecutionController = None
    _batch_config: CSPBenchConfig = None

    def run(
        self,
        task: ExperimentTask,
        dataset_obj: Dataset,
        alg: AlgParams,
    ) -> BaseStatus:
        """Execute the experiment following the new step-by-step process.
        
        Args:
            task: Experiment task configuration
            dataset_obj: Dataset to run experiments on
            alg: Algorithm parameters configuration
            
        Returns:
            Final execution status
        """

        repetitions = max(1, int(getattr(task, "repetitions", 1)))
        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        dataset_id = getattr(dataset_obj, "id", None)

        global_seed = getattr(self._batch_config.system, "global_seed", None)
        if global_seed is None:
            global_seed = alg.params.get("seed", None)

        logger.info(f"Starting experiment for {alg.name} on dataset {dataset_id} with {repetitions} repetitions")

        # STEP 1: Identify all necessary executions
        logger.info("STEP 1: Identifying necessary executions")
        all_needed_executions = self._identify_needed_executions(
            task, dataset_id, alg, repetitions
        )
        
        # STEP 2: List already completed executions 
        logger.info("STEP 2: Checking already completed executions")
        completed_executions = self._get_completed_executions(all_needed_executions)
        
        # STEP 3: Determine pending executions
        logger.info("STEP 3: Determining pending executions")
        pending_executions = self._get_pending_executions(
            all_needed_executions, completed_executions
        )
        
        if not pending_executions:
            logger.info("No pending executions found - all have been completed")
            return BaseStatus.COMPLETED
            
        logger.info(f"Found {len(pending_executions)} pending executions out of {len(all_needed_executions)} total")
        
        # STEP 4: Persist pending executions in batch
        logger.info("STEP 4: Persisting pending executions in batch")
        self._persist_pending_executions(pending_executions)
        
        # STEP 5: Execute pending executions
        logger.info("STEP 5: Executing pending executions")
        results = self._execute_pending_executions(
            pending_executions, strings, alphabet, global_seed
        )

        status = self._execution_controller.check_status()
        if status != BaseStatus.RUNNING:
            return status

        # Evaluate results (considers BaseStatus enum)
        if not results:
            return BaseStatus.COMPLETED
            
        flat_results = [list(item.values())[0] for item in results]
        failed_repetitions = [
            r
            for r in flat_results
            if r.get("status")
            not in (BaseStatus.COMPLETED, BaseStatus.COMPLETED.value, "completed")
        ]
        final_result = BaseStatus.ERROR if failed_repetitions else BaseStatus.COMPLETED
        return final_result

    def _identify_needed_executions(
        self, task: ExperimentTask, dataset_id: str, alg: AlgParams, repetitions: int
    ) -> list[dict[str, Any]]:
        """STEP 1: Identify all executions that need to be performed.
        
        Args:
            task: Experiment task configuration
            dataset_id: Dataset identifier
            alg: Algorithm parameters
            repetitions: Number of repetitions needed
            
        Returns:
            List of execution data dictionaries
        """
        needed_executions = []
        
        config_id = (
            getattr(self._combination_store, "_preset_id", None) or "default"
        )
        config_id = str(config_id).replace(":", "_")
        
        for r in range(1, repetitions + 1):
            experiment_unit_id = (
                f"experiment:{task.id}:{dataset_id}:{config_id}:{alg.name}:{r}"
            )
            
            execution_data = {
                "unit_id": experiment_unit_id,
                "sequencia": r,
                "task_id": task.id,
                "dataset_id": dataset_id,
                "config_id": config_id,
                "algorithm_name": alg.name,
                "algorithm_params": alg.params
            }
            needed_executions.append(execution_data)
            
        logger.debug(f"Identified {len(needed_executions)} necessary executions")
        return needed_executions

    def _get_completed_executions(self, all_executions: list[dict[str, Any]]) -> set[str]:
        """STEP 2: List all executions that already have final status.
        
        Args:
            all_executions: List of all execution data dictionaries
            
        Returns:
            Set of unit_ids for completed executions
        """
        # Final statuses we consider as "completed"
        final_statuses = [
            BaseStatus.COMPLETED.value,
            "completed",
            BaseStatus.FAILED.value,
            "failed", 
            BaseStatus.ERROR.value,
            "error",
        ]
        
        # Fetch all executions from current combination at once
        all_existing_executions = self._combination_store.get_executions()
        
        # Filter only executions with final statuses
        completed_executions = [
            ex for ex in all_existing_executions 
            if ex.get("status") in final_statuses
        ]
        
        # Create a set with unit_ids of already completed executions
        completed_unit_ids = {ex["unit_id"] for ex in completed_executions}
        
        # Filter only those that are in our list of necessary executions
        needed_unit_ids = {ex["unit_id"] for ex in all_executions}
        relevant_completed = completed_unit_ids.intersection(needed_unit_ids)
        
        logger.debug(f"Found {len(relevant_completed)} already completed executions out of {len(completed_executions)} with final status, from {len(all_existing_executions)} total in database")
        return relevant_completed

    def _get_pending_executions(
        self, all_executions: list[dict[str, Any]], completed: set[str]
    ) -> list[dict[str, Any]]:
        """STEP 3: Compare and return only executions that were not performed.
        
        Args:
            all_executions: List of all execution data dictionaries
            completed: Set of completed unit_ids
            
        Returns:
            List of pending execution data dictionaries
        """
        pending = []
        
        for execution in all_executions:
            if execution["unit_id"] not in completed:
                pending.append(execution)
                
        logger.debug(f"Determined {len(pending)} pending executions")
        return pending

    def _persist_pending_executions(self, pending_executions: list[dict[str, Any]]) -> None:
        """STEP 4: Persist the list of pending executions in batch.
        
        Args:
            pending_executions: List of pending execution data dictionaries
        """
        if not pending_executions:
            return
            
        # Prepare data for batch insertion
        executions_to_insert = []
        for execution in pending_executions:
            execution_data = {
                "unit_id": execution["unit_id"],
                "sequencia": execution["sequencia"],
                # combination_id will be added automatically by submit_executions
            }
            executions_to_insert.append(execution_data)
            
        # Use batch insertion for better performance
        inserted_count = self._combination_store.submit_executions(executions_to_insert)
        logger.info(f"Persisted {inserted_count} pending executions in batch")

    def _execute_pending_executions(
        self, 
        pending_executions: list[dict[str, Any]], 
        strings: list[str], 
        alphabet: str,
        global_seed: int | None
    ) -> list[dict[str, Any]]:
        """STEP 5: Execute the list of pending executions (sequential or parallel).
        
        Args:
            pending_executions: List of pending execution data dictionaries
            strings: Input strings for algorithms
            alphabet: Alphabet used
            global_seed: Global seed for random number generation
            
        Returns:
            List of execution result dictionaries
        """
        if not pending_executions:
            return []
            
        results: list[dict[str, Any]] = []

        try:
            max_workers = self._execution_controller.get_worker_config()["cpu"][
                "max_workers"
            ]
        except Exception:
            max_workers = 1

        if max_workers > 1:
            # Parallel execution
            logger.info(f"Executing {len(pending_executions)} executions in parallel with {max_workers} workers")
            results = self._execute_parallel(pending_executions, strings, alphabet, global_seed)
        else:
            # Sequential execution
            logger.info(f"Executing {len(pending_executions)} executions sequentially")
            results = self._execute_sequential(pending_executions, strings, alphabet, global_seed)
            
        return results

    def _execute_parallel(
        self, 
        pending_executions: list[dict[str, Any]], 
        strings: list[str], 
        alphabet: str,
        global_seed: int | None
    ) -> list[dict[str, Any]]:
        """Execute pending executions in parallel.
        
        Args:
            pending_executions: List of pending execution data dictionaries
            strings: Input strings for algorithms
            alphabet: Alphabet used
            global_seed: Global seed for random number generation
            
        Returns:
            List of execution result dictionaries
        """
        results = []
        
        # Get CPU configuration for worker processes
        cpu_config = self._execution_controller.create_worker_config()
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        try:
            max_workers = self._execution_controller.get_worker_config()["cpu"][
                "max_workers"
            ]
        except Exception:
            max_workers = 1

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            
            for execution in pending_executions:
                if self._execution_controller.check_status() != BaseStatus.RUNNING:
                    break

                # Calculate seed for this repetition
                seed = None
                if global_seed is not None:
                    seed = global_seed
                    if self._batch_config.system.seed_increment:
                        seed += execution["sequencia"]

                future = executor.submit(
                    _worker_exec,
                    execution["unit_id"],
                    seed,
                    strings,
                    alphabet,
                    self._batch_config.system.distance_method,
                    self._batch_config.system.enable_distance_cache,
                    execution["algorithm_params"],
                    self._combination_store.work_id,
                    self._execution_controller.internal_jobs,
                    execution["algorithm_name"],
                    worker_config.get("cpu"),
                )
                futures.append((execution["sequencia"], future))

            # Collect results
            for sequencia, future in futures:
                result, unit_id = future.result()
                results.append({unit_id: result})
                
        return results

    def _execute_sequential(
        self, 
        pending_executions: list[dict[str, Any]], 
        strings: list[str], 
        alphabet: str,
        global_seed: int | None
    ) -> list[dict[str, Any]]:
        """Execute pending executions sequentially.
        
        Args:
            pending_executions: List of pending execution data dictionaries
            strings: Input strings for algorithms
            alphabet: Alphabet used
            global_seed: Global seed for random number generation
            
        Returns:
            List of execution result dictionaries
        """
        results = []
        
        # Get CPU configuration for sequential execution
        cpu_config = self._execution_controller.create_worker_config()
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        for execution in pending_executions:
            status = self._execution_controller.check_status()
            if status != BaseStatus.RUNNING:
                break

            # Calculate seed for this repetition
            seed = None
            if global_seed is not None:
                seed = global_seed
                if self._batch_config.system.seed_increment:
                    seed += execution["sequencia"]

            result, _ = _worker_exec(
                execution["unit_id"],
                seed,
                strings,
                alphabet,
                self._batch_config.system.distance_method,
                self._batch_config.system.enable_distance_cache,
                execution["algorithm_params"],
                self._combination_store.work_id,
                self._execution_controller.internal_jobs,
                execution["algorithm_name"],
                worker_config.get("cpu"),
            )
            results.append({execution["unit_id"]: result})
            
        return results
