"""OptimizationExecutor: executes OptimizationTask with Optuna.

Refactored to follow ExperimentExecutor pattern:
- Each trial is stored as an execution in the executions table
- Implements resource control and status checking
- Supports trial parallelization
- Full integration with persistence system
"""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Optional

try:
    import optuna
    from optuna.pruners import MedianPruner, SuccessiveHalvingPruner
    from optuna.samplers import CmaEsSampler, RandomSampler, TPESampler

    OptunaTrial = optuna.trial.Trial
    OptunaStudy = optuna.Study
except ImportError:
    optuna = None
    TPESampler = RandomSampler = CmaEsSampler = None
    MedianPruner = SuccessiveHalvingPruner = None
    OptunaTrial = Any
    OptunaStudy = Any

from src.application.ports.repositories import AbstractExecutionEngine
from src.domain.config import AlgParams, OptimizationTask
from src.domain.dataset import Dataset
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger

from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.OptimizationExecutor")


def _trial_worker(
    optimization_unit_id: str,
    trial_number: int,
    trial_params: dict[str, Any],
    strings: list[str],
    alphabet: str,
    distance_method: str,
    use_cache: bool,
    base_params: dict[str, Any],
    internal_jobs: int,
    algorithm_name: str,
    cpu_config: dict[str, Any] | None = None,
    work_id: str | None = None,
) -> tuple[dict[str, Any], str]:
    """Function executed in subprocess for an optimization trial (needs to be picklable).
    
    Args:
        optimization_unit_id: Unique identifier for this optimization unit
        trial_number: Trial sequence number
        trial_params: Parameters for this trial
        strings: Input strings for the algorithm
        alphabet: Alphabet used
        distance_method: Distance calculation method
        use_cache: Whether to use distance cache
        base_params: Base algorithm parameters
        internal_jobs: Number of internal jobs for algorithms
        algorithm_name: Name of algorithm to execute
        cpu_config: CPU configuration for worker process
        work_id: Work identifier for status checking
        
    Returns:
        Tuple of (trial_result_dict, optimization_unit_id)
    """
    import psutil

    from src.domain.distance import create_distance_calculator
    from src.domain.status import BaseStatus as _BaseStatus
    from src.infrastructure.logging_config import get_logger
    from src.infrastructure.monitoring.persistence_monitor import (
        PersistenceMonitor,
    )
    from src.infrastructure.persistence.work_state.core import (  # local import
        WorkPersistence,
    )
    from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
        ExecutionScopedPersistence,
    )
    from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
        WorkScopedPersistence,
    )

    logger = get_logger("CSPBench.TrialWorker")

    # Apply CPU configuration to this worker process
    if cpu_config:
        try:
            current_process = psutil.Process()

            # Apply CPU affinity if specified
            if cpu_config.get("affinity"):
                current_process.cpu_affinity(cpu_config["affinity"])
                logger.info(
                    f"[TRIAL-WORKER] CPU affinity set to: {cpu_config['affinity']}"
                )

            # Apply CPU priority if max_workers is limited
            if cpu_config.get("exclusive_cores", False):
                # Lower priority for worker processes when using exclusive cores
                current_process.nice(5)  # Slightly lower priority than main process
                logger.info("[TRIAL-WORKER] CPU priority adjusted for exclusive cores")

        except Exception as e:
            logger.warning(f"[TRIAL-WORKER] Cannot apply CPU configuration: {e}")

    # Apply memory limits if specified in cpu_config
    if cpu_config and cpu_config.get("max_memory_gb"):
        try:
            import resource

            max_memory_bytes = int(cpu_config["max_memory_gb"] * 1024 * 1024 * 1024)
            resource.setrlimit(resource.RLIMIT_AS, (max_memory_bytes, max_memory_bytes))
            logger.info(
                f"[TRIAL-WORKER] Memory limit set to {cpu_config['max_memory_gb']}GB"
            )
        except Exception as e:
            logger.warning(f"[TRIAL-WORKER] Cannot apply memory limit: {e}")

    # Recriar store e wrappers using factory methods
    store = WorkPersistence()
    work_scoped = WorkScopedPersistence(work_id, store)
    execution_store = work_scoped.for_execution(optimization_unit_id)

    # Controller with work_id for status checks
    dummy_controller = ExecutionController(work_id=work_id)
    dummy_controller._internal_jobs = internal_jobs  # type: ignore[attr-defined]

    # Apply CPU configuration manually if provided (since we can't pass ResourcesConfig to subprocess)
    if cpu_config:
        dummy_controller._exclusive_cores = cpu_config.get("exclusive_cores", False)
        if "max_workers" in cpu_config:
            dummy_controller._max_workers = cpu_config["max_workers"]
            dummy_controller._current_workers = cpu_config["max_workers"]

    # Create monitor with controller for cancellation checks
    monitor = PersistenceMonitor(execution_store, execution_controller=dummy_controller)

    # Combine base parameters with trial parameters
    final_params = {**base_params, **trial_params}

    distance_calc = create_distance_calculator(
        distance_method=distance_method, strings=strings, use_cache=use_cache
    )

    try:
        execution_store.update_execution_status(_BaseStatus.RUNNING)
        monitor.on_progress(0.0, f"Starting trial {trial_number} of algorithm")

        result = run_algorithm(
            algorithm_name=algorithm_name,
            strings=strings,
            alphabet=alphabet,
            distance_calculator=distance_calc,
            execution_controller=dummy_controller,
            monitor=monitor,
            seed=final_params.get("seed"),
            params=final_params,
        )

        monitor.on_progress(1.0, f"Trial {trial_number} finished")
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

        # Return trial result with objective for Optuna
        trial_result = dict(result)
        trial_result["trial_number"] = trial_number
        trial_result["trial_params"] = trial_params
        trial_result["objective"] = objective

        return trial_result, optimization_unit_id

    except Exception as exc:  # pragma: no cover - defensive
        execution_store.update_execution_status("failed", result={"error": str(exc)})
        return {
            "status": _BaseStatus.FAILED,
            "error": str(exc),
            "trial_number": trial_number,
            "trial_params": trial_params,
        }, optimization_unit_id


class OptimizationExecutor(AbstractExecutionEngine):
    """Executor for optimization tasks using Optuna."""

    def run(
        self,
        task: OptimizationTask,
        dataset_obj: Dataset,
        alg: AlgParams,
    ) -> BaseStatus:
        """Execute optimization using Optuna.
        
        Args:
            task: Optimization task configuration
            dataset_obj: Dataset to optimize on
            alg: Algorithm parameters configuration
            
        Returns:
            Final execution status
        """

        if optuna is None:
            logger.error("Optuna is not installed. Run: pip install optuna")
            return BaseStatus.FAILED

        try:
            # Extract configuration
            config = task.config or {}
            trials = config.get("trials", 50)
            direction = config.get("direction", "minimize")
            timeout_per_trial = config.get("timeout_per_trial", 300)

            # Create unique study name based on current combination
            # Use combination information: task_id, dataset_id, preset_id, algorithm_id
            combination_info = self._get_combination_info()
            study_name = self._create_study_name(task, combination_info)

            logger.info(f"Starting optimization for combination: {study_name}")
            logger.info(f"Config: trials={trials}, direction={direction}")

            strings = dataset_obj.sequences
            alphabet = dataset_obj.alphabet
            dataset_id = getattr(dataset_obj, "id", None)

            # Get base parameters from algorithm configuration
            base_params = dict(alg.params or {})

            # Get optimization parameters configuration
            optimization_params = task.parameters or {}

            # Setup Optuna study with unique name per combination
            study = self._create_optuna_study(study_name, direction, config)

            logger.info(f"Optuna study created: {study_name}")

            logger.info(f"Starting optimization: {trials} trials, direction: {direction}")

            try:
                max_workers = self._execution_controller.get_worker_config()["cpu"][
                    "max_workers"
                ]
            except Exception:
                max_workers = 1

            completed_trials = 0
            failed_trials = 0

            logger.info(f"Configuration: max_workers={max_workers}")

            if max_workers > 1:
                # Parallel execution
                completed_trials, failed_trials = self._run_parallel_optimization(
                    study,
                    task,
                    dataset_obj,
                    alg,
                    trials,
                    timeout_per_trial,
                    strings,
                    alphabet,
                    dataset_id,
                    base_params,
                    optimization_params,
                    max_workers,
                )
            else:
                # Sequential execution
                completed_trials, failed_trials = self._run_sequential_optimization(
                    study,
                    task,
                    dataset_obj,
                    alg,
                    trials,
                    timeout_per_trial,
                    strings,
                    alphabet,
                    dataset_id,
                    base_params,
                    optimization_params,
                )

            # Check final status
            if self._execution_controller.check_status() != BaseStatus.RUNNING:
                return self._execution_controller.check_status()

            # Determine final result
            if completed_trials == 0:
                logger.error("No trial was completed successfully")
                return BaseStatus.FAILED
            elif failed_trials > completed_trials:
                logger.warning(
                    f"More trials failed ({failed_trials}) than completed ({completed_trials})"
                )
                return BaseStatus.ERROR
            else:
                logger.info(
                    f"Optimization completed: {completed_trials} trials completed, {failed_trials} failed"
                )
                return BaseStatus.COMPLETED

        except Exception as e:
            logger.error(f"Error during optimization: {e}")
            logger.error(f"Error type: {type(e)}")
            logger.error("Full stack trace:", exc_info=True)
            return BaseStatus.FAILED

    def _create_optuna_study(
        self, study_name: str, direction: str, config: dict
    ) -> Any:
        """Create an Optuna study with specified configuration.
        
        Args:
            study_name: Name for the study
            direction: Optimization direction ('minimize' or 'maximize')
            config: Configuration dictionary with sampler/pruner settings
            
        Returns:
            Configured Optuna study object
        """

        # Configure sampler
        sampler_name = config.get("sampler", "TPESampler")
        if sampler_name == "TPESampler":
            sampler = TPESampler(
                n_startup_trials=config.get("n_startup_trials", 10),
                n_ei_candidates=config.get("n_ei_candidates", 24),
                multivariate=config.get("multivariate", False),
            )
        elif sampler_name == "RandomSampler":
            sampler = RandomSampler(seed=config.get("seed"))
        elif sampler_name == "CmaEsSampler":
            sampler = CmaEsSampler(
                n_startup_trials=config.get("n_startup_trials", 1),
                restart_strategy=config.get("restart_strategy", "ipop"),
            )
        else:
            sampler = TPESampler()

        # Configure pruner
        pruner_name = config.get("pruner", "MedianPruner")
        if pruner_name == "MedianPruner":
            pruner = MedianPruner(
                n_startup_trials=config.get("n_startup_trials", 5),
                n_warmup_steps=config.get("n_warmup_steps", 10),
                interval_steps=config.get("interval_steps", 5),
            )
        elif pruner_name == "SuccessiveHalvingPruner":
            pruner = SuccessiveHalvingPruner(
                min_resource=config.get("min_resource", 1),
                reduction_factor=config.get("reduction_factor", 4),
                min_early_stopping_rate=config.get("min_early_stopping_rate", 0),
            )
        else:
            pruner = MedianPruner()

        # Create study
        # Configure storage baseado na configuração
        storage_config = config.get("storage")
        storage = self._configure_storage(storage_config, study_name)

        # Suppress verbose Optuna logs
        optuna_logger = logging.getLogger("optuna")
        original_level = optuna_logger.level
        optuna_logger.setLevel(logging.WARNING)

        try:
            study = optuna.create_study(
                study_name=study_name,
                direction=direction,
                sampler=sampler,
                pruner=pruner,
                storage=storage,
                load_if_exists=True,
            )
        finally:
            # Restore original logger level
            optuna_logger.setLevel(original_level)

        logger.info(f"Optuna study created: {study_name}")
        return study

    def _generate_trial_params(
        self, trial: Any, optimization_params: dict
    ) -> dict[str, Any]:
        """Generate trial parameters based on optimization configuration.
        
        Args:
            trial: Optuna trial object
            optimization_params: Parameter space configuration
            
        Returns:
            Generated parameters for this trial
        """
        params = {}

        # Handle case where optimization_params might contain algorithm name as key
        # Check if first level contains algorithm names (like 'BLF-GA', 'CSC')
        first_key = (
            next(iter(optimization_params.keys())) if optimization_params else None
        )
        if first_key and isinstance(optimization_params[first_key], dict):
            # Check if this looks like algorithm-specific config
            first_value = optimization_params[first_key]
            if isinstance(first_value, dict) and "type" not in first_value:
                # This is likely an algorithm name, extract the actual parameters
                # For now, take the first algorithm's parameters
                optimization_params = first_value

        for param_name, param_config in optimization_params.items():
            param_type = param_config.get("type", "uniform")

            if param_type == "categorical":
                params[param_name] = trial.suggest_categorical(
                    param_name, param_config["choices"]
                )
            elif param_type == "int":
                params[param_name] = trial.suggest_int(
                    param_name,
                    param_config["low"],
                    param_config["high"],
                    step=param_config.get("step", 1),
                )
            elif param_type == "float":
                params[param_name] = trial.suggest_float(
                    param_name,
                    param_config["low"],
                    param_config["high"],
                    step=param_config.get("step", None),
                )
            elif param_type == "uniform":
                params[param_name] = trial.suggest_float(
                    param_name, param_config["low"], param_config["high"]
                )
            elif param_type == "loguniform":
                params[param_name] = trial.suggest_loguniform(
                    param_name, param_config["low"], param_config["high"]
                )
            else:
                logger.warning(f"Unknown parameter type: {param_type}")

        return params

    def _run_parallel_optimization(
        self,
        study,
        task,
        dataset_obj,
        alg,
        trials,
        timeout_per_trial,
        strings,
        alphabet,
        dataset_id,
        base_params,
        optimization_params,
        max_workers,
    ) -> tuple[int, int]:
        """Execute optimization in parallel.
        
        Args:
            study: Optuna study object
            task: Optimization task configuration
            dataset_obj: Dataset object
            alg: Algorithm parameters
            trials: Number of trials to run
            timeout_per_trial: Timeout per trial in seconds
            strings: Input strings
            alphabet: Alphabet used
            dataset_id: Dataset identifier
            base_params: Base algorithm parameters
            optimization_params: Optimization parameter configuration
            max_workers: Maximum number of worker processes
            
        Returns:
            Tuple of (completed_trials, failed_trials)
        """

        # Get CPU configuration for worker processes
        cpu_config = self._execution_controller.create_worker_config()

        # Add memory configuration to be passed to workers
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        completed_trials = 0
        failed_trials = 0

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []

            for trial_num in range(trials):
                if self._execution_controller.check_status() != BaseStatus.RUNNING:
                    break

                # Create a trial to get parameters
                trial = study.ask()
                trial_params = self._generate_trial_params(trial, optimization_params)

                # Build standardized unit_id: type:task:dataset:config:name:trial_X
                config_id = (
                    getattr(self._combination_store, "_preset_id", None) or "default"
                )
                config_id = str(config_id).replace(":", "_")
                optimization_unit_id = f"optimization:{task.id}:{dataset_id}:{config_id}:{alg.name}:{trial.number}"

                # Check if this trial already exists and is completed
                ex = self._combination_store.get_executions(
                    unit_id=optimization_unit_id
                )
                if ex and ex[0]["status"] in (
                    BaseStatus.COMPLETED.value,
                    "completed",
                    BaseStatus.FAILED.value,
                    "failed",
                ):
                    continue

                self._combination_store.submit_execution(
                    unit_id=optimization_unit_id, sequencia=trial.number
                )

                future = executor.submit(
                    _trial_worker,
                    optimization_unit_id,
                    trial.number,
                    trial_params,
                    strings,
                    alphabet,
                    self._batch_config.system.distance_method,
                    self._batch_config.system.enable_distance_cache,
                    base_params,
                    self._execution_controller.internal_jobs,
                    alg.name,
                    worker_config.get("cpu"),
                    self._combination_store.work_id,
                )
                futures.append((trial, future))

            # Collect results
            for trial, future in futures:
                try:
                    result, unit_id = future.result(timeout=timeout_per_trial)

                    # Report result back to Optuna
                    if "objective" in result and result["objective"] is not None:
                        study.tell(trial, result["objective"])
                        completed_trials += 1
                    else:
                        study.tell(trial, float("inf"))  # Failed trial
                        failed_trials += 1

                except Exception as e:
                    logger.error(f"Trial {trial.number} failed: {e}")
                    study.tell(trial, float("inf"))  # Failed trial
                    failed_trials += 1

        return completed_trials, failed_trials

    def _run_sequential_optimization(
        self,
        study,
        task,
        dataset_obj,
        alg,
        trials,
        timeout_per_trial,
        strings,
        alphabet,
        dataset_id,
        base_params,
        optimization_params,
    ) -> tuple[int, int]:
        """Execute optimization sequentially.
        
        Args:
            study: Optuna study object
            task: Optimization task configuration
            dataset_obj: Dataset object
            alg: Algorithm parameters
            trials: Number of trials to run
            timeout_per_trial: Timeout per trial in seconds
            strings: Input strings
            alphabet: Alphabet used
            dataset_id: Dataset identifier
            base_params: Base algorithm parameters
            optimization_params: Optimization parameter configuration
            
        Returns:
            Tuple of (completed_trials, failed_trials)
        """

        # Get CPU configuration for sequential execution
        cpu_config = self._execution_controller.create_worker_config()

        # Add memory configuration to be passed to workers
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        completed_trials = 0
        failed_trials = 0

        for trial_num in range(trials):
            status = self._execution_controller.check_status()
            if status != BaseStatus.RUNNING:
                break

            # Create a trial to get parameters
            trial = study.ask()
            trial_params = self._generate_trial_params(trial, optimization_params)

            # Build standardized unit_id: type:task:dataset:config:name:trial_X
            config_id = (
                getattr(self._combination_store, "_preset_id", None) or "default"
            )
            config_id = str(config_id).replace(":", "_")
            optimization_unit_id = f"optimization:{task.id}:{dataset_id}:{config_id}:{alg.name}:{trial.number}"

            # Check if this trial already exists and is completed
            ex = self._combination_store.get_executions(unit_id=optimization_unit_id)
            if ex and ex[0]["status"] in (
                BaseStatus.COMPLETED.value,
                "completed",
                BaseStatus.FAILED.value,
                "failed",
            ):
                continue

            self._combination_store.submit_execution(
                unit_id=optimization_unit_id, sequencia=trial.number
            )

            try:
                result, _ = _trial_worker(
                    optimization_unit_id,
                    trial.number,
                    trial_params,
                    strings,
                    alphabet,
                    self._batch_config.system.distance_method,
                    self._batch_config.system.enable_distance_cache,
                    base_params,
                    self._execution_controller.internal_jobs,
                    alg.name,
                    worker_config.get("cpu"),
                    self._combination_store.work_id,
                )

                # Report result back to Optuna
                if "objective" in result and result["objective"] is not None:
                    study.tell(trial, result["objective"])
                    completed_trials += 1
                else:
                    study.tell(trial, float("inf"))  # Failed trial
                    failed_trials += 1

            except Exception as e:
                logger.error(f"Trial {trial.number} failed: {e}")
                study.tell(trial, float("inf"))  # Failed trial
                failed_trials += 1

        return completed_trials, failed_trials

    def _configure_storage(self, storage_config: Any, study_name: str) -> Optional[str]:
        """Configure Optuna storage based on provided configuration.

        Args:
            storage_config: Storage configuration (string, boolean or None)
            study_name: Study name for automatic filename generation

        Returns:
            SQLite connection string or None for in-memory storage
        """
        if not storage_config:
            return None  # In-memory storage

        # For any storage value (True or string), save in work directory
        from pathlib import Path
        from src.infrastructure.utils.path_utils import get_output_base_directory

        # Use current work's output directory (where results are stored)
        work_output_dir = get_output_base_directory() / self._combination_store.work_id
        optuna_dir = work_output_dir / "optuna"
        optuna_dir.mkdir(parents=True, exist_ok=True)

        # If specific string, use as filename
        if isinstance(storage_config, str) and not storage_config.startswith("sqlite:"):
            db_filename = storage_config
        else:
            # To keep all studies of a task in the same file,
            # use task_id instead of individual study_name
            task_id = self._combination_store._task_id or "unknown_task"
            db_filename = f"{task_id}_optuna_studies.db"

        storage_path = optuna_dir / db_filename
        return f"sqlite:///{storage_path}"

    def _get_combination_info(self) -> dict[str, str]:
        """Get current combination information.
        
        Returns:
            Dictionary with combination identifiers
        """
        return {
            "task_id": self._combination_store._task_id or "unknown",
            "dataset_id": self._combination_store._dataset_id or "unknown",
            "preset_id": self._combination_store._preset_id or "unknown",
            "algorithm_id": self._combination_store._algorithm_id or "unknown",
        }

    def _create_study_name(self, task: Any, combination_info: dict[str, str]) -> str:
        """Create unique study name based on combination.
        
        Args:
            task: Task configuration object
            combination_info: Dictionary with combination identifiers
            
        Returns:
            Unique study name string
        """
        task_id = combination_info["task_id"]
        dataset_id = combination_info["dataset_id"]
        preset_id = combination_info["preset_id"]
        algorithm_id = combination_info["algorithm_id"]

        # Format: task_id|dataset_id|preset_id|algorithm_id
        study_name = f"{task_id}|{dataset_id}|{preset_id}|{algorithm_id}"
        return study_name
