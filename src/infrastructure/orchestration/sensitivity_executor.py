"""SensitivityExecutor: executes SensitivityTask with parameter sampling.

Reimplemented to follow ExperimentExecutor and OptimizationExecutor pattern:
- Each sample is stored as an execution in the executions table
- Implements resource control and status checking
- Sample parallelization
- Full integration with persistence system
- Sensitivity analysis based on result variance
"""

from __future__ import annotations

import random
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import variance
from typing import Any

from src.application.ports.repositories import AbstractExecutionEngine
from src.domain.config import AlgParams, CSPBenchConfig, SensitivityTask
from src.domain.dataset import Dataset
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
    WorkScopedPersistence,
)

from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.SensitivityExecutor")


def _sample_worker(
    sensitivity_unit_id: str,
    sample_number: int,
    sample_params: dict[str, Any],
    strings: list[str],
    alphabet: str,
    distance_method: str,
    use_cache: bool,
    internal_jobs: int,
    algorithm_name: str,
    cpu_config: dict[str, Any] | None = None,
    work_id: str | None = None,
) -> tuple[dict[str, Any], str]:
    """Function executed in subprocess for a sensitivity analysis sample (needs to be picklable).
    
    Args:
        sensitivity_unit_id: Unique identifier for this sensitivity unit
        sample_number: Sample sequence number
        sample_params: Parameters for this sample
        strings: Input strings for the algorithm
        alphabet: Alphabet used
        distance_method: Distance calculation method
        use_cache: Whether to use distance cache
        internal_jobs: Number of internal jobs for algorithms
        algorithm_name: Name of algorithm to execute
        cpu_config: CPU configuration for worker process
        work_id: Work identifier for status checking
        
    Returns:
        Tuple of (sample_result_dict, sensitivity_unit_id)
    """
    import psutil

    from src.domain.distance import create_distance_calculator
    from src.domain.status import BaseStatus as _BaseStatus
    from src.infrastructure.execution_control import ExecutionController
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

    logger = get_logger("CSPBench.SampleWorker")

    # Apply CPU configuration to this worker process
    if cpu_config:
        try:
            current_process = psutil.Process()

            # Apply CPU affinity if specified
            if cpu_config.get("affinity"):
                current_process.cpu_affinity(cpu_config["affinity"])
                logger.info(
                    f"[SAMPLE-WORKER] CPU affinity set to: {cpu_config['affinity']}"
                )

            # Apply CPU priority if max_workers is limited
            if cpu_config.get("exclusive_cores", False):
                # Lower priority for worker processes when using exclusive cores
                current_process.nice(5)  # Slightly lower priority than main process
                logger.info("[SAMPLE-WORKER] CPU priority adjusted for exclusive cores")

        except Exception as e:
            logger.warning(f"[SAMPLE-WORKER] Cannot apply CPU configuration: {e}")

    # Apply memory limits if specified in cpu_config
    if cpu_config and cpu_config.get("max_memory_gb"):
        try:
            import resource

            max_memory_bytes = int(cpu_config["max_memory_gb"] * 1024 * 1024 * 1024)
            resource.setrlimit(resource.RLIMIT_AS, (max_memory_bytes, max_memory_bytes))
            logger.info(
                f"[SAMPLE-WORKER] Memory limit set to {cpu_config['max_memory_gb']}GB"
            )
        except Exception as e:
            logger.warning(f"[SAMPLE-WORKER] Cannot apply memory limit: {e}")

    # Recreate store and wrappers using factory methods
    store = WorkPersistence()
    work_scoped = WorkScopedPersistence(work_id, store)
    execution_store = work_scoped.for_execution(sensitivity_unit_id)

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

    distance_calc = create_distance_calculator(
        distance_method=distance_method, strings=strings, use_cache=use_cache
    )

    try:
        execution_store.update_execution_status(_BaseStatus.RUNNING)
        monitor.on_progress(
            0.0, f"Starting sample {sample_number} of sensitivity analysis"
        )

        result = run_algorithm(
            algorithm_name=algorithm_name,
            strings=strings,
            alphabet=alphabet,
            distance_calculator=distance_calc,
            execution_controller=dummy_controller,
            monitor=monitor,
            seed=sample_params.get("seed"),
            params=sample_params,
        )

        monitor.on_progress(1.0, f"Sample {sample_number} finished")
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

        # Return sample result with objective for sensitivity analysis
        sample_result = dict(result)
        sample_result["sample_number"] = sample_number
        sample_result["sample_params"] = sample_params
        sample_result["objective"] = objective

        return sample_result, sensitivity_unit_id

    except Exception as exc:  # pragma: no cover - defensive
        execution_store.update_execution_status("failed", result={"error": str(exc)})
        return {
            "status": _BaseStatus.FAILED,
            "error": str(exc),
            "sample_number": sample_number,
            "sample_params": sample_params,
        }, sensitivity_unit_id


class SensitivityExecutor(AbstractExecutionEngine):
    """Executor for sensitivity analysis tasks."""

    def __init__(
        self,
        combination_store: CombinationScopedPersistence,
        execution_controller: ExecutionController,
        batch_config: CSPBenchConfig,
    ):
        """Initialize sensitivity executor.
        
        Args:
            combination_store: Scoped persistence for combination data
            execution_controller: Controller for execution limits and status
            batch_config: Batch configuration object
        """
        super().__init__(combination_store, execution_controller, batch_config)

    def run(
        self,
        task: SensitivityTask,
        dataset_obj: Dataset,
        alg: AlgParams,
    ) -> BaseStatus:
        """Execute sensitivity analysis by parameter sampling.
        
        Args:
            task: Sensitivity analysis task configuration
            dataset_obj: Dataset to analyze
            alg: Algorithm parameters configuration
            
        Returns:
            Final execution status
        """

        try:
            # Extract configuration
            config = task.config or {}
            samples = config.get("samples", 32)  # Default increased for better analysis
            method = task.method or "morris"  # Default sensitivity method

            logger.info("Starting sensitivity analysis for combination")
            logger.info(f"Config: samples={samples}, method={method}")

            strings = dataset_obj.sequences
            alphabet = dataset_obj.alphabet
            dataset_id = getattr(dataset_obj, "id", None)

            # Get base parameters from algorithm configuration
            base_params = dict(alg.params or {})

            # Get sensitivity parameters configuration
            # Check if parameters are nested by algorithm name
            sensitivity_params = task.parameters or {}

            # If parameters are nested by algorithm name, extract the correct subset
            if alg.name in sensitivity_params:
                sensitivity_params = sensitivity_params[alg.name]
                logger.debug(
                    f"Using algorithm-specific parameters for {alg.name}: {list(sensitivity_params.keys())}"
                )
            elif sensitivity_params and isinstance(
                next(iter(sensitivity_params.values())), dict
            ):
                # Check if any of the top-level keys match algorithm names
                for key, value in sensitivity_params.items():
                    if isinstance(value, dict) and "type" not in value:
                        # This looks like an algorithm name with nested parameters
                        if key == alg.name:
                            sensitivity_params = value
                            logger.debug(
                                f"Found nested parameters for algorithm {alg.name}"
                            )
                            break
                else:
                    # No algorithm-specific parameters found, use all parameters
                    logger.debug(
                        "No algorithm-specific parameters found, using all parameters"
                    )

            # Validate that we have parameters to analyze
            if not sensitivity_params:
                logger.warning(
                    f"No sensitivity parameters found for algorithm {alg.name}"
                )
                return BaseStatus.ERROR

            # Get global seed from system config
            global_seed = getattr(self._batch_config.system, "global_seed", None)
            if global_seed is None:
                global_seed = base_params.get("seed", None)

            rng = random.Random(global_seed)

            logger.info(
                f"Starting sensitivity analysis: {samples} samples, method: {method}"
            )

            try:
                max_workers = self._execution_controller.get_worker_config()["cpu"][
                    "max_workers"
                ]
            except Exception:
                max_workers = 1

            completed_samples = 0
            failed_samples = 0
            results: list[dict[str, Any]] = []

            logger.info(f"Configuration: max_workers={max_workers}")

            if max_workers > 1:
                # Parallel execution
                completed_samples, failed_samples, results = (
                    self._run_parallel_sensitivity(
                        task,
                        dataset_obj,
                        alg,
                        samples,
                        base_params,
                        sensitivity_params,
                        strings,
                        alphabet,
                        dataset_id,
                        max_workers,
                        rng,
                    )
                )
            else:
                # Sequential execution
                completed_samples, failed_samples, results = (
                    self._run_sequential_sensitivity(
                        task,
                        dataset_obj,
                        alg,
                        samples,
                        base_params,
                        sensitivity_params,
                        strings,
                        alphabet,
                        dataset_id,
                        rng,
                    )
                )

            # Check final status
            if self._execution_controller.check_status() != BaseStatus.RUNNING:
                return self._execution_controller.check_status()

            # Perform sensitivity analysis if we have enough results
            if (
                completed_samples >= 2
            ):  # Need at least 2 samples for variance calculation
                self._perform_sensitivity_analysis(results, sensitivity_params, method)
                logger.info(
                    f"Sensitivity analysis completed: {completed_samples} samples completed, {failed_samples} failed"
                )
                return BaseStatus.COMPLETED
            elif completed_samples == 0:
                logger.error("No sample was completed successfully")
                return BaseStatus.FAILED
            else:
                logger.warning(
                    f"Too few samples for valid analysis: {completed_samples} completed, {failed_samples} failed"
                )
                return BaseStatus.ERROR

        except Exception as e:
            logger.error(f"Error during sensitivity analysis: {e}")
            logger.error(f"Error type: {type(e)}")
            logger.error("Full stack trace:", exc_info=True)
            return BaseStatus.FAILED

    def _generate_sample_params(
        self,
        base_params: dict[str, Any],
        sensitivity_params: dict[str, Any],
        rng: random.Random,
    ) -> dict[str, Any]:
        """Generate sample parameters based on sensitivity configuration.
        
        Args:
            base_params: Base algorithm parameters
            sensitivity_params: Sensitivity parameter configuration
            rng: Random number generator
            
        Returns:
            Generated parameters for this sample
        """
        sample_params = dict(base_params)

        for param_name, param_config in sensitivity_params.items():
            param_type = param_config.get("type", "uniform")

            if param_type == "categorical":
                if "values" in param_config:
                    sample_params[param_name] = rng.choice(param_config["values"])
                elif "choices" in param_config:
                    sample_params[param_name] = rng.choice(param_config["choices"])
            elif param_type in ("int", "integer"):
                if "bounds" in param_config:
                    low, high = param_config["bounds"]
                    sample_params[param_name] = rng.randint(low, high)
                elif "low" in param_config and "high" in param_config:
                    sample_params[param_name] = rng.randint(
                        param_config["low"], param_config["high"]
                    )
            elif param_type == "float":
                if "bounds" in param_config:
                    low, high = param_config["bounds"]
                    sample_params[param_name] = rng.uniform(low, high)
                elif "low" in param_config and "high" in param_config:
                    sample_params[param_name] = rng.uniform(
                        param_config["low"], param_config["high"]
                    )
            elif param_type == "uniform":
                if "bounds" in param_config:
                    low, high = param_config["bounds"]
                    sample_params[param_name] = rng.uniform(low, high)
                elif "low" in param_config and "high" in param_config:
                    sample_params[param_name] = rng.uniform(
                        param_config["low"], param_config["high"]
                    )
            elif param_type == "loguniform":
                import math

                if "bounds" in param_config:
                    low, high = param_config["bounds"]
                    log_low = math.log(low)
                    log_high = math.log(high)
                    sample_params[param_name] = math.exp(rng.uniform(log_low, log_high))
                elif "low" in param_config and "high" in param_config:
                    log_low = math.log(param_config["low"])
                    log_high = math.log(param_config["high"])
                    sample_params[param_name] = math.exp(rng.uniform(log_low, log_high))
            else:
                logger.warning(f"Unknown parameter type: {param_type}")

        return sample_params

    def _run_parallel_sensitivity(
        self,
        task,
        dataset_obj,
        alg,
        samples,
        base_params,
        sensitivity_params,
        strings,
        alphabet,
        dataset_id,
        max_workers,
        rng,
    ) -> tuple[int, int, list[dict[str, Any]]]:
        """Execute sensitivity analysis in parallel.
        
        Args:
            task: Sensitivity task configuration
            dataset_obj: Dataset object
            alg: Algorithm parameters
            samples: Number of samples to run
            base_params: Base algorithm parameters
            sensitivity_params: Sensitivity parameter configuration
            strings: Input strings
            alphabet: Alphabet used
            dataset_id: Dataset identifier
            max_workers: Maximum number of worker processes
            rng: Random number generator
            
        Returns:
            Tuple of (completed_samples, failed_samples, results)
        """

        # Get CPU configuration for worker processes
        cpu_config = self._execution_controller.create_worker_config()

        # Add memory configuration to be passed to workers
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        completed_samples = 0
        failed_samples = 0
        results = []

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []

            for sample_num in range(samples):
                if self._execution_controller.check_status() != BaseStatus.RUNNING:
                    break

                # Generate sample parameters
                sample_params = self._generate_sample_params(
                    base_params, sensitivity_params, rng
                )

                # Build standardized unit_id: type:task:dataset:config:name:sample_X
                config_id = (
                    getattr(self._combination_store, "_preset_id", None) or "default"
                )
                config_id = str(config_id).replace(":", "_")
                sensitivity_unit_id = f"sensitivity:{task.id}:{dataset_id}:{config_id}:{alg.name}:{sample_num}"

                # Check if this sample already exists and is completed
                ex = self._combination_store.get_executions(unit_id=sensitivity_unit_id)
                if ex and ex[0]["status"] in (
                    BaseStatus.COMPLETED.value,
                    "completed",
                    BaseStatus.FAILED.value,
                    "failed",
                ):
                    continue

                self._combination_store.submit_execution(
                    unit_id=sensitivity_unit_id, sequencia=sample_num
                )

                future = executor.submit(
                    _sample_worker,
                    sensitivity_unit_id,
                    sample_num,
                    sample_params,
                    strings,
                    alphabet,
                    self._batch_config.system.distance_method,
                    self._batch_config.system.enable_distance_cache,
                    self._execution_controller.internal_jobs,
                    alg.name,
                    worker_config.get("cpu"),
                    self._combination_store.work_id,
                )
                futures.append(future)

            # Collect results
            for future in futures:
                try:
                    result, unit_id = future.result()

                    if "objective" in result and result["objective"] is not None:
                        results.append(result)
                        completed_samples += 1
                    else:
                        failed_samples += 1

                except Exception as e:
                    logger.error(f"Sample failed: {e}")
                    failed_samples += 1

        return completed_samples, failed_samples, results

    def _run_sequential_sensitivity(
        self,
        task,
        dataset_obj,
        alg,
        samples,
        base_params,
        sensitivity_params,
        strings,
        alphabet,
        dataset_id,
        rng,
    ) -> tuple[int, int, list[dict[str, Any]]]:
        """Execute sensitivity analysis sequentially.
        
        Args:
            task: Sensitivity task configuration
            dataset_obj: Dataset object
            alg: Algorithm parameters
            samples: Number of samples to run
            base_params: Base algorithm parameters
            sensitivity_params: Sensitivity parameter configuration
            strings: Input strings
            alphabet: Alphabet used
            dataset_id: Dataset identifier
            rng: Random number generator
            
        Returns:
            Tuple of (completed_samples, failed_samples, results)
        """

        # Get CPU configuration for sequential execution
        cpu_config = self._execution_controller.create_worker_config()

        # Add memory configuration to be passed to workers
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        completed_samples = 0
        failed_samples = 0
        results = []

        for sample_num in range(samples):
            status = self._execution_controller.check_status()
            if status != BaseStatus.RUNNING:
                break

            # Generate sample parameters
            sample_params = self._generate_sample_params(
                base_params, sensitivity_params, rng
            )

            # Build standardized unit_id: type:task:dataset:config:name:sample_X
            config_id = (
                getattr(self._combination_store, "_preset_id", None) or "default"
            )
            config_id = str(config_id).replace(":", "_")
            sensitivity_unit_id = f"sensitivity:{task.id}:{dataset_id}:{config_id}:{alg.name}:{sample_num}"

            # Check if this sample already exists and is completed
            ex = self._combination_store.get_executions(unit_id=sensitivity_unit_id)
            if ex and ex[0]["status"] in (
                BaseStatus.COMPLETED.value,
                "completed",
                BaseStatus.FAILED.value,
                "failed",
            ):
                continue

            self._combination_store.submit_execution(
                unit_id=sensitivity_unit_id, sequencia=sample_num
            )

            try:
                result, _ = _sample_worker(
                    sensitivity_unit_id,
                    sample_num,
                    sample_params,
                    strings,
                    alphabet,
                    self._batch_config.system.distance_method,
                    self._batch_config.system.enable_distance_cache,
                    self._execution_controller.internal_jobs,
                    alg.name,
                    worker_config.get("cpu"),
                    self._combination_store.work_id,
                )

                if "objective" in result and result["objective"] is not None:
                    results.append(result)
                    completed_samples += 1
                else:
                    failed_samples += 1

            except Exception as e:
                logger.error(f"Sample {sample_num} failed: {e}")
                failed_samples += 1

        return completed_samples, failed_samples, results

    def _perform_sensitivity_analysis(
        self,
        results: list[dict[str, Any]],
        sensitivity_params: dict[str, Any],
        method: str,
    ) -> dict[str, float]:
        """Perform sensitivity analysis based on sample results.
        
        Args:
            results: List of sample execution results
            sensitivity_params: Sensitivity parameter configuration
            method: Analysis method to use
            
        Returns:
            Dictionary of parameter sensitivity scores
        """

        if len(results) < 2:
            logger.warning("Too few results for sensitivity analysis")
            return {}

        # Extract objectives from results
        objectives = [
            r.get("objective", float("inf"))
            for r in results
            if r.get("objective") is not None
        ]

        if len(objectives) < 2:
            logger.warning("Too few valid objectives for sensitivity analysis")
            return {}

        # Simple variance-based sensitivity analysis
        param_scores: dict[str, float] = {}
        overall_variance = variance(objectives) if len(objectives) > 1 else 0.0

        logger.info(f"Overall objective variance: {overall_variance}")

        # For each parameter, calculate sensitivity based on parameter variation impact
        for param_name in sensitivity_params.keys():
            param_values = []
            param_objectives = []

            for result in results:
                sample_params = result.get("sample_params", {})
                if param_name in sample_params and result.get("objective") is not None:
                    param_value = sample_params[param_name]
                    # Only include numeric values for correlation analysis
                    if isinstance(param_value, (int, float)) and not isinstance(param_value, bool):
                        param_values.append(float(param_value))
                        param_objectives.append(float(result["objective"]))

            if len(param_values) >= 2:
                # Calculate correlation between parameter values and objectives
                try:
                    import numpy as np

                    # Ensure both arrays are numpy arrays with numeric dtype
                    param_values_array = np.array(param_values, dtype=np.float64)
                    param_objectives_array = np.array(param_objectives, dtype=np.float64)
                    
                    correlation = np.corrcoef(param_values_array, param_objectives_array)[0, 1]
                    param_scores[param_name] = (
                        abs(correlation) if not np.isnan(correlation) else 0.0
                    )
                except ImportError:
                    # Fallback to simple variance calculation if numpy not available
                    param_scores[param_name] = overall_variance
                except (ValueError, TypeError) as e:
                    logger.warning(f"Error calculating correlation for {param_name}: {e}")
                    param_scores[param_name] = overall_variance
            else:
                param_scores[param_name] = 0.0

        logger.info(f"Calculated sensitivity scores: {param_scores}")

        # Store sensitivity analysis results
        self._store_sensitivity_results(param_scores, method, len(results))

        return param_scores

    def _store_sensitivity_results(
        self, param_scores: dict[str, float], method: str, num_samples: int
    ) -> None:
        """Store sensitivity analysis results.
        
        Args:
            param_scores: Parameter sensitivity scores
            method: Analysis method used
            num_samples: Number of samples analyzed
        """

        sensitivity_results = {
            "method": method,
            "num_samples": num_samples,
            "param_scores": param_scores,
            "analysis_date": time.time(),
        }

        # Store in combination-level metadata or events
        try:
            # Use the combination store wrapper's generic_event method
            self._combination_store.generic_event(
                unit_id="sensitivity_analysis",  # Use a descriptive unit_id for this event type
                event_type="progress",
                message=f"Sensitivity analysis completed using {method} method",
                context=sensitivity_results,
            )
            logger.info("Sensitivity analysis results stored")
        except Exception as e:
            logger.warning(
                f"Error storing sensitivity analysis results: {e}"
            )
