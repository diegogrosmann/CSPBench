"""PipelineRunner with control support (pause/cancel/resume) via WorkManager."""

from __future__ import annotations

import os
import time
from dataclasses import asdict
from typing import Any

from src.application.ports.repositories import ExecutionEngine
from src.application.services.dataset_service import load_dataset
from src.domain.config import (
    CSPBenchConfig,
    ExperimentTask,
    OptimizationTask,
    SensitivityTask,
)
from src.domain.dataset import Dataset
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import (
    ExecutionController,
    ExecutionLimitError,
)
from src.infrastructure.logging_config import get_logger
from src.infrastructure.persistence.work_state import (
    WorkScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)

from .experiment_executor import ExperimentExecutor
from .optimization_executor import OptimizationExecutor
from .sensitivity_executor import SensitivityExecutor

# Create module logger
logger = get_logger("CSPBench.PipelineRunner")


class PipelineRunner:
    """Main pipeline runner with support for pause/cancel/resume operations."""

    def __init__(self, work_store: WorkScopedPersistence):
        """Initialize pipeline runner.

        Args:
            work_store: Work-scoped persistence wrapper
        """
        self.work_store: WorkScopedPersistence = work_store
        self.work_id = work_store.work_id
        # ExecutionController will be initialized in run() with proper resources config
        self.execution_controller: ExecutionController = None
        self.datasets_cache: dict[str, Dataset] = {}

        logger.info(f"PipelineRunner initialized for work_id: {self.work_id}")

    def run(self, config: CSPBenchConfig) -> None:
        """Execute the pipeline following the new reorganized flow.

        Args:
            config: Complete CSPBench configuration
        """
        logger.info("Starting reorganized pipeline execution")
        logger.info(
            f"Configuration: {len(config.tasks.items)} tasks of type {config.tasks.type}"
        )

        # Initialize ExecutionController with resources configuration
        self.execution_controller = ExecutionController(
            work_id=self.work_id, resources=config.resources
        )

        # Apply resource limits to main process
        if config.resources:
            if config.resources.memory:
                self.execution_controller.apply_memory_limits()
            if config.resources.cpu:
                self.execution_controller.apply_cpu_limits()

        try:
            self.work_store.update_work_status(BaseStatus.RUNNING)

            # Phase 1: Generate all datasets first
            logger.info("=== PHASE 1: Generating datasets ===")
            self.execution_controller.check_batch_timeout()  # Check timeout before each phase
            self._generate_all_datasets(config)

            # Phase 2: Generate all combinations
            logger.info("=== PHASE 2: Generating combinations ===")
            self.execution_controller.check_batch_timeout()  # Check timeout before each phase
            self._generate_pipeline_combinations(config)

            # Phase 3: Execute combinations
            logger.info("=== PHASE 3: Executing combinations ===")
            self.execution_controller.check_batch_timeout()  # Check timeout before each phase
            status = self._execute_combinations(config)

            # Wait for all executions to finish before updating work status
            logger.info("Waiting for all executions in progress to complete...")
            self._wait_for_all_executions_to_complete()

            # Check if work was paused or canceled during execution
            current_work_status = self.work_store.get_work_status()

            if current_work_status in [BaseStatus.PAUSED, BaseStatus.CANCELED]:
                logger.info(
                    f"Pipeline was {current_work_status.lower()} during execution"
                )
                # Don't update status - keep as is (paused/canceled)
                self.final_status_to_set = None
            elif status == BaseStatus.RUNNING:
                logger.error("Pipeline still running after completion")
                self.final_status_to_set = BaseStatus.FAILED
            else:
                logger.info("Pipeline executed successfully")
                self.final_status_to_set = status

        except ExecutionLimitError as timeout_exc:
            logger.error(f"Pipeline interrupted by timeout: {timeout_exc}")
            if self.work_store:
                self.work_store.work_error(timeout_exc)
                self.work_store.update_work_status(
                    BaseStatus.CANCELED, error=str(timeout_exc)
                )
        except Exception as e:
            logger.error(f"Error during pipeline execution: {e}", exc_info=True)
            if self.work_store:
                self.work_store.work_error(e)
                self.work_store.update_work_status(BaseStatus.FAILED, error=str(e))
            raise e
        finally:
            # Cleanup execution controller
            if self.execution_controller:
                self.execution_controller.cleanup()

            # Check if final export should be done (before updating final status)
            final_work_status = self.work_store.get_work_status()
            should_export = final_work_status in [
                BaseStatus.COMPLETED,
                BaseStatus.ERROR,
                BaseStatus.FAILED,
            ]

            # If we have a pending status to set, also check if it requires export
            if hasattr(self, "final_status_to_set") and self.final_status_to_set:
                should_export = should_export or self.final_status_to_set in [
                    BaseStatus.COMPLETED,
                    BaseStatus.ERROR,
                    BaseStatus.FAILED,
                ]

            if should_export:
                logger.info(f"Starting final export (status: {final_work_status})")
                # Finalization export (raw db + manifest + full_results + summary)
                try:
                    from src.infrastructure.export.finalization_service import (
                        FinalizationConfig,
                        FinalizationService,
                    )

                    # Get output directory from work row (output_path field) if present
                    output_dir = None
                    try:
                        output_dir = self.work_store.get_work_output_path(
                            self.work_id
                        )  # may return Path or None
                    except Exception:
                        output_dir = None
                    if output_dir is None:
                        # fallback: derive from environment OUTPUT_BASE_DIRECTORY
                        from pathlib import Path

                        base = os.environ.get("OUTPUT_BASE_DIRECTORY", "./data/outputs")
                        output_dir = Path(base) / self.work_id
                    output_dir.mkdir(parents=True, exist_ok=True)

                    # Use WorkScopedPersistence instead of accessing db_path
                    cfg = FinalizationConfig(
                        work_id=self.work_id,
                        output_dir=output_dir,
                        tool_version=FinalizationService.detect_tool_version(),
                    )
                    FinalizationService(cfg, work_store=self.work_store).run()
                    logger.info("Final export completed successfully")
                except Exception as fe:  # pragma: no cover - defensive
                    logger.error(f"Error during export finalization: {fe}")
            else:
                logger.info(f"Final export skipped due to status: {final_work_status}")

            # Now update the final status after export is complete
            if hasattr(self, "final_status_to_set") and self.final_status_to_set:
                logger.info(f"Setting final work status to: {self.final_status_to_set}")
                self.work_store.update_work_status(self.final_status_to_set)

    def _generate_all_datasets(self, config: CSPBenchConfig) -> None:
        """Phase 1: Generate all necessary datasets.

        Args:
            config: CSPBench configuration containing dataset definitions
        """
        # Collect all unique dataset IDs used in tasks
        used_dataset_ids = set()

        for task in config.tasks.items:
            # task.datasets is now List[str] containing dataset IDs
            used_dataset_ids.update(task.datasets)

        logger.info(
            f"Found {len(used_dataset_ids)} unique dataset IDs used: {sorted(used_dataset_ids)}"
        )

        # Check if all referenced IDs exist in configuration
        missing_datasets = used_dataset_ids - set(config.datasets.keys())
        if missing_datasets:
            raise ValueError(
                f"Dataset IDs referenced but not defined: {sorted(missing_datasets)}"
            )

        # Process each used dataset
        for dataset_id in sorted(used_dataset_ids):
            logger.info(f"Processing dataset: {dataset_id}")

            try:
                # Get dataset configuration from config.datasets dictionary
                dataset_config = config.datasets[dataset_id]
                self._process_dataset(dataset_config)

                logger.info(f"Dataset {dataset_id} processed.")

            except Exception as e:
                logger.error(f"Error processing dataset {dataset_id}: {e}")
                raise

    def _generate_pipeline_combinations(self, config: CSPBenchConfig) -> None:
        """Phase 2: Generate all possible pipeline combinations.

        Args:
            config: CSPBench configuration containing tasks and algorithms
        """
        try:
            requeued = self.work_store.init_combination()
            if requeued:
                logger.info(
                    "Running/paused/canceled combinations restarted to 'queued'. Checking for new combinations..."
                )
                return  # Don't generate new combinations if there are ongoing ones

            # Always generate combinations: both for new cases and to add new ones after reset
            combinations = []

            for task in config.tasks.items:
                total_sequences = None

                if isinstance(task, ExperimentTask):
                    total_sequences = task.repetitions
                elif isinstance(task, OptimizationTask):
                    total_sequences = (
                        task.config.get("trials", 50) if task.config else 50
                    )
                elif isinstance(task, SensitivityTask):
                    total_sequences = (
                        task.config.get("samples", 100) if task.config else 100
                    )

                # task.datasets is now List[str] containing dataset IDs
                for dataset_id in task.datasets:
                    # task.algorithms is now List[str] containing algorithm preset IDs
                    for preset_id in task.algorithms:
                        # Get preset from config.algorithms dictionary
                        preset = config.algorithms.get(preset_id)
                        if not preset:
                            self.work_store.preset_error(
                                preset_id, "Algorithm not found in configuration"
                            )
                            logger.warning(
                                f"Algorithm preset '{preset_id}' not found in configuration"
                            )
                            continue

                        for alg in preset.items:
                            combination = {
                                "task_id": task.id,
                                "dataset_id": dataset_id,
                                "preset_id": preset_id,
                                "algorithm_id": alg.name,
                                "mode": task.type,
                                "total_sequences": total_sequences,
                            }
                            combinations.append(combination)

            logger.info(f"Generated {len(combinations)} combinations for processing")

            # Submit all combinations (INSERT OR IGNORE ensures duplicates are ignored)
            if combinations:
                inserted_count = self.work_store.submit_combinations(combinations)
                logger.info(f"{inserted_count} new combinations inserted into database")
            else:
                logger.warning("No valid combinations were generated")
        except Exception as e:
            self.work_store.combination_error("N/A", e)
            raise e

    def _execute_combinations(self, config: CSPBenchConfig) -> BaseStatus:
        """Phase 3: Execute all pending combinations.

        Args:
            config: CSPBench configuration

        Returns:
            Final execution status
        """
        if not self.work_store:
            logger.error("Work store not available for execution")
            return BaseStatus.FAILED

        execution_statuses = []  # List to store execution statuses

        while True:
            # If control indicates no longer RUNNING, stop
            if self.execution_controller.check_status() != BaseStatus.RUNNING:
                return self.execution_controller.check_status()

            # Get next pending combination
            combination = self.work_store.get_next_pending_combination()
            logger.debug(
                "[PipelineRunner] Next combination returned=%s",
                combination,
            )
            if not combination:
                logger.info("All combinations have been processed")
                break

            work_combination = self.work_store.for_combination(combination["id"])

            # Execute combination and store status
            status = self._execute_single_combination(
                combination, config, work_combination
            )
            execution_statuses.append(status)

        # Check status list and determine final result
        if not execution_statuses:
            # No combination was executed
            return BaseStatus.COMPLETED

        # If failure exists, return failure
        if BaseStatus.FAILED in execution_statuses:
            return BaseStatus.FAILED

        # If error exists, return error
        if BaseStatus.ERROR in execution_statuses:
            return BaseStatus.ERROR

        # If all completed, return completed
        if all(status == BaseStatus.COMPLETED for status in execution_statuses):
            return BaseStatus.COMPLETED

        # If something else, return error
        return BaseStatus.ERROR

    def _execute_single_combination(
        self,
        combination: dict[str, Any],
        config: CSPBenchConfig,
        work_combination: CombinationScopedPersistence,
    ) -> BaseStatus:
        """Execute a single combination.

        Args:
            combination: Combination data dictionary
            config: CSPBench configuration
            work_combination: Combination-scoped persistence wrapper

        Returns:
            Execution status for this combination
        """
        task_id = combination["task_id"]
        dataset_id = combination["dataset_id"]
        preset_id = combination["preset_id"]
        algorithm_id = combination["algorithm_id"]

        logger.info(
            f"Executing combination: {task_id}/{dataset_id}/{preset_id}/{algorithm_id}"
        )

        try:
            # Mark as running
            work_combination.update_combination_status(BaseStatus.RUNNING)

            # Find configuration objects
            task = self._find_task(config, task_id)
            if not task:
                raise ValueError(f"Task not found: {task_id}")

            alg = self._find_algorithm(config, preset_id, algorithm_id)
            if not alg:
                raise ValueError(f"Algorithm not found: {algorithm_id}")

            dataset_obj = self.datasets_cache.get(dataset_id)
            if not dataset_obj:
                raise ValueError(f"Dataset not found in cache: {dataset_id}")

            # Execute algorithm
            status = self._execute_algorithm(
                work_combination, task, alg, dataset_obj, config
            )

            work_combination.update_combination_status(status)

            logger.info(
                f"Combination {task_id}/{dataset_id}/{preset_id}/{algorithm_id} completed: {status}"
            )

            return status

        except Exception as e:
            logger.error(f"Error executing combination: {e}")
            work_combination.update_combination_status(BaseStatus.FAILED)
            work_combination.record_error(e)
            if self.work_store:
                self.work_store.algorithm_error(algorithm_id, e)
            return BaseStatus.FAILED

    def _execute_algorithm(
        self, work_combination, task, alg, dataset_obj, config
    ) -> BaseStatus:
        """Execute a single algorithm with proper result capture.

        Args:
            work_combination: Combination-scoped persistence wrapper
            task: Task configuration object
            alg: Algorithm parameters object
            dataset_obj: Dataset object to process
            config: Complete CSPBench configuration

        Returns:
            Algorithm execution status
        """
        try:
            engine: ExecutionEngine = self._get_executor(work_combination, config, task)

            result_status = engine.run(
                task=task,
                dataset_obj=dataset_obj,
                alg=alg,
            )

            return result_status

        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise e

    def _process_dataset(self, dataset_config) -> None:
        """Process a dataset (cache or generate new).

        Args:
            dataset_config: Dataset configuration object
        """
        dataset_id = getattr(dataset_config, "id", None)

        # Check cache first
        if self.work_store and dataset_id:
            try:
                dataset_obj = self.work_store.get_dataset(dataset_id)
                if dataset_obj is None or not dataset_obj.sequences:
                    logger.warning(f"Dataset {dataset_id} found in cache but is empty")
                else:
                    self.datasets_cache[dataset_id] = dataset_obj
                    return
            except Exception as e:
                logger.warning(f"Error loading dataset from cache {dataset_id}: {e}")

        # No cache or empty dataset, generate new
        try:
            # Generate/load dataset normally
            logger.info(f"Generating dataset: {dataset_id}")
            resolver, params = load_dataset(dataset_config)

            meta = {
                "batch_params": asdict(dataset_config) if dataset_config else {},
                "dataset_params": params if params else {},
                "dataset_statistics": resolver.get_statistics() if resolver else {},
            }

            self.work_store.submit_dataset(
                id=dataset_id, dataset_obj=resolver, meta=meta
            )
            self.datasets_cache[dataset_id] = resolver

        except Exception as e:
            logger.error(f"Error generating/saving dataset {dataset_id}: {e}")
            if self.work_store:
                self.work_store.dataset_error(dataset_id, e)
            raise ValueError(f"Error generating/saving dataset {dataset_id}: {e}")

    def _find_task(self, config: CSPBenchConfig, task_id: str):
        """Find a task by ID.

        Args:
            config: CSPBench configuration
            task_id: Task identifier to find

        Returns:
            Task object or None if not found
        """
        for task in config.tasks.items:
            if task.id == task_id:
                return task
        return None

    def _find_algorithm(
        self, config: CSPBenchConfig, preset_id: str, algorithm_id: str
    ):
        """Find an algorithm within a preset using the new structure.

        Args:
            config: CSPBench configuration
            preset_id: Algorithm preset identifier
            algorithm_id: Algorithm identifier within the preset

        Returns:
            Algorithm object or None if not found
        """
        preset = config.algorithms.get(preset_id)
        if not preset:
            return None

        for alg in preset.items:
            if alg.name == algorithm_id:
                return alg
        return None

    def _get_executor(self, work_combination, batch_config, task) -> ExecutionEngine:
        """Factory method to get appropriate execution engine based on task type.

        Args:
            work_combination: Combination-scoped persistence wrapper
            batch_config: Batch configuration object
            task: Task configuration object

        Returns:
            Appropriate execution engine instance
        """
        if isinstance(task, ExperimentTask):
            return ExperimentExecutor(
                combination_store=work_combination,
                execution_controller=self.execution_controller,
                batch_config=batch_config,
            )
        elif isinstance(task, OptimizationTask):
            return OptimizationExecutor(
                combination_store=work_combination,
                execution_controller=self.execution_controller,
                batch_config=batch_config,
            )
        elif isinstance(task, SensitivityTask):
            return SensitivityExecutor(
                combination_store=work_combination,
                execution_controller=self.execution_controller,
                batch_config=batch_config,
            )
        else:
            raise ValueError(f"Unsupported task type: {type(task)}")

    def _wait_for_all_executions_to_complete(self, max_wait_time: int = 60) -> None:
        """Wait for all executions in progress to finish before finalizing pipeline.

        Args:
            max_wait_time: Maximum wait time in seconds
        """
        logger.info("Checking executions in progress...")

        start_time = time.time()
        check_interval = 2.0  # Check every 2 seconds

        while time.time() - start_time < max_wait_time:
            # Search for executions still in progress
            running_executions = self.work_store.get_running_executions()

            if not running_executions:
                logger.info("All executions have been completed")
                return

            logger.info(
                f"Waiting for {len(running_executions)} executions in progress to complete..."
            )

            # Log running executions for debug
            for execution in running_executions[:5]:  # Show only first 5
                logger.debug(
                    f"Execution in progress: {execution.get('unit_id', 'unknown')}"
                )

            time.sleep(check_interval)

        # If we got here, timeout was reached
        remaining_executions = self.work_store.get_running_executions()
        if remaining_executions:
            logger.warning(
                f"Timeout reached waiting for {len(remaining_executions)} executions to complete. "
                f"Some executions may be left in inconsistent state."
            )

            # Force timeout on pending executions
            for execution in remaining_executions:
                try:
                    execution_id = execution.get("id")
                    unit_id = execution.get("unit_id", "unknown")
                    logger.warning(f"Forcing timeout for execution: {unit_id}")

                    # Use execution_scoped to update status
                    execution_store = self.work_store.for_execution(unit_id)
                    execution_store.update_execution_status(
                        status=BaseStatus.FAILED.value,
                        result={
                            "error": "Execution timeout during pipeline finalization"
                        },
                    )
                except Exception as e:
                    logger.error(f"Error forcing timeout on execution {unit_id}: {e}")
        else:
            logger.info("All executions were completed during wait period")
