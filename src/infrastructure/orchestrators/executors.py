"""
Main Executor - Task Router

Implements standardized interface delegating execution to specialized orchestrators.
"""

import time
from typing import Any, Dict, List

from src.application.ports import ExecutorInterface
from src.domain import Dataset
from src.domain.errors import AlgorithmExecutionError
from src.infrastructure.logging_config import get_logger


class Executor(ExecutorInterface):
    """Main executor that delegates tasks to specialized orchestrators."""

    def __init__(self):
        """Initialize executor as pure router."""
        self._logger = get_logger(__name__)
        self._current_batch_config = None

    def set_batch_config(self, batch_config: Dict[str, Any]) -> None:
        """Set current batch configuration (delegated to orchestrators)."""
        self._current_batch_config = batch_config
        self._logger.debug(f"Batch configuration set: {type(batch_config)}")

    def execute_batch(
        self, batch_config: Dict[str, Any], monitoring_service=None
    ) -> List[Dict[str, Any]]:
        """
        Execute algorithm batch delegating to ExecutionOrchestrator.

        Args:
            batch_config: Batch configuration
            monitoring_service: Optional monitoring service

        Returns:
            List[Dict[str, Any]]: List of execution results
        """
        try:
            from src.infrastructure import (
                DomainAlgorithmRegistry,
                FileDatasetRepository,
            )
            from src.infrastructure.orchestrators.execution_orchestrator import (
                ExecutionOrchestrator,
            )

            # Configure dependencies
            algorithm_registry = DomainAlgorithmRegistry()
            dataset_repository = FileDatasetRepository("./datasets")

            # Create execution orchestrator
            orchestrator = ExecutionOrchestrator(
                algorithm_registry=algorithm_registry,
                dataset_repository=dataset_repository,
                monitoring_service=monitoring_service,
            )

            # Delegate execution
            return orchestrator.execute_batch(batch_config, monitoring_service)

        except Exception as e:
            self._logger.error(f"Error in batch execution: {e}")
            raise AlgorithmExecutionError(f"Error in batch execution: {e}") from e

    def execute_optimization(
        self,
        algorithm_name: str,
        dataset: Dataset,
        optimization_config: Dict[str, Any],
        monitoring_service=None,
        config_index: int = 1,
        total_configs: int = 1,
        dataset_index: int = 1,
        total_datasets: int = 1,
        dataset_name: str | None = None,
    ) -> Dict[str, Any]:
        """
        Execute hyperparameter optimization delegating to OptimizationOrchestrator.

        Args:
            algorithm_name: Algorithm name
            dataset: Dataset for optimization
            optimization_config: Optimization configuration
            monitoring_service: Optional monitoring service
            config_index: Current configuration index
            total_configs: Total configurations
            dataset_index: Current dataset index
            total_datasets: Total datasets
            dataset_name: Original dataset name (optional)

        Returns:
            Dict[str, Any]: Optimization result
        """
        try:
            from src.infrastructure import (
                DomainAlgorithmRegistry,
                FileDatasetRepository,
            )
            from src.infrastructure.orchestrators.optimization_orchestrator import (
                OptimizationOrchestrator,
            )

            # Configure dependencies
            algorithm_registry = DomainAlgorithmRegistry()
            dataset_repository = FileDatasetRepository("./datasets")

            # Use original dataset name if provided, otherwise use temporary
            if dataset_name:
                dataset_id = dataset_name
            else:
                dataset_id = f"temp_optimization_{int(time.time())}"

            # Save dataset temporarily in repository for orchestrator
            dataset_repository.save(dataset, dataset_id)

            # Create complete configuration for orchestrator
            full_config = {
                "algorithm": algorithm_name,
                "dataset": dataset_id,
                "optimization": optimization_config,
                "export": optimization_config.get("export", {"enabled": True}),
                "monitoring": optimization_config.get("monitoring", {}),
                "resources": optimization_config.get("resources", {}),
                "plots": optimization_config.get("plots", {}),
            }

            # Create orchestrator
            orchestrator = OptimizationOrchestrator(
                algorithm_registry=algorithm_registry,
                dataset_repository=dataset_repository,
                config=full_config,
                monitoring_service=monitoring_service,
                config_index=config_index,
                total_configs=total_configs,
                dataset_index=dataset_index,
                total_datasets=total_datasets,
            )

            # Execute optimization
            results = orchestrator.run_optimization()

            # Clean up temporary dataset only if it was temporary
            if not dataset_name:
                try:
                    dataset_repository.delete(dataset_id)
                except:
                    pass

            return results

        except Exception as e:
            self._logger.error(f"Error in optimization: {e}")
            raise AlgorithmExecutionError(f"Error in optimization: {e}") from e

    def execute_sensitivity_analysis(
        self,
        algorithm_name: str,
        dataset: Dataset,
        sensitivity_config: Dict[str, Any],
        monitoring_service=None,
    ) -> Dict[str, Any]:
        """
        Execute sensitivity analysis delegating to SensitivityOrchestrator.

        Args:
            algorithm_name: Algorithm name
            dataset: Dataset for analysis
            sensitivity_config: Analysis configuration
            monitoring_service: Optional monitoring service

        Returns:
            Dict[str, Any]: Sensitivity analysis result
        """
        try:
            from src.infrastructure import DomainAlgorithmRegistry
            from src.infrastructure.orchestrators.sensitivity_orchestrator import (
                SensitivityOrchestrator,
            )

            # Configure registry
            algorithm_registry = DomainAlgorithmRegistry()

            # Create sensitivity orchestrator
            orchestrator = SensitivityOrchestrator(
                algorithm_registry, self, monitoring_service=monitoring_service
            )

            # Execute analysis
            results = orchestrator.execute_sensitivity_analysis(
                algorithm_name, dataset, sensitivity_config
            )

            return results

        except Exception as e:
            self._logger.error(f"Error in sensitivity analysis: {e}")
            raise AlgorithmExecutionError(
                f"Error in sensitivity analysis: {e}"
            ) from e

    def get_execution_status(self, execution_id: str) -> str:
        """
        Get status of a specific execution.

        Args:
            execution_id: Execution ID

        Returns:
            str: Execution status (running, completed, failed, cancelled)
        """
        # For now, return a generic status
        # TODO: Implement execution tracking system
        return "completed"

    def cancel_execution(self, execution_id: str) -> bool:
        """
        Cancel a running execution.

        Args:
            execution_id: ID of execution to cancel

        Returns:
            bool: True if cancellation was successful
        """
        # For now, return False
        # TODO: Implement cancellation system
        return False
