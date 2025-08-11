"""
Standardized Interface for Executors

Defines the common contract that all executors must implement,
ensuring consistency in execution of different types of tasks.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, List

from src.domain import Dataset


class ExecutorInterface(ABC):
    """Standardized interface for CSP algorithm executors."""

    @abstractmethod
    def execute_batch(
        self,
        batch_config: Dict[str, Any],
        monitoring_service=None,
        session_manager=None,
    ) -> List[Dict[str, Any]]:
        """
        Execute batch of algorithms (including single executions).

        Args:
            batch_config: Batch configuration
            monitoring_service: Optional monitoring service

        Returns:
            List[Dict[str, Any]]: List of execution results
        """
        pass

    @abstractmethod
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
    ) -> Dict[str, Any]:
        """
        Execute hyperparameter optimization.

        Args:
            algorithm_name: Algorithm name
            dataset: Dataset for optimization
            optimization_config: Optimization configuration
            monitoring_service: Optional monitoring service
            config_index: Current configuration index
            total_configs: Total configurations
            dataset_index: Current dataset index
            total_datasets: Total datasets

        Returns:
            Dict[str, Any]: Optimization result containing:
                - best_params: best parameters found
                - best_value: best value achieved
                - n_trials: number of trials executed
                - optimization_history: optimization history
        """
        pass

    @abstractmethod
    def execute_sensitivity_analysis(
        self,
        algorithm_name: str,
        dataset: Dataset,
        sensitivity_config: Dict[str, Any],
        monitoring_service=None,
        *,
        task_index: int = 1,
        total_tasks: int = 1,
        dataset_index: int = 1,
        total_datasets: int = 1,
        config_index: int = 1,
        total_configs: int = 1,
        algorithm_index: int = 1,
        total_algorithms: int = 1,
    ) -> Dict[str, Any]:
        """
        Execute sensitivity analysis.

        Args:
            algorithm_name: Algorithm name
            dataset: Dataset for analysis
            sensitivity_config: Analysis configuration
            monitoring_service: Optional monitoring service

        Returns:
            Dict[str, Any]: Analysis result containing:
                - sensitivity_indices: sensitivity indices
                - parameter_rankings: parameter rankings
                - analysis_method: method used
                - samples_executed: number of samples executed
        """
        pass

    @abstractmethod
    def get_execution_status(self, execution_id: str) -> str:
        """
        Get status of a specific execution.

        Args:
            execution_id: Execution ID

        Returns:
            str: Execution status (running, completed, failed, cancelled)
        """
        pass

    @abstractmethod
    def cancel_execution(self, execution_id: str) -> bool:
        """
        Cancel a running execution.

        Args:
            execution_id: Execution ID to cancel

        Returns:
            bool: True if cancellation was successful
        """
        pass

    def set_batch_config(self, batch_config: Dict[str, Any]) -> None:
        """
        Set current batch configuration (optional method).

        Args:
            batch_config: Batch configuration
        """
        pass
