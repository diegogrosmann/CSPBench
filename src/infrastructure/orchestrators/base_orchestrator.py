"""
Base for Orchestrators

Provides common functionality for all orchestrators,
including standardized integration with monitoring system.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict

from src.infrastructure.logging_config import get_logger


class BaseOrchestrator(ABC):
    """Abstract base for all system orchestrators."""

    def __init__(self, monitoring_service=None):
        """
        Initialize base orchestrator.

        Args:
            monitoring_service: Optional monitoring service
        """
        self.monitoring_service = monitoring_service
        self._logger = get_logger(__name__)

    @abstractmethod
    def execute(self, **kwargs) -> Dict[str, Any]:
        """
        Execute specific orchestration.

        Returns:
            Dict[str, Any]: Orchestration result
        """
        pass

    def _report_progress(self, progress: float, message: str) -> None:
        """
        Report progress to monitoring system.

        Args:
            progress: Progress percentage (0-100)
            message: Descriptive progress message
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.report_progress(progress, message)
            except Exception as e:
                self._logger.warning(f"Error reporting progress: {e}")

    def _update_task_data(self, task_type: str, **data) -> None:
        """
        Update task-specific data in monitoring.

        Args:
            task_type: Task type (execution, optimization, sensitivity)
            **data: Task-specific data
        """
        if not self.monitoring_service:
            return

        try:
            if task_type == "execution":
                self.monitoring_service.update_execution_data(**data)
            elif task_type == "optimization":
                self.monitoring_service.update_optimization_data(**data)
            elif task_type == "sensitivity":
                self.monitoring_service.update_sensitivity_data(**data)
            else:
                self._logger.warning(f"Unknown task type: {task_type}")
        except Exception as e:
            self._logger.warning(f"Error updating task data: {e}")

    def _log_error(self, error: Exception, context: str = "") -> None:
        """
        Log error with appropriate context.

        Args:
            error: Exception that occurred
            context: Additional error context
        """
        error_msg = "Error in orchestration"
        if context:
            error_msg += f" ({context})"
        error_msg += f": {error}"

        self._logger.error(error_msg)

    def _log_info(self, message: str) -> None:
        """
        Log information message.

        Args:
            message: Message to log
        """
        self._logger.info(message)

    def _log_debug(self, message: str) -> None:
        """
        Log debug information message.

        Args:
            message: Message to log
        """
        self._logger.debug(message)
