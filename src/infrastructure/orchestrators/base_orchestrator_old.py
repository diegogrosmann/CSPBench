"""
Base for Orchestrators

Provides common functionality for all orchestrators,
including standardized integration with the new monitoring system.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

from src.infrastructure.logging_config import get_logger


class BaseOrchestrator(ABC):
    """Abstract base for all system orchestrators."""

    def __init__(self, monitoring_service=None):
        """
        Initialize base orchestrator.

        Args:
            monitoring_service: Optional monitoring service (using new system)
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

    def _emit_execution_started(self, execution_name: str, total_items: int, metadata: Optional[Dict[str, Any]] = None) -> None:
        """
        Emit execution started event.

        Args:
            execution_name: Name of the execution
            total_items: Total items to process
            metadata: Optional metadata
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.start_execution(execution_name, total_items, metadata)
            except Exception as e:
                self._logger.warning(f"Error starting execution in monitoring: {e}")

    def _emit_execution_progress(self, current_item: int, total_items: int, item_name: str, 
                                progress_percent: float, message: str = "", 
                                context: Optional[Dict[str, Any]] = None) -> None:
        """
        Emit execution progress event.

        Args:
            current_item: Current item number
            total_items: Total items
            item_name: Name of current item
            progress_percent: Progress percentage (0-100)
            message: Optional progress message
            context: Optional context data
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.update_execution_progress(
                    current_item=current_item,
                    total_items=total_items,
                    item_name=item_name,
                    progress_percent=progress_percent,
                    message=message,
                    context=context
                )
            except Exception as e:
                self._logger.warning(f"Error updating execution progress in monitoring: {e}")

    def _emit_execution_finished(self, success: bool, total_processed: int, 
                                results: Optional[Dict[str, Any]] = None, 
                                error_message: Optional[str] = None) -> None:
        """
        Emit execution finished event.

        Args:
            success: Whether execution was successful
            total_processed: Total items processed
            results: Optional results data
            error_message: Optional error message if failed
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.finish_execution(success, total_processed, results, error_message)
            except Exception as e:
                self._logger.warning(f"Error finishing execution in monitoring: {e}")

    def _emit_algorithm_progress(self, algorithm_name: str, progress_percent: float, 
                                message: str = "", item_id: Optional[str] = None, 
                                context: Optional[Dict[str, Any]] = None) -> None:
        """
        Emit algorithm progress event.

        Args:
            algorithm_name: Name of the algorithm
            progress_percent: Progress percentage (0-100)
            message: Optional progress message
            item_id: Optional item ID
            context: Optional context data
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.report_algorithm_progress(
                    algorithm_name=algorithm_name,
                    progress_percent=progress_percent,
                    message=message,
                    item_id=item_id,
                    context=context
                )
            except Exception as e:
                self._logger.warning(f"Error reporting algorithm progress in monitoring: {e}")

    def _emit_error(self, error_message: str, error_type: str = "generic", 
                   context: Optional[Dict[str, Any]] = None) -> None:
        """
        Emit error event.

        Args:
            error_message: Error message
            error_type: Type of error
            context: Optional context data
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.report_error(error_message, error_type, context)
            except Exception as e:
                self._logger.warning(f"Error reporting error in monitoring: {e}")

    # Backward compatibility
    def _report_progress(self, progress: float, message: str) -> None:
        """
        Legacy method for backward compatibility.

        Args:
            progress: Progress percentage (0-100)
            message: Progress message
        """
        # Map to new execution progress system if possible
        self._emit_execution_progress(
            current_item=int(progress),
            total_items=100,
            item_name="Processing",
            progress_percent=progress,
            message=message
        )

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

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
                # Use the new monitoring service API
                self.monitoring_service.report_progress(progress, message)
            except Exception as e:
                self._logger.warning(f"Error reporting progress: {e}")

    def _start_execution(self, execution_name: str, total_items: int, metadata: Optional[Dict[str, Any]] = None) -> None:
        """
        Start an execution block.
        
        Args:
            execution_name: Name of the execution
            total_items: Total number of items to process
            metadata: Optional metadata
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.start_execution(execution_name, total_items, metadata)
            except Exception as e:
                self._logger.warning(f"Error starting execution: {e}")

    def _update_execution_progress(self, current_item: int, total_items: int, item_name: str, 
                                  progress_percent: float, message: str = "", 
                                  context: Optional[Dict[str, Any]] = None) -> None:
        """
        Update execution progress.
        
        Args:
            current_item: Current item number
            total_items: Total items
            item_name: Name of current item
            progress_percent: Progress percentage (0-100)
            message: Optional progress message
            context: Optional context data
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.update_execution_progress(
                    current_item, total_items, item_name, progress_percent, message, context
                )
            except Exception as e:
                self._logger.warning(f"Error updating execution progress: {e}")

    def _finish_execution(self, success: bool, total_processed: int, results: Optional[Dict[str, Any]] = None,
                         error_message: Optional[str] = None) -> None:
        """
        Finish the current execution.
        
        Args:
            success: Whether the execution was successful
            total_processed: Total items processed
            results: Optional results data
            error_message: Optional error message if failed
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.finish_execution(success, total_processed, results, error_message)
            except Exception as e:
                self._logger.warning(f"Error finishing execution: {e}")

    def _report_algorithm_progress(self, algorithm_name: str, progress_percent: float, message: str = "",
                                  item_id: Optional[str] = None, context: Optional[Dict[str, Any]] = None) -> None:
        """
        Report algorithm progress.
        
        Args:
            algorithm_name: Name of the algorithm
            progress_percent: Progress percentage (0-100)
            message: Optional progress message
            item_id: Optional item ID
            context: Optional context data
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.report_algorithm_progress(
                    algorithm_name, progress_percent, message, item_id, context
                )
            except Exception as e:
                self._logger.warning(f"Error reporting algorithm progress: {e}")

    def _report_error(self, error_message: str, error_type: str = "generic", 
                     context: Optional[Dict[str, Any]] = None) -> None:
        """
        Report an error.
        
        Args:
            error_message: Error message
            error_type: Type of error
            context: Optional context data
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.report_error(error_message, error_type, context)
            except Exception as e:
                self._logger.warning(f"Error reporting error: {e}")

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
