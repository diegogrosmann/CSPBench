"""Simplified monitoring service using progress broker."""

from typing import Any, Dict, Optional
import logging

from .progress_broker import ProgressBroker
from .progress_events import TaskType, DisplayEvent


class MonitoringService:
    """
    Simplified monitoring service that uses the progress broker.

    This service is injected into orchestrators and provides a simple API
    for reporting progress without knowing about display concerns.
    """

    def __init__(
        self, progress_broker: ProgressBroker, session_id: Optional[str] = None
    ):
        """
        Initialize monitoring service.

        Args:
            progress_broker: Progress broker for event distribution
            session_id: Optional session ID for tracking
        """
        self._broker = progress_broker
        self._session_id = session_id
        self._logger = logging.getLogger(__name__)
        self._active_task: Optional[str] = None
        self._active_execution: Optional[str] = None

    def start_task(
        self,
        task_type: TaskType,
        task_name: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Start a new task.

        Args:
            task_type: Type of task (EXECUTION, OPTIMIZATION, SENSITIVITY)
            task_name: Name of the task
            metadata: Optional metadata
        """
        self._active_task = task_name
        self._broker.emit_task_started(
            task_type=task_type,
            task_name=task_name,
            metadata=metadata,
            session_id=self._session_id,
        )
        self._logger.info(f"Task started: {task_name} ({task_type.value})")

    def notify_execution_started(
        self, execution_name: str, metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Notify that a new execution has started.

        Args:
            execution_name: Name of the execution
            metadata: Optional metadata
        """
        self._active_execution = execution_name
        self._broker.emit_execution_started(
            execution_name=execution_name,
            total_items=1,  # Default for now
            metadata=metadata,
        )
        self._logger.info(f"Execution started: {execution_name}")

    def finish_task(
        self,
        success: bool,
        results: Optional[Dict[str, Any]] = None,
        error_message: Optional[str] = None,
    ) -> None:
        """
        Finish the current task.

        Args:
            success: Whether the task was successful
            results: Optional results data
            error_message: Optional error message if failed
        """
        if not self._active_task:
            self._logger.warning("Trying to finish task but no active task")
            return

        self._broker.emit_task_finished(
            task_type=TaskType.EXECUTION,  # We can infer this or store it
            task_name=self._active_task,
            success=success,
            results=results,
            error_message=error_message,
            session_id=self._session_id,
        )
        self._logger.info(f"Task finished: {self._active_task} (success: {success})")
        self._active_task = None

    def start_execution(
        self,
        execution_name: str,
        total_items: int,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Start an execution block.

        Args:
            execution_name: Name of the execution
            total_items: Total number of items to process
            metadata: Optional metadata
        """
        self._active_execution = execution_name
        self._broker.emit_execution_started(
            execution_name=execution_name,
            total_items=total_items,
            metadata=metadata,
            session_id=self._session_id,
        )
        self._logger.debug(f"Execution started: {execution_name} ({total_items} items)")

    def update_execution_progress(
        self,
        current_item: int,
        total_items: int,
        item_name: str,
        progress_percent: float,
        message: str = "",
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
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
        if not self._active_execution:
            self._logger.warning(
                "Trying to update execution progress but no active execution"
            )
            return

        self._broker.emit_execution_progress(
            execution_name=self._active_execution,
            current_item=current_item,
            total_items=total_items,
            item_name=item_name,
            progress_percent=progress_percent,
            message=message,
            context=context,
            session_id=self._session_id,
        )

    def finish_execution(
        self,
        success: bool,
        total_processed: int,
        results: Optional[Dict[str, Any]] = None,
        error_message: Optional[str] = None,
    ) -> None:
        """
        Finish the current execution.

        Args:
            success: Whether the execution was successful
            total_processed: Total items processed
            results: Optional results data
            error_message: Optional error message if failed
        """
        if not self._active_execution:
            self._logger.warning("Trying to finish execution but no active execution")
            return

        self._broker.emit_execution_finished(
            execution_name=self._active_execution,
            success=success,
            total_processed=total_processed,
            results=results,
            error_message=error_message,
            session_id=self._session_id,
        )
        self._logger.debug(
            f"Execution finished: {self._active_execution} (success: {success})"
        )
        self._active_execution = None

    def report_algorithm_progress(
        self,
        algorithm_name: str,
        progress_percent: float,
        message: str = "",
        item_id: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Report algorithm progress.

        Args:
            algorithm_name: Name of the algorithm
            progress_percent: Progress percentage (0-100)
            message: Optional progress message
            item_id: Optional item ID
            context: Optional context data
        """
        self._broker.emit_algorithm_progress(
            algorithm_name=algorithm_name,
            progress_percent=progress_percent,
            message=message,
            item_id=item_id,
            context=context,
            session_id=self._session_id,
        )

    def report_algorithm_finished(
        self,
        algorithm_name: str,
        success: bool,
        best_string: str = "",
        max_distance: int = -1,
        execution_time: float = 0.0,
        metadata: Optional[Dict[str, Any]] = None,
        repetition_number: Optional[int] = None,
    ) -> None:
        """
        Report algorithm completion.

        Args:
            algorithm_name: Name of the algorithm
            success: Whether the algorithm completed successfully
            best_string: Best solution found
            max_distance: Maximum distance achieved
            execution_time: Time taken for execution
            metadata: Optional metadata from algorithm
            repetition_number: Optional repetition number for individual tracking
        """
        self._broker.emit_algorithm_finished(
            algorithm_name=algorithm_name,
            success=success,
            best_string=best_string,
            max_distance=max_distance,
            execution_time=execution_time,
            metadata=metadata,
            repetition_number=repetition_number,
            session_id=self._session_id,
        )
        self._logger.debug(
            f"Algorithm finished: {algorithm_name} (success: {success})"
        )

    def report_warning(
        self,
        algorithm_name: str,
        warning_message: str,
        item_id: Optional[str] = None,
        repetition_number: Optional[int] = None,
        execution_context: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Report an algorithm warning.

        Args:
            algorithm_name: Name of the algorithm
            warning_message: Warning message
            item_id: Optional item ID
            repetition_number: Optional repetition number for individual tracking
            execution_context: Optional execution context data
        """
        self._broker.emit_warning(
            algorithm_name=algorithm_name,
            warning_message=warning_message,
            item_id=item_id,
            repetition_number=repetition_number,
            execution_context=execution_context or {},
            session_id=self._session_id,
        )
        self._logger.warning(f"Algorithm warning ({algorithm_name}): {warning_message}")

    def report_error(
        self,
        error_message: str,
        error_type: str = "generic",
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Report an error.

        Args:
            error_message: Error message
            error_type: Type of error
            context: Optional context data
        """
        self._broker.emit_error(
            error_message=error_message,
            error_type=error_type,
            context=context,
            session_id=self._session_id,
        )
        self._logger.error(f"Error reported: {error_message}")

    # Unified display API
    def emit_event(self, event: DisplayEvent) -> None:
        """Publish a unified display event to the broker/UI.

        Keeps a single callback for processing, optimization and analysis.
        """
        try:
            self._broker.emit_display_event(event)
        except Exception as e:
            self._logger.warning(f"Failed to emit display event: {e}")

    # Compatibility methods for monitor interface
    def update_hierarchy(
        self, level, level_id: str, progress: float, message: str = "", data=None
    ) -> None:
        """Compatibility method - maps to new update_execution_progress."""
        if self._active_execution:
            self.update_execution_progress(
                current_item=int(progress),
                total_items=100,
                item_name=level_id,
                progress_percent=progress,
                message=message,
                context=data,
            )

    def finish_monitoring(self, results: Optional[Dict[str, Any]] = None) -> None:
        """Compatibility method - maps to finish_task."""
        success = bool(results and results.get("results"))
        self.finish_task(success, results)

    def start_item(
        self, item_id: str, item_type: str = "repetition", context=None, metadata=None
    ) -> None:
        """Compatibility method for starting items."""
        # No-op for now, just log
        self._logger.debug(f"Started item: {item_id} ({item_type})")

    def finish_item(
        self,
        item_id: str,
        success: bool = True,
        result=None,
        error: Optional[str] = None,
    ) -> None:
        """Compatibility method for finishing items."""
        # No-op for now, just log
        if not success and error:
            self.report_error(error)
        self._logger.debug(f"Finished item: {item_id} (success: {success})")

    def update_item(
        self, item_id: str, progress: float, message: str = "", context=None
    ) -> None:
        """Compatibility method for updating items."""
        # Map to algorithm progress if possible
        self.report_algorithm_progress(
            algorithm_name=item_id,
            progress_percent=progress,
            message=message,
            context=context,
        )

    def algorithm_callback(
        self,
        algorithm_name: str,
        progress: float,
        message: str = "",
        item_id: Optional[str] = None,
    ) -> None:
        """Compatibility method for algorithm callbacks."""
        self.report_algorithm_progress(algorithm_name, progress, message, item_id)
