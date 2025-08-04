"""Simplified monitoring service using progress broker."""

from typing import Any, Dict, Optional
import logging

from .progress_broker import ProgressBroker
from .progress_events import TaskType


class MonitoringService:
    """
    Simplified monitoring service that uses the progress broker.
    
    This service is injected into orchestrators and provides a simple API
    for reporting progress without knowing about display concerns.
    """
    
    def __init__(self, progress_broker: ProgressBroker, session_id: Optional[str] = None):
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
    
    def start_task(self, task_type: TaskType, task_name: str, metadata: Optional[Dict[str, Any]] = None) -> None:
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
            session_id=self._session_id
        )
        self._logger.info(f"Task started: {task_name} ({task_type.value})")
    
    def notify_execution_started(self, execution_name: str, metadata: Optional[Dict[str, Any]] = None) -> None:
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
            metadata=metadata
        )
        self._logger.info(f"Execution started: {execution_name}")
    
    def finish_task(self, success: bool, results: Optional[Dict[str, Any]] = None, 
                   error_message: Optional[str] = None) -> None:
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
            session_id=self._session_id
        )
        self._logger.info(f"Task finished: {self._active_task} (success: {success})")
        self._active_task = None
    
    def start_execution(self, execution_name: str, total_items: int, metadata: Optional[Dict[str, Any]] = None) -> None:
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
            session_id=self._session_id
        )
        self._logger.debug(f"Execution started: {execution_name} ({total_items} items)")
    
    def update_execution_progress(self, current_item: int, total_items: int, item_name: str, 
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
        if not self._active_execution:
            self._logger.warning("Trying to update execution progress but no active execution")
            return
            
        self._broker.emit_execution_progress(
            execution_name=self._active_execution,
            current_item=current_item,
            total_items=total_items,
            item_name=item_name,
            progress_percent=progress_percent,
            message=message,
            context=context,
            session_id=self._session_id
        )
    
    def finish_execution(self, success: bool, total_processed: int, results: Optional[Dict[str, Any]] = None,
                        error_message: Optional[str] = None) -> None:
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
            session_id=self._session_id
        )
        self._logger.debug(f"Execution finished: {self._active_execution} (success: {success})")
        self._active_execution = None
    
    def report_algorithm_progress(self, algorithm_name: str, progress_percent: float, message: str = "",
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
        self._broker.emit_algorithm_progress(
            algorithm_name=algorithm_name,
            progress_percent=progress_percent,
            message=message,
            item_id=item_id,
            context=context,
            session_id=self._session_id
        )
    
    def report_error(self, error_message: str, error_type: str = "generic", 
                    context: Optional[Dict[str, Any]] = None) -> None:
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
            session_id=self._session_id
        )
        self._logger.error(f"Error reported: {error_message}")
    
    # Legacy compatibility methods (to be removed gradually)
    def report_progress(self, progress: float, message: str = "") -> None:
        """Legacy method for backward compatibility."""
        if self._active_execution:
            self.update_execution_progress(
                current_item=int(progress),
                total_items=100,
                item_name="Progress",
                progress_percent=progress,
                message=message
            )
    
    def algorithm_callback(self, algorithm_name: str, progress: float, message: str = "", 
                          item_id: Optional[str] = None) -> None:
        """Legacy method for backward compatibility."""
        self.report_algorithm_progress(algorithm_name, progress, message, item_id)
    
    def show_error(self, error_message: str) -> None:
        """Legacy method for backward compatibility."""
        self.report_error(error_message)
    
    def start_monitoring(self, task_type: TaskType, task_name: str, batch_config: Optional[Dict[str, Any]] = None) -> None:
        """Legacy method for backward compatibility."""
        self.start_task(task_type, task_name, batch_config)
    
    def finish_monitoring(self, results: Optional[Dict[str, Any]] = None) -> None:
        """Legacy method for backward compatibility."""
        success = bool(results and results.get("results"))
        self.finish_task(success, results)
    
    def update_hierarchy(self, level, level_id: str, progress: float, message: str = "", data=None) -> None:
        """Legacy method for backward compatibility."""
        # Map old hierarchy calls to new execution progress
        if self._active_execution:
            self.update_execution_progress(
                current_item=int(progress),
                total_items=100,
                item_name=level_id,
                progress_percent=progress,
                message=message,
                context=data
            )
    
    def start_item(self, item_id: str, item_type: str = "repetition", context=None, metadata=None) -> None:
        """Legacy method for backward compatibility."""
        # No-op for now, just log
        self._logger.debug(f"Started item: {item_id} ({item_type})")
    
    def finish_item(self, item_id: str, success: bool = True, result=None, error: Optional[str] = None) -> None:
        """Legacy method for backward compatibility."""
        # No-op for now, just log
        if not success and error:
            self.report_error(error)
        self._logger.debug(f"Finished item: {item_id} (success: {success})")
    
    def update_item(self, item_id: str, progress: float, message: str = "", context=None) -> None:
        """Legacy method for backward compatibility."""
        # Map to algorithm progress if possible
        self.report_algorithm_progress(
            algorithm_name=item_id,
            progress_percent=progress,
            message=message,
            context=context
        )
