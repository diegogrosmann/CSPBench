"""Web display adapter for progress events."""

from typing import Optional, Any, Dict, TYPE_CHECKING
import logging

if TYPE_CHECKING:
    from src.application.monitoring.progress_events import ProgressEvent, AlgorithmProgressEvent, ExecutionProgressEvent


class WebDisplay:
    """
    Display adapter for web interface.
    
    Listens to progress events and updates the web session manager
    with structured progress data for the web interface.
    """
    
    def __init__(self, session_manager, session_id: str):
        """
        Initialize web display.
        
        Args:
            session_manager: Web session manager
            session_id: Session ID to update
        """
        self._session_manager = session_manager
        self._session_id = session_id
        self._logger = logging.getLogger(__name__)
        
        # Track current state
        self._current_task = None
        self._current_execution = None
        self._execution_progress: Dict[str, Any] = {}
        
    def handle_event(self, event: Any) -> None:
        """
        Handle a progress event.
        
        Args:
            event: Progress event to handle
        """
        try:
            event_type = type(event).__name__
            
            if hasattr(self, f'_handle_{event_type.lower().replace("event", "")}'):
                handler = getattr(self, f'_handle_{event_type.lower().replace("event", "")}')
                handler(event)
            else:
                self._handle_generic(event)
                
        except Exception as e:
            self._logger.error(f"Error handling web event {type(event).__name__}: {e}")
    
    def _handle_taskstarted(self, event) -> None:
        """Handle TaskStartedEvent."""
        self._current_task = {
            "type": event.task_type.value,
            "name": event.task_name,
            "status": "running",
            "started_at": event.timestamp.isoformat() if event.timestamp else None,
            "metadata": event.metadata
        }
        
        self._update_session({
            "status": "running",
            "message": f"Starting {event.task_type.value}: {event.task_name}",
            "current_task": self._current_task
        })
        
        self._add_log("INFO", f"Starting {event.task_type.value}: {event.task_name}")
    
    def _handle_taskfinished(self, event) -> None:
        """Handle TaskFinishedEvent."""
        if self._current_task:
            self._current_task.update({
                "status": "completed" if event.success else "failed",
                "finished_at": event.timestamp.isoformat() if event.timestamp else None,
                "results": event.results
            })
        
        status = "completed" if event.success else "failed"
        message = f"Task {status}: {event.task_name}"
        if event.error_message:
            message += f" - {event.error_message}"
        
        self._update_session({
            "status": status,
            "message": message,
            "current_task": self._current_task,
            "completed_at": event.timestamp.isoformat() if event.timestamp else None
        })
        
        log_level = "SUCCESS" if event.success else "ERROR"
        self._add_log(log_level, message)
    
    def _handle_executionstarted(self, event) -> None:
        """Handle ExecutionStartedEvent."""
        self._current_execution = {
            "name": event.execution_name,
            "total_items": event.total_items,
            "current_item": 0,
            "progress": 0.0,
            "status": "running",
            "started_at": event.timestamp.isoformat() if event.timestamp else None
        }
        
        self._execution_progress = {
            "current_execution": 1,
            "total_executions": 1,  # This should be updated from higher level context
            "current_dataset": 1,
            "total_datasets": 1,
            "current_algorithm": 1,
            "total_algorithms": 1,
            "current_run": 1,
            "total_runs": 1,
            "overall_progress": 0
        }
        
        self._update_session({
            "current_execution": self._current_execution,
            "progress": self._execution_progress
        })
        
        self._add_log("INFO", f"Execution started: {event.execution_name} ({event.total_items} items)")
    
    def _handle_executionprogress(self, event) -> None:
        """Handle ExecutionProgressEvent."""
        if self._current_execution:
            self._current_execution.update({
                "current_item": event.current_item,
                "item_name": event.item_name,
                "progress": event.progress_percent,
                "message": event.message
            })
        
        # Update hierarchical progress if context is available
        if event.context:
            self._execution_progress.update(self._map_context_to_progress(event.context))
        
        self._execution_progress["overall_progress"] = int(event.progress_percent)
        
        self._update_session({
            "current_execution": self._current_execution,
            "progress": self._execution_progress
        })
    
    def _handle_executionfinished(self, event) -> None:
        """Handle ExecutionFinishedEvent."""
        if self._current_execution:
            self._current_execution.update({
                "status": "completed" if event.success else "failed",
                "total_processed": event.total_processed,
                "finished_at": event.timestamp.isoformat() if event.timestamp else None
            })
        
        message = f"Execution {event.execution_name}: "
        message += "completed" if event.success else "failed"
        message += f" ({event.total_processed} items)"
        
        self._update_session({
            "current_execution": self._current_execution
        })
        
        log_level = "SUCCESS" if event.success else "ERROR"
        self._add_log(log_level, message)
    
    def _handle_algorithmprogress(self, event) -> None:
        """Handle AlgorithmProgressEvent."""
        # Update current execution with algorithm details
        if self._current_execution:
            self._current_execution.update({
                "algorithm": event.algorithm_name,
                "algorithm_progress": event.progress_percent,
                "algorithm_message": event.message
            })
        
        self._update_session({
            "current_execution": self._current_execution
        })
        
        # Send algorithm log
        self._add_log("INFO", f"[{event.algorithm_name}] {event.message} ({event.progress_percent:.1f}%)",
                     source=event.algorithm_name)
    
    def _handle_error(self, event) -> None:
        """Handle ErrorEvent."""
        self._add_log("ERROR", f"Error ({event.error_type}): {event.error_message}")
        
        self._update_session({
            "status": "error",
            "message": event.error_message
        })
    
    def _handle_generic(self, event) -> None:
        """Handle generic events."""
        self._add_log("DEBUG", f"{type(event).__name__}: {getattr(event, 'message', 'Event occurred')}")
    
    def _map_context_to_progress(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Map context data to hierarchical progress structure.
        
        Args:
            context: Context data from event
            
        Returns:
            Mapped progress data
        """
        progress_data = {}
        
        # Map from internal data if available
        if context:
            # Execution level
            progress_data["current_execution"] = context.get("execution_index", context.get("config_index", 1))
            progress_data["total_executions"] = context.get("total_executions", context.get("total_configs", 1))
            
            # Dataset level
            progress_data["current_dataset"] = context.get("dataset_index", 1)
            progress_data["total_datasets"] = context.get("total_datasets", 1)
            
            # Algorithm level
            progress_data["current_algorithm"] = context.get("algorithm_index", context.get("algorithm_config_index", 1))
            progress_data["total_algorithms"] = context.get("total_algorithms", 1)
            
            # Run level
            progress_data["current_run"] = context.get("current_repetition", context.get("current_run", 1))
            progress_data["total_runs"] = context.get("total_repetitions", context.get("total_runs", 1))
        
        return progress_data
    
    def _update_session(self, data: Dict[str, Any]) -> None:
        """
        Update session with data.
        
        Args:
            data: Data to update in session
        """
        try:
            self._session_manager.update_session(self._session_id, data)
        except Exception as e:
            self._logger.error(f"Error updating web session: {e}")
    
    def _add_log(self, level: str, message: str, source: Optional[str] = None) -> None:
        """
        Add log to session.
        
        Args:
            level: Log level
            message: Log message
            source: Optional source of the log
        """
        try:
            self._session_manager.add_log(self._session_id, level, message, source=source)
        except Exception as e:
            self._logger.error(f"Error adding log to web session: {e}")
