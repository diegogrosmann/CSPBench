"""Progress broker for centralizing and distributing progress events."""

from typing import Any, Callable, Dict, List, Optional
from threading import Lock
import logging

from .progress_events import (
    ProgressEvent,
    TaskStartedEvent,
    TaskFinishedEvent,
    ExecutionStartedEvent,
    ExecutionProgressEvent,
    ExecutionFinishedEvent,
    AlgorithmProgressEvent,
    ErrorEvent,
    DisplayEvent,
)


class ProgressBroker:
    """
    Centralized broker for handling progress events.

    Receives events from orchestrators and distributes them to registered listeners.
    This decouples orchestrators from display components.
    """

    def __init__(self):
        """Initialize the progress broker."""
        self._listeners: Dict[str, List[Callable[[ProgressEvent], None]]] = {}
        self._global_listeners: List[Callable[[ProgressEvent], None]] = []
        self._lock = Lock()
        self._logger = logging.getLogger(__name__)

    def subscribe(
        self, event_type: str, listener: Callable[[ProgressEvent], None]
    ) -> None:
        """
        Subscribe a listener to a specific event type.

        Args:
            event_type: Type of event to listen for (e.g., 'task_started', 'algorithm_progress')
            listener: Callable that will receive the event
        """
        with self._lock:
            if event_type not in self._listeners:
                self._listeners[event_type] = []
            self._listeners[event_type].append(listener)
            self._logger.debug(f"Subscribed listener to {event_type}")

    def subscribe_all(self, listener: Callable[[ProgressEvent], None]) -> None:
        """
        Subscribe a listener to all event types.

        Args:
            listener: Callable that will receive all events
        """
        with self._lock:
            self._global_listeners.append(listener)
            self._logger.debug("Subscribed global listener")

    def unsubscribe(
        self, event_type: str, listener: Callable[[ProgressEvent], None]
    ) -> None:
        """
        Unsubscribe a listener from a specific event type.

        Args:
            event_type: Type of event to stop listening for
            listener: Listener to remove
        """
        with self._lock:
            if event_type in self._listeners:
                try:
                    self._listeners[event_type].remove(listener)
                    self._logger.debug(f"Unsubscribed listener from {event_type}")
                except ValueError:
                    pass

    def unsubscribe_all(self, listener: Callable[[ProgressEvent], None]) -> None:
        """
        Unsubscribe a listener from all events.

        Args:
            listener: Listener to remove
        """
        with self._lock:
            try:
                self._global_listeners.remove(listener)
                self._logger.debug("Unsubscribed global listener")
            except ValueError:
                pass

    def emit(self, event: ProgressEvent) -> None:
        """
        Emit an event to all relevant listeners.

        Args:
            event: Event to emit
        """
        event_type = type(event).__name__.lower().replace("event", "")

        with self._lock:
            # Send to specific event type listeners
            if event_type in self._listeners:
                for listener in self._listeners[event_type]:
                    try:
                        listener(event)
                    except Exception as e:
                        self._logger.error(
                            f"Error in specific listener for {event_type}: {e}"
                        )

            # Send to global listeners
            for listener in self._global_listeners:
                try:
                    listener(event)
                except Exception as e:
                    self._logger.error(f"Error in global listener: {e}")

    def emit_task_started(
        self,
        task_type,
        task_name: str,
        metadata: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit TaskStartedEvent."""
        event = TaskStartedEvent(
            task_type=task_type,
            task_name=task_name,
            metadata=metadata or {},
            session_id=session_id,
        )
        self.emit(event)

    def emit_task_finished(
        self,
        task_type,
        task_name: str,
        success: bool,
        results: Optional[Dict[str, Any]] = None,
        error_message: Optional[str] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit TaskFinishedEvent."""
        event = TaskFinishedEvent(
            task_type=task_type,
            task_name=task_name,
            success=success,
            results=results or {},
            error_message=error_message,
            session_id=session_id,
        )
        self.emit(event)

    def emit_execution_started(
        self,
        execution_name: str,
        total_items: int,
        metadata: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit ExecutionStartedEvent."""
        event = ExecutionStartedEvent(
            execution_name=execution_name,
            total_items=total_items,
            metadata=metadata or {},
            session_id=session_id,
        )
        self.emit(event)

    def emit_execution_progress(
        self,
        execution_name: str,
        current_item: int,
        total_items: int,
        item_name: str,
        progress_percent: float,
        message: str = "",
        context: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit ExecutionProgressEvent."""
        event = ExecutionProgressEvent(
            execution_name=execution_name,
            current_item=current_item,
            total_items=total_items,
            item_name=item_name,
            progress_percent=progress_percent,
            message=message,
            context=context or {},
            session_id=session_id,
        )
        self.emit(event)

    def emit_execution_finished(
        self,
        execution_name: str,
        success: bool,
        total_processed: int,
        results: Optional[Dict[str, Any]] = None,
        error_message: Optional[str] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit ExecutionFinishedEvent."""
        event = ExecutionFinishedEvent(
            execution_name=execution_name,
            success=success,
            total_processed=total_processed,
            results=results or {},
            error_message=error_message,
            session_id=session_id,
        )
        self.emit(event)

    def emit_algorithm_progress(
        self,
        algorithm_name: str,
        progress_percent: float,
        message: str = "",
        item_id: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit AlgorithmProgressEvent."""
        event = AlgorithmProgressEvent(
            algorithm_name=algorithm_name,
            progress_percent=progress_percent,
            message=message,
            item_id=item_id,
            context=context or {},
            session_id=session_id,
        )
        self.emit(event)

    def emit_algorithm_finished(
        self,
        algorithm_name: str,
        success: bool,
        best_string: str = "",
        max_distance: int = -1,
        execution_time: float = 0.0,
        metadata: Optional[Dict[str, Any]] = None,
        repetition_number: Optional[int] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit AlgorithmFinishedEvent."""
        from .progress_events import AlgorithmFinishedEvent

        event = AlgorithmFinishedEvent(
            algorithm_name=algorithm_name,
            success=success,
            best_string=best_string,
            max_distance=max_distance,
            execution_time=execution_time,
            metadata=metadata or {},
            repetition_number=repetition_number,
            session_id=session_id,
        )
        self.emit(event)

    def emit_warning(
        self,
        algorithm_name: str,
        warning_message: str,
        item_id: Optional[str] = None,
        repetition_number: Optional[int] = None,
        execution_context: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit WarningEvent."""
        from .progress_events import WarningEvent

        event = WarningEvent(
            algorithm_name=algorithm_name,
            warning_message=warning_message,
            item_id=item_id,
            repetition_number=repetition_number,
            execution_context=execution_context or {},
            session_id=session_id,
        )
        self.emit(event)

    def emit_error(
        self,
        error_message: str,
        error_type: str = "generic",
        context: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> None:
        """Convenience method to emit ErrorEvent."""
        event = ErrorEvent(
            error_message=error_message,
            error_type=error_type,
            context=context or {},
            session_id=session_id,
        )
        self.emit(event)

    def emit_display_event(self, event: DisplayEvent) -> None:
        """Emit a unified display event for any phase.

        The display/UI should subscribe to this single stream to render processing,
        optimization and analysis using the same component.
        """
        self.emit(event)

    def clear_listeners(self) -> None:
        """Clear all listeners."""
        with self._lock:
            self._listeners.clear()
            self._global_listeners.clear()
            self._logger.debug("Cleared all listeners")
