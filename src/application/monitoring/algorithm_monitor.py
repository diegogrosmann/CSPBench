"""
Algorithm Monitor - Bridge between CSP algorithms and monitoring system.

Provides integration between algorithm callbacks and the centralized monitoring system,
ensuring all algorithm execution details are captured and propagated as events.
"""

import time
import uuid
from typing import Any, Callable, Dict, Optional
import logging

from .progress_broker import ProgressBroker
from .progress_events import (
    AlgorithmCallbackEvent,
    AlgorithmWarningEvent,
    AlgorithmHistoryEvent,
    RepetitionStartedEvent,
    RepetitionProgressEvent,
    RepetitionFinishedEvent,
    ExecutionLevel,
    ExecutionPhase,
)

logger = logging.getLogger(__name__)


class AlgorithmMonitor:
    """
    Bridge between algorithm callbacks and monitoring system.

    Creates monitored callbacks that capture all algorithm execution details
    and propagate them as structured events through the monitoring system.
    """

    def __init__(
        self,
        broker: ProgressBroker,
        session_id: str,
        algorithm_name: str,
        item_id: Optional[str] = None,
        repetition_id: Optional[str] = None,
    ):
        """
        Initialize algorithm monitor.

        Args:
            broker: Progress broker for event emission
            session_id: Current session ID
            algorithm_name: Name of the algorithm being monitored
            item_id: Optional item ID being processed
            repetition_id: Optional repetition/trial/sample ID
        """
        self.broker = broker
        self.session_id = session_id
        self.algorithm_name = algorithm_name
        self.item_id = item_id
        self.repetition_id = repetition_id or str(uuid.uuid4())

        # Tracking state
        self.start_time = time.time()
        self.last_progress = 0.0
        self.history_entries = []

        logger.debug(
            f"Created AlgorithmMonitor for {algorithm_name} "
            f"(item: {item_id}, repetition: {repetition_id})"
        )

    def create_progress_callback(self) -> Callable[[str, float], None]:
        """
        Create monitored progress callback for algorithm.

        Returns:
            Callback function that captures progress and emits events
        """

        def progress_callback(message: str, progress: float = 0.0) -> None:
            try:
                # Emit algorithm callback event
                event = AlgorithmCallbackEvent(
                    session_id=self.session_id,
                    algorithm_name=self.algorithm_name,
                    callback_type="progress",
                    message=message,
                    progress_percent=progress,
                    item_id=self.item_id,
                    repetition_id=self.repetition_id,
                    context={
                        "elapsed_time": time.time() - self.start_time,
                        "progress_delta": progress - self.last_progress,
                    },
                )
                self.broker.emit(event)

                # Update repetition progress if significant change
                if abs(progress - self.last_progress) >= 1.0:  # 1% threshold
                    repetition_event = RepetitionProgressEvent(
                        session_id=self.session_id,
                        execution_level=ExecutionLevel.REPETITION,
                        level_id=self.repetition_id,
                        item_name=self.item_id or "unknown",
                        algorithm_name=self.algorithm_name,
                        repetition_number=1,  # Can be set externally
                        progress_percent=progress,
                        message=message,
                        phase=ExecutionPhase.RUNNING,
                        metadata={
                            "callback_source": "algorithm_progress",
                            "elapsed_time": time.time() - self.start_time,
                        },
                    )
                    self.broker.emit(repetition_event)
                    self.last_progress = progress

                logger.debug(f"{self.algorithm_name}: {message} ({progress:.1f}%)")

            except Exception as e:
                logger.error(
                    f"Error in progress callback for {self.algorithm_name}: {e}"
                )

        return progress_callback

    def create_warning_callback(self) -> Callable[[str], None]:
        """
        Create monitored warning callback for algorithm.

        Returns:
            Callback function that captures warnings and emits events
        """

        def warning_callback(message: str) -> None:
            try:
                # Emit algorithm warning event
                event = AlgorithmWarningEvent(
                    session_id=self.session_id,
                    algorithm_name=self.algorithm_name,
                    warning_message=message,
                    item_id=self.item_id,
                    repetition_id=self.repetition_id,
                    context={
                        "elapsed_time": time.time() - self.start_time,
                        "warning_count": len(
                            [
                                e
                                for e in self.history_entries
                                if e.get("type") == "warning"
                            ]
                        )
                        + 1,
                    },
                )
                self.broker.emit(event)

                logger.warning(f"{self.algorithm_name}: {message}")

            except Exception as e:
                logger.error(
                    f"Error in warning callback for {self.algorithm_name}: {e}"
                )

        return warning_callback

    def report_history_entry(self, iteration: int, **data) -> None:
        """
        Report algorithm history entry.

        Args:
            iteration: Current iteration number
            **data: History data (fitness, best solution, etc.)
        """
        try:
            # Store locally
            entry = {
                "iteration": iteration,
                "timestamp": time.time(),
                "elapsed_time": time.time() - self.start_time,
                **data,
            }
            self.history_entries.append(entry)

            # Emit history event
            event = AlgorithmHistoryEvent(
                session_id=self.session_id,
                algorithm_name=self.algorithm_name,
                iteration=iteration,
                history_data=entry,
                item_id=self.item_id,
                repetition_id=self.repetition_id,
            )
            self.broker.emit(event)

            logger.debug(
                f"{self.algorithm_name}: History entry for iteration {iteration}"
            )

        except Exception as e:
            logger.error(f"Error reporting history for {self.algorithm_name}: {e}")

    def start_repetition(
        self, repetition_number: int = 1, total_repetitions: int = 1
    ) -> None:
        """
        Signal start of repetition/trial/sample.

        Args:
            repetition_number: Current repetition number
            total_repetitions: Total number of repetitions
        """
        try:
            event = RepetitionStartedEvent(
                session_id=self.session_id,
                execution_level=ExecutionLevel.REPETITION,
                level_id=self.repetition_id,
                item_name=self.item_id or "unknown",
                algorithm_name=self.algorithm_name,
                repetition_number=repetition_number,
                total_repetitions=total_repetitions,
                metadata={
                    "start_time": self.start_time,
                    "monitor_id": id(self),
                },
            )
            self.broker.emit(event)

            logger.info(
                f"Started repetition {repetition_number}/{total_repetitions} for {self.algorithm_name}"
            )

        except Exception as e:
            logger.error(f"Error starting repetition for {self.algorithm_name}: {e}")

    def finish_repetition(
        self,
        success: bool = True,
        result: Optional[Dict[str, Any]] = None,
        error_message: Optional[str] = None,
        repetition_number: int = 1,
    ) -> None:
        """
        Signal completion of repetition/trial/sample.

        Args:
            success: Whether repetition completed successfully
            result: Result data from the repetition
            error_message: Error message if failed
            repetition_number: Current repetition number
        """
        try:
            execution_time = time.time() - self.start_time

            event = RepetitionFinishedEvent(
                session_id=self.session_id,
                execution_level=ExecutionLevel.REPETITION,
                level_id=self.repetition_id,
                item_name=self.item_id or "unknown",
                algorithm_name=self.algorithm_name,
                repetition_number=repetition_number,
                success=success,
                result=result,
                error_message=error_message,
                execution_time=execution_time,
                metadata={
                    "history_entries": len(self.history_entries),
                    "final_progress": self.last_progress,
                },
            )
            self.broker.emit(event)

            status = "successfully" if success else f"with error: {error_message}"
            logger.info(
                f"Finished repetition {repetition_number} for {self.algorithm_name} {status}"
            )

        except Exception as e:
            logger.error(f"Error finishing repetition for {self.algorithm_name}: {e}")

    def get_summary(self) -> Dict[str, Any]:
        """
        Get monitoring summary for this algorithm execution.

        Returns:
            Summary dictionary with execution details
        """
        return {
            "algorithm_name": self.algorithm_name,
            "item_id": self.item_id,
            "repetition_id": self.repetition_id,
            "elapsed_time": time.time() - self.start_time,
            "last_progress": self.last_progress,
            "history_entries": len(self.history_entries),
            "session_id": self.session_id,
        }


class AlgorithmMonitorFactory:
    """Factory for creating algorithm monitors with consistent configuration."""

    def __init__(self, broker: ProgressBroker, session_id: str):
        """
        Initialize monitor factory.

        Args:
            broker: Progress broker for event emission
            session_id: Current session ID
        """
        self.broker = broker
        self.session_id = session_id
        self.monitors = {}  # Track active monitors

    def create_monitor(
        self,
        algorithm_name: str,
        item_id: Optional[str] = None,
        repetition_id: Optional[str] = None,
    ) -> AlgorithmMonitor:
        """
        Create new algorithm monitor.

        Args:
            algorithm_name: Name of the algorithm
            item_id: Optional item ID
            repetition_id: Optional repetition ID

        Returns:
            Configured algorithm monitor
        """
        monitor = AlgorithmMonitor(
            broker=self.broker,
            session_id=self.session_id,
            algorithm_name=algorithm_name,
            item_id=item_id,
            repetition_id=repetition_id,
        )

        # Track monitor
        monitor_key = f"{algorithm_name}_{item_id}_{repetition_id}"
        self.monitors[monitor_key] = monitor

        return monitor

    def get_monitor(
        self, algorithm_name: str, item_id: str, repetition_id: str
    ) -> Optional[AlgorithmMonitor]:
        """Get existing monitor by key."""
        monitor_key = f"{algorithm_name}_{item_id}_{repetition_id}"
        return self.monitors.get(monitor_key)

    def cleanup_monitor(
        self, algorithm_name: str, item_id: str, repetition_id: str
    ) -> None:
        """Remove monitor from tracking."""
        monitor_key = f"{algorithm_name}_{item_id}_{repetition_id}"
        self.monitors.pop(monitor_key, None)

    def get_all_monitors(self) -> Dict[str, AlgorithmMonitor]:
        """Get all active monitors."""
        return self.monitors.copy()
