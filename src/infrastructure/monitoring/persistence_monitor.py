"""PersistenceMonitor

Monitor responsible for persisting execution events (progress, warnings, errors)
using `ExecutionScopedPersistence`.

Layer: infrastructure.
Depends only on the domain protocol `AlgorithmMonitor` and the execution-scoped
persistence wrapper.

Features:
- on_progress: records progress (with optional throttling)
- on_warning: records warning associated with the unit
- on_error: records error (or warning if only message)
- is_cancelled: allows external cancellation via ExecutionController or internal flag

Throttling: avoids excessive progress recording based on minimum percentage delta
and/or minimum time interval. Always records initial 0% (if invoked) and final 100%.
"""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any, Optional

from src.domain.monitoring import AlgorithmMonitor
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
    ExecutionScopedPersistence,
)

if TYPE_CHECKING:
    from src.infrastructure.execution_control import ExecutionController


class PersistenceMonitor(AlgorithmMonitor):
    """
    Monitor that persists events to the store via `ExecutionScopedPersistence`.

    Provides persistence for algorithm monitoring events with configurable
    throttling to prevent excessive database writes while ensuring important
    progress milestones are captured.

    Attributes:
        _exec_store: Execution-scoped persistence wrapper
        _min_delta_pct: Minimum progress delta to persist
        _min_interval_s: Minimum interval between persistences
        _last_progress: Last recorded progress value
        _last_time: Last persistence timestamp
        _cancelled: Internal cancellation flag
        _execution_controller: Optional execution controller for status checks
        _logger: Logger instance for this class

    Example:
        >>> exec_store = ExecutionScopedPersistence(work_id="123", unit_id="algo1")
        >>> monitor = PersistenceMonitor(exec_store, min_delta_pct=1.0, min_interval_s=0.5)
        >>> monitor.on_progress(50.0, "Processing data...")
        >>> monitor.on_warning("Low memory detected")
    """

    __slots__ = (
        "_exec_store",
        "_min_delta_pct",
        "_min_interval_s",
        "_last_progress",
        "_last_time",
        "_cancelled",
        "_execution_controller",
        "_logger",
    )

    def __init__(
        self,
        exec_store: ExecutionScopedPersistence,
        *,
        min_delta_pct: float = 0.01,
        min_interval_s: float = 0.00,
        execution_controller: Optional[ExecutionController] = None,
    ) -> None:
        """
        Initialize PersistenceMonitor with execution store and throttling settings.

        Args:
            exec_store: Execution-scoped persistence wrapper for event storage
            min_delta_pct: Minimum progress delta to persist (default 1.0%)
            min_interval_s: Minimum interval between persistences (default 0.5s)
            execution_controller: Optional controller for status checks

        Note:
            Progress throttling helps reduce database load while ensuring
            significant progress changes are captured.
        """
        self._exec_store = exec_store
        self._min_delta_pct = min_delta_pct
        self._min_interval_s = min_interval_s
        self._last_progress = -float("inf")
        self._last_time = 0.0
        self._cancelled = False
        self._execution_controller = execution_controller
        self._logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    @classmethod
    def with_execution_controller(
        cls,
        exec_store: ExecutionScopedPersistence,
        execution_controller: ExecutionController,
        *,
        min_delta_pct: float = 0.01,
        min_interval_s: float = 0.00,
    ) -> PersistenceMonitor:
        """
        Convenience factory method to create PersistenceMonitor with ExecutionController.

        Creates a PersistenceMonitor instance configured with an ExecutionController
        for enhanced cancellation detection and resource management integration.

        Args:
            exec_store: Execution-scoped persistence store
            execution_controller: Controller for status checks and resource management
            min_delta_pct: Minimum progress delta to persist
            min_interval_s: Minimum time interval between persistencies

        Returns:
            PersistenceMonitor: Configured PersistenceMonitor instance

        Example:
            >>> controller = ExecutionController("work_123")
            >>> exec_store = ExecutionScopedPersistence("work_123", "unit_1")
            >>> monitor = PersistenceMonitor.with_execution_controller(exec_store, controller)
        """
        return cls(
            exec_store=exec_store,
            min_delta_pct=min_delta_pct,
            min_interval_s=min_interval_s,
            execution_controller=execution_controller,
        )

    # ---------------------------------------------------------------------
    # AlgorithmMonitor API Implementation
    # ---------------------------------------------------------------------

    def on_progress(self, progress: float, message: str, /, **data: Any) -> None:
        """
        Persist progress event with optional throttling (non-critical operation).

        Records progress to the execution store if throttling conditions are met.
        Progress persistence errors are logged but don't interrupt execution.

        Args:
            progress: Progress percentage (0.0 to 100.0)
            message: Progress description message
            **data: Additional progress data (currently unused)

        Note:
            This is a non-critical operation - errors are logged but execution continues.
            Throttling prevents excessive database writes while capturing key milestones.
        """
        try:
            now = time.time()
            pct_ok = (progress - self._last_progress) >= self._min_delta_pct
            time_ok = (now - self._last_time) >= self._min_interval_s

            if pct_ok and time_ok:
                self._exec_store.add_progress(progress, message)
                self._last_progress = progress
                self._last_time = now
        except Exception as e:
            # Progress errors should not be critical - just log and continue
            self._logger.warning(
                "Non-critical error persisting progress (continuing execution): %s (progress=%.2f, message=%s)",
                e,
                progress,
                message,
                exc_info=False,  # Reduced verbosity for progress errors
            )

    def on_warning(self, message: str, /, **data: Any) -> None:
        """
        Persist warning event.

        Records a warning associated with the current execution unit.
        Warnings indicate non-fatal issues that should be tracked.

        Args:
            message: Warning message describing the issue
            **data: Additional context data for the warning

        Raises:
            Exception: Re-raises any persistence errors after logging
        """
        try:
            context = data if data else None
            # Use the execution scoped unit_id as the unit_id parameter
            self._exec_store.unit_warning(self._exec_store._unit_id, message, context)
        except Exception as e:
            self._logger.error(
                "Error persisting warning: %s (message=%s, data=%s)",
                e,
                message,
                data,
                exc_info=True,
            )

    def on_error(
        self, message: str, exc: Exception | None = None, /, **data: Any
    ) -> None:
        """
        Persist error event (or warning if only message).

        Records either a full error (if exception provided) or treats
        it as a categorized warning (if only message provided).

        Args:
            message: Error/warning message
            exc: Optional exception object for full error recording
            **data: Additional context data

        Note:
            If exc is None, the event is recorded as a warning with error categorization.
            If exc is provided, it's recorded as a full error event.
        """
        try:
            if exc is not None:
                # Record full error
                self._exec_store.unit_error(self._exec_store._unit_id, exc)
            else:
                # No concrete exception -> treat as categorized warning
                context = (
                    {"as_error_message": message, **data}
                    if data
                    else {"as_error_message": message}
                )
                self._exec_store.unit_warning(
                    self._exec_store._unit_id, message, context
                )
        except Exception as e:
            self._logger.error(
                "Error persisting error/warning: %s (message=%s, exc=%s, data=%s)",
                e,
                message,
                exc,
                data,
                exc_info=True,
            )

    def is_cancelled(self) -> bool:
        """
        Check if execution is cancelled via ExecutionController or internal flag.

        Checks both internal cancellation flag and ExecutionController status.
        Only CANCELED status from ExecutionController is considered cancellation.

        Returns:
            bool: True if execution should be cancelled, False otherwise

        Note:
            COMPLETED/FAILED work status doesn't trigger cancellation since
            individual executions may still be in progress when pipeline finishes.
        """
        # Check internal flag first (explicit cancellation by user/system)
        if self._cancelled:
            return True

        # Check ExecutionController if available
        if self._execution_controller:
            try:
                from src.domain.status import BaseStatus

                status = self._execution_controller.check_status()
                # Only CANCELED status means cancelled
                # Note: We don't cancel on COMPLETED/FAILED work status because
                # individual executions may still be in progress when the pipeline
                # finishes and updates the work status
                if status == BaseStatus.CANCELED:
                    return True
            except Exception as e:
                self._logger.warning(
                    "Error checking cancellation status via ExecutionController: %s",
                    e,
                    exc_info=True,
                )

        # Default to not cancelled
        return False

    # ------------------------------------------------------------------
    # Extension Methods
    # ------------------------------------------------------------------

    def cancel(self) -> None:
        """
        Signal external cancellation.

        Sets internal cancellation flag to True, which will be detected
        by subsequent calls to is_cancelled().
        """
        self._cancelled = True

    # Convenience for algorithms that want direct access (avoid extra import)
    @property
    def execution(self) -> ExecutionScopedPersistence:
        """
        Return the execution-scoped persistence store.

        Provides direct access to the underlying persistence store for
        algorithms that need additional persistence operations beyond
        the standard monitoring interface.

        Returns:
            ExecutionScopedPersistence: The execution store instance
        """
        return self._exec_store
