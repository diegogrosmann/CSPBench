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
    """Monitor that persists events to the store via `ExecutionScopedPersistence`.

    Args:
        exec_store: Execution-scoped persistence wrapper.
        min_delta_pct: Minimum progress delta to persist (default 1.0).
        min_interval_s: Minimum interval between persistences (default 0.5s).
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

        Args:
            exec_store: Execution-scoped persistence store
            execution_controller: Controller for status checks and resource management
            min_delta_pct: Minimum progress delta to persist
            min_interval_s: Minimum time interval between persistencies

        Returns:
            Configured PersistenceMonitor instance
        """
        return cls(
            exec_store=exec_store,
            min_delta_pct=min_delta_pct,
            min_interval_s=min_interval_s,
            execution_controller=execution_controller,
        )

    # ---------------------------------------------------------------------
    # API AlgorithmMonitor
    # ---------------------------------------------------------------------
    def on_progress(
        self, progress: float, message: str, /, **data: Any
    ) -> None:
        """Persist progress event with optional throttling (non-critical operation)."""
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
                exc_info=False  # Reduced verbosity for progress errors
            )

    def on_warning(self, message: str, /, **data: Any) -> None:
        """Persist warning event."""
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
        """Persist error event (or warning if only message)."""
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
                self._exec_store.unit_warning(self._exec_store._unit_id, message, context)
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
        """Check if execution is cancelled via ExecutionController or internal flag."""
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
    # Extension
    # ------------------------------------------------------------------
    def cancel(self) -> None:
        """Signal external cancellation."""
        self._cancelled = True

    # Convenience for algorithms that want direct access (avoid extra import)
    @property
    def execution(self) -> ExecutionScopedPersistence:
        """Return the execution-scoped persistence store."""
        return self._exec_store

    def cancel(self) -> None:
        """Signal external cancellation."""
        self._cancelled = True
