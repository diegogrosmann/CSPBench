"""
Monitoring Interfaces for CSP Algorithms.

Pure domain module containing the monitoring protocol for algorithm execution.
This module defines the contract for monitoring algorithm progress, warnings,
and errors without any infrastructure dependencies.

The monitoring interface follows the domain-driven design principle, providing
a clean contract that can be implemented by various monitoring strategies
including null monitors, persistence monitors, composite monitors, and others.

Key Design Principles:
    - Protocol-based design for flexibility
    - Non-intrusive monitoring (best-effort, no exceptions)
    - Cancellation support for long-running operations
    - Extensible interface for future enhancements
    - Minimal dependencies (pure domain code)

Future extensions can add new methods to the protocol. To maintain backward
compatibility, provide defaults in implementations or use checks with hasattr.
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class AlgorithmMonitor(Protocol):
    """
    Minimal contract for algorithm execution monitoring.

    Defines the interface for monitoring algorithm execution including progress
    reporting, warning/error logging, and cancellation checking. All methods
    should be implemented as best-effort operations that do not raise exceptions
    to avoid interrupting algorithm flow.

    The protocol is designed to be lightweight and non-intrusive, allowing
    algorithms to report events without worrying about monitoring failures.
    Implementations should handle all errors gracefully and continue operation.

    Common Parameters:
        progress (float): Progress percentage from 0.0 to 100.0. Values may
            exceed 100.0 in special cases where algorithms report multi-phase
            progress or provide estimates.
        message (str): Short, descriptive event message for human consumption.
        **data (Any): Additional serializable payload data for structured
            logging and analysis. Should contain only simple types that can
            be easily serialized (str, int, float, bool, list, dict).

    Implementation Notes:
        - All methods should be non-blocking and fast
        - Implementations must not raise exceptions
        - Progress values should be monotonic when possible
        - Messages should be concise and informative
        - Data payloads should be kept minimal for performance
    """

    # Core monitoring methods --------------------------------------------------
    def on_progress(self, progress: float, message: str, /, **data: Any) -> None:
        """
        Report algorithm progress event.

        Called by algorithms to report execution progress. Implementations
        may apply throttling (time-based or minimum progress delta) to avoid
        flooding with progress updates.

        Args:
            progress (float): Progress percentage (0.0 to 100.0, may exceed).
            message (str): Short progress description.
            **data (Any): Additional progress metadata.

        Note:
            This is a high-frequency method that may be called many times
            during algorithm execution. Implementations should be efficient
            and consider throttling strategies.
        """
        ...  # pragma: no cover

    def on_warning(self, message: str, /, **data: Any) -> None:
        """
        Report non-fatal warning event.

        Called by algorithms to report warnings that don't prevent execution
        from continuing but should be logged for analysis and debugging.

        Args:
            message (str): Warning message description.
            **data (Any): Additional warning context data.

        Examples:
            - Parameter values outside recommended ranges
            - Performance degradation warnings
            - Data quality issues that don't prevent processing
        """
        ...  # pragma: no cover

    # Optional monitoring methods ----------------------------------------------
    def on_error(
        self, message: str, exc: Exception | None = None, /, **data: Any
    ) -> None:
        """
        Report error event (optional method).

        Called by algorithms to report errors that may or may not be fatal.
        The presence of an exception object indicates a more serious error
        than a simple message-only report.

        Args:
            message (str): Error message description.
            exc (Optional[Exception]): Exception object if available.
            **data (Any): Additional error context data.

        Note:
            This method is optional. Implementations may treat errors as
            warnings or implement specific error handling strategies.
        """
        ...  # pragma: no cover

    def is_cancelled(self) -> bool:
        """
        Check if execution has been cancelled externally.

        Allows algorithms to query cancellation status during long-running
        operations and abort gracefully when requested. This enables
        responsive cancellation without forceful process termination.

        Returns:
            bool: True if execution should be cancelled, False otherwise.

        Default Behavior:
            Returns False (not cancelled) by default.

        Usage Pattern:
            Algorithms should check this method periodically during long
            loops or expensive operations and break gracefully when True.
        """
        return False  # pragma: no cover
