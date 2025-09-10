"""
Monitoring Interfaces for CSP Algorithms.

Pure domain module containing ONLY the monitoring protocol (`AlgorithmMonitor`).
Concrete implementations (Null, Composite, Throttled, persistence adapters, etc.)
should live outside this unit to keep the domain minimal and without unnecessary
dependencies.

Future extensions can add new methods to the protocol; to maintain backward
compatibility, provide defaults here or use checks with `hasattr`.
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class AlgorithmMonitor(Protocol):
    """
    Minimal contract for algorithm execution monitoring.

    All methods should be best-effort: implementations should **not** raise
    exceptions to avoid interrupting the algorithm flow.

    Common Parameters:
        progress (float): 0.0 to 100.0 representing percentage (or can extrapolate >100 if makes sense in special cases).
        message (str): Short event description.
        **data: Additional serializable/simple payload.
    """

    # Main methods --------------------------------------------------
    def on_progress(self, progress: float, message: str, /, **data: Any) -> None:
        """
        Progress event.

        Recommendation: Implementations can apply throttling (time or minimum
        progress delta) to avoid flooding.
        """
        ...  # pragma: no cover

    def on_warning(self, message: str, /, **data: Any) -> None:
        """Non-fatal warning event."""
        ...  # pragma: no cover

    # Optional methods ---------------------------------------------------
    def on_error(
        self, message: str, exc: Exception | None = None, /, **data: Any
    ) -> None:
        """Error event (optional)."""
        ...  # pragma: no cover

    def is_cancelled(self) -> bool:
        """
        Indicate if execution was canceled externally.

        Algorithms can query in long loops and abort gracefully.
        Default: False.
        """
        return False  # pragma: no cover
