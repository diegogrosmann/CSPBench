"""
Work Status Management.

Defines standardized status enumeration and management for work execution
across the entire CSPBench system, ensuring consistent status handling
and proper state transitions.

This module provides the core status management infrastructure including:
- Standardized status enumeration with clear semantics
- Status validation and normalization utilities
- Status transition rules and validation
- Helper functions for status categorization
- Type-safe status operations

The status system supports the complete work lifecycle from submission
through completion, with proper categorization into incomplete and
final states to support different operational needs.

Status Flow:
    QUEUED -> RUNNING -> {COMPLETED, FAILED, ERROR}
           -> PAUSED -> RUNNING
           -> CANCELED

Categories:
    - Incomplete: Work has not reached a terminal state
    - Final: Work has completed (successfully or with errors)
"""

from enum import Enum


class BaseStatus(str, Enum):
    """
    Standardized work execution statuses across the entire system.

    This enumeration defines all possible states a work item can be in
    during its lifecycle, from submission to completion. Each status
    has clear semantics and defined transition rules.

    Status Categories:
        Incomplete States: QUEUED, RUNNING, PAUSED, CANCELED
            - Work has not reached a final conclusion
            - May be resumed or restarted

        Final States: COMPLETED, FAILED, ERROR
            - Work has reached a terminal state
            - Represents different completion outcomes

    Status Descriptions:
        QUEUED: Work submitted and waiting for execution
        RUNNING: Work currently being executed
        PAUSED: Work temporarily suspended, can be resumed
        CANCELED: Work canceled before completion
        COMPLETED: Work finished successfully
        FAILED: Work failed and could not complete
        ERROR: Work completed but with errors
    """

    QUEUED = "queued"  # Work in queue, waiting to start
    RUNNING = "running"  # Work currently executing
    PAUSED = "paused"  # Work paused, can be resumed
    CANCELED = "canceled"  # Work canceled before completion
    COMPLETED = "completed"  # Work finished successfully
    FAILED = "failed"  # Failure prevented execution from proceeding
    ERROR = "error"  # Error during execution, but proceeded to end

    @property
    def is_final(self) -> bool:
        """
        Check if status represents a finalized (terminal) state.

        Terminal states indicate that work has reached a conclusion
        and will not proceed further without explicit restart.

        Returns:
            bool: True if status is terminal (COMPLETED, FAILED, ERROR).
        """
        return self in FINAL_STATUSES

    @property
    def is_incomplete(self) -> bool:
        """
        Check if status represents an incomplete (non-finalized) state.

        Incomplete states indicate that work has not yet reached a
        final conclusion and may continue or be resumed.

        Returns:
            bool: True if status is incomplete (QUEUED, RUNNING, PAUSED, CANCELED).
        """
        return self in INCOMPLETE_STATUSES


# Finalized (terminal) statuses - work has reached an end
# These statuses indicate that work was completed in some way
FINAL_STATUSES: set[BaseStatus] = {
    BaseStatus.COMPLETED,  # Success - work executed completely
    BaseStatus.FAILED,  # Failure - execution prevented from proceeding
    BaseStatus.ERROR,  # Error - execution proceeded but with problems
}

# Incomplete (non-finalized) statuses - work has not reached an end
# These statuses indicate that work has not yet been completed
INCOMPLETE_STATUSES: set[BaseStatus] = {
    BaseStatus.QUEUED,  # Waiting - work in queue for execution
    BaseStatus.RUNNING,  # Active - work being executed at the moment
    BaseStatus.PAUSED,  # Suspended - work paused, can be resumed
    BaseStatus.CANCELED,  # Canceled - work interrupted before end
}

# Status transition rules - defines allowed transitions between statuses
ALLOWEDSTATUS: dict[BaseStatus, set[BaseStatus]] = {
    BaseStatus.QUEUED: {BaseStatus.RUNNING, BaseStatus.CANCELED, BaseStatus.PAUSED},
    BaseStatus.RUNNING: {
        BaseStatus.PAUSED,
        BaseStatus.COMPLETED,
        BaseStatus.FAILED,
        BaseStatus.CANCELED,
        BaseStatus.ERROR,
    },
    BaseStatus.PAUSED: {BaseStatus.RUNNING, BaseStatus.CANCELED, BaseStatus.QUEUED},
    BaseStatus.CANCELED: {BaseStatus.QUEUED, BaseStatus.PAUSED},
}


def normalize_status(value: str | BaseStatus) -> str:
    """
    Normalize a status to validated lowercase string.

    Accepts BaseStatus instances or string values and normalizes them
    to consistent lowercase string format with validation against
    known status values.

    Args:
        value (Union[str, BaseStatus]): Status as BaseStatus enum or
            string (case-insensitive).

    Returns:
        str: Normalized string (one of BaseStatus values).

    Raises:
        ValueError: If value doesn't represent a valid status.

    Examples:
        >>> normalize_status(BaseStatus.RUNNING)
        'running'
        >>> normalize_status('COMPLETED')
        'completed'
        >>> normalize_status('  Queued  ')
        'queued'
    """
    if isinstance(value, BaseStatus):
        return value.value
    if not isinstance(value, str):  # fallback
        value = str(value)
    norm = value.strip().lower()
    if norm not in {s.value for s in BaseStatus}:
        raise ValueError(f"Invalid status: {value}")
    return norm
