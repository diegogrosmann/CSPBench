"""
Work Status Management.

Defines standardized status enumeration and management for work execution
across the entire CSPBench system, ensuring consistent status handling
and proper state transitions.

This module provides:
- Standardized status enumeration
- Status validation and normalization
- Status transition rules
- Helper functions for status management
"""

from enum import Enum


class BaseStatus(str, Enum):
    """
    Standardized work execution statuses across the entire system.
    
    This enumeration defines all possible states a work item can be in
    during its lifecycle, from submission to completion.
    
    Status Categories:
    - Incomplete: QUEUED, RUNNING, PAUSED, CANCELED
    - Final: COMPLETED, FAILED, ERROR
    """

    QUEUED = "queued"      # Work in queue, waiting to start
    RUNNING = "running"    # Work currently executing
    PAUSED = "paused"      # Work paused, can be resumed
    CANCELED = "canceled"  # Work canceled before completion
    COMPLETED = "completed"  # Work finished successfully
    FAILED = "failed"      # Failure prevented execution from proceeding
    ERROR = "error"        # Error during execution, but proceeded to end

    @property
    def is_final(self) -> bool:
        """
        Indicate if status represents a finalized (terminal) state.
        
        Returns:
            bool: True if status is terminal
        """
        return self in FINAL_STATUSES

    @property
    def is_incomplete(self) -> bool:
        """
        Indicate if status represents an incomplete (non-finalized) state.
        
        Returns:
            bool: True if status is incomplete
        """
        return self in INCOMPLETE_STATUSES


# Finalized (terminal) statuses - work has reached an end
# These statuses indicate that work was completed in some way
FINAL_STATUSES: set[BaseStatus] = {
    BaseStatus.COMPLETED,  # Success - work executed completely
    BaseStatus.FAILED,     # Failure - execution prevented from proceeding
    BaseStatus.ERROR,      # Error - execution proceeded but with problems
}

# Incomplete (non-finalized) statuses - work has not reached an end
# These statuses indicate that work has not yet been completed
INCOMPLETE_STATUSES: set[BaseStatus] = {
    BaseStatus.QUEUED,     # Waiting - work in queue for execution
    BaseStatus.RUNNING,    # Active - work being executed at the moment
    BaseStatus.PAUSED,     # Suspended - work paused, can be resumed
    BaseStatus.CANCELED,   # Canceled - work interrupted before end
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


def normalize_status(value) -> str:
    """
    Normalize a status value to lowercase string.

    Accepts:
    - BaseStatus instances (returns .value)
    - Strings (converts to lowercase and strip)
    - Other types (converts to string and applies lowercase)

    Args:
        value: Value to be normalized (BaseStatus, str, or other)

    Returns:
        str: Normalized string in lowercase

    Raises:
        ValueError: If normalized value doesn't correspond to a valid status
    """
    if isinstance(value, BaseStatus):
        return value.value

    if isinstance(value, str):
        normalized = value.strip().lower()
    else:
        normalized = str(value).strip().lower()

    # Validate if it's a known status
    valid_values = {s.value for s in BaseStatus}
    if normalized not in valid_values:
        raise ValueError(f"Invalid status: {value}. Valid: {sorted(valid_values)}")

    return normalized


def normalize_status(value: str | BaseStatus) -> str:
    """
    Normalize a status to validated lowercase string.

    Args:
        value: Status as BaseStatus or string (case-insensitive)
        
    Returns:
        str: Normalized string (one of BaseStatus values)
        
    Raises:
        ValueError: If value doesn't represent a valid status
    """
    if isinstance(value, BaseStatus):
        return value.value
    if not isinstance(value, str):  # fallback
        value = str(value)
    norm = value.strip().lower()
    if norm not in {s.value for s in BaseStatus}:
        raise ValueError(f"Invalid status: {value}")
    return norm
    norm = value.strip().lower()
    if norm not in {s.value for s in BaseStatus}:
        raise ValueError(f"Invalid status: {value}")
    return norm
