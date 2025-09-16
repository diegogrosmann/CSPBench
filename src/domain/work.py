"""
Domain Work Entity.

This module contains the WorkItem domain entity representing a work execution
unit in the CSPBench system. It handles the complete lifecycle of work items
including status management, configuration storage, and persistence operations.

The WorkItem entity encapsulates all information needed to track and manage
work execution from submission through completion, including configuration,
status transitions, timing information, and execution metadata.

Features:
    - Status transition management with validation
    - Configuration serialization and deserialization
    - Timestamp tracking for lifecycle events
    - Error handling and metadata storage
    - Dictionary conversion for persistence operations
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from typing import Any, Dict, Optional

from src.domain.status import ALLOWEDSTATUS, BaseStatus

from .config import CSPBenchConfig

# Alias to maintain compatibility
WorkStatus = BaseStatus


@dataclass
class WorkItem:
    """
    Domain entity representing a work execution unit.

    Represents a single work item in the CSPBench system, containing
    configuration, status, timing information, and execution metadata.
    This entity handles status transitions, configuration management,
    and serialization to/from dictionary format for persistence.

    The WorkItem follows domain-driven design principles, encapsulating
    business logic for work lifecycle management while remaining independent
    of infrastructure concerns.

    Attributes:
        id (str): Unique work identifier.
        config (CSPBenchConfig): CSPBench configuration for execution.
        status (BaseStatus): Current execution status.
        created_at (float): Creation timestamp (Unix time).
        updated_at (float): Last update timestamp (Unix time).
        output_path (str): Path for work output files.
        error (Optional[str]): Error message if execution failed.
        extra (Optional[Dict[str, Any]]): Additional metadata dictionary.
    """

    id: str
    config: CSPBenchConfig
    status: BaseStatus
    created_at: float
    updated_at: float
    output_path: str
    error: Optional[str] = None
    extra: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        """
        Validate and convert config to CSPBenchConfig if needed.

        Ensures that the config attribute is always a proper CSPBenchConfig
        instance, converting from dict if necessary. This provides flexibility
        in work item creation while maintaining type safety.

        Raises:
            TypeError: If config is not CSPBenchConfig or dict.
        """
        if isinstance(self.config, dict):
            self.config = CSPBenchConfig.from_dict(self.config)
        elif not isinstance(self.config, CSPBenchConfig):
            raise TypeError(
                f"config must be CSPBenchConfig or dict, got {type(self.config)}"
            )

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> WorkItem:
        """
        Create WorkItem from dictionary representation.

        Factory method that creates a WorkItem instance from a dictionary
        representation, typically used when loading from persistence storage.
        Handles JSON deserialization and type conversion for all fields.

        Args:
            data (Dict[str, Any]): Dictionary containing work item data with keys:
                - id: Work identifier
                - config_json: Configuration as JSON string or dict
                - status: Status string value
                - created_at: Creation timestamp
                - updated_at: Update timestamp
                - output_path: Output path string
                - error: Optional error message
                - extra_json: Extra metadata as JSON string or dict

        Returns:
            WorkItem: Created work item instance.

        Raises:
            TypeError: If config data is not in expected format.
            ValueError: If required fields are missing or invalid.
        """
        config_data = data.get("config_json")

        # Handle config JSON if it's a string
        if isinstance(config_data, str):
            config_data = json.loads(config_data)

        if isinstance(config_data, dict):
            config = CSPBenchConfig.from_dict(config_data)
        elif isinstance(config_data, CSPBenchConfig):
            config = config_data
        else:
            raise TypeError(
                f"config must be CSPBenchConfig or dict, got {type(config_data)}"
            )

        # Handle extra JSON if it's a string
        extra = data.get("extra_json")
        if isinstance(extra, str):
            extra = json.loads(extra)

        # Timestamps should be floats now
        created_at = float(data["created_at"])
        updated_at = float(data["updated_at"])

        return cls(
            id=data["id"],
            config=config,
            status=BaseStatus(data["status"]),
            created_at=created_at,
            updated_at=updated_at,
            output_path=data["output_path"],
            error=data.get("error"),
            extra=extra,
        )

    # --- Internal helpers ---
    def _touch(self) -> None:
        """
        Update the last modified timestamp.

        Updates the updated_at field with the current timestamp to track
        when the work item was last modified.
        """
        self.updated_at = time.time()

    def can_transition(self, new: BaseStatus) -> bool:
        """
        Check if transition to new status is allowed.

        Validates whether a status transition is permitted based on the
        current status and the defined transition rules in ALLOWEDSTATUS.

        Args:
            new (BaseStatus): Target status for transition.

        Returns:
            bool: True if transition is allowed, False otherwise.
        """
        return new in ALLOWEDSTATUS.get(self.status, set())

    def _set_status(self, new: BaseStatus) -> bool:
        """
        Attempt to set new status with validation.

        Internal method that attempts to change the work item status
        after validating that the transition is allowed. Updates the
        timestamp if the transition is successful.

        Args:
            new (BaseStatus): Target status.

        Returns:
            bool: True if status change was successful, False otherwise.
        """
        if not self.can_transition(new):
            return False
        self.status = new
        self._touch()
        return True

    # --- Status transitions ---
    def mark_running(self) -> bool:
        """
        Mark work item as running.

        Transitions the work item to RUNNING status if the current status
        allows this transition.

        Returns:
            bool: True if transition successful, False otherwise.
        """
        return self._set_status(BaseStatus.RUNNING)

    def mark_finished(self) -> bool:
        """
        Mark work item as finished/completed.

        Transitions the work item to COMPLETED status if the current status
        allows this transition.

        Returns:
            bool: True if transition successful, False otherwise.
        """
        ok = self._set_status(BaseStatus.COMPLETED)
        return ok

    def mark_error(self, error: str) -> bool:
        """
        Mark work item as failed with error message.

        Transitions the work item to FAILED status and stores the provided
        error message if the current status allows this transition.

        Args:
            error (str): Error message to store.

        Returns:
            bool: True if transition successful, False otherwise.
        """
        ok = self._set_status(BaseStatus.FAILED)
        if ok:
            self.error = error
        return ok

    def pause(self) -> bool:
        """
        Pause work item execution.

        Transitions the work item to PAUSED status if the current status
        allows this transition.

        Returns:
            bool: True if transition successful, False otherwise.
        """
        return self._set_status(BaseStatus.PAUSED)

    def resume(self) -> bool:
        """
        Resume paused work item.

        Transitions the work item from PAUSED to RUNNING status if the
        current status allows this transition.

        Returns:
            bool: True if transition successful, False otherwise.
        """
        return self._set_status(BaseStatus.RUNNING)

    def cancel(self) -> bool:
        """
        Cancel work item execution.

        Transitions the work item to CANCELED status if the current status
        allows this transition.

        Returns:
            bool: True if transition successful, False otherwise.
        """
        return self._set_status(BaseStatus.CANCELED)

    def restart(self) -> bool:
        """
        Restart work item (reset to QUEUED status).

        Resets the work item to QUEUED status, allowing it to be executed
        again. Only allows restart from terminal states (COMPLETED, FAILED)
        and non-terminal states (PAUSED, CANCELED).

        Returns:
            bool: True if restart successful, False otherwise.
        """
        if self.status in (
            BaseStatus.COMPLETED,
            BaseStatus.FAILED,
            BaseStatus.PAUSED,
            BaseStatus.CANCELED,
        ):
            self.status = BaseStatus.QUEUED
            self._touch()
            return True
        return False

    # --- Representation ---
    def to_dict(self) -> dict[str, Any]:
        """
        Convert WorkItem to dictionary with serialized fields.

        Creates a dictionary representation suitable for persistence,
        with proper JSON serialization of complex fields like configuration
        and extra metadata.

        Returns:
            Dict[str, Any]: Dictionary representation of work item with keys:
                - id: Work identifier
                - config_json: Configuration as JSON string
                - status: Status value string
                - created_at: Creation timestamp as float
                - updated_at: Update timestamp as float
                - output_path: Output path string
                - error: Error message or None
                - extra_json: Extra metadata as JSON string
        """
        return {
            "id": self.id,
            "config_json": json.dumps(self.config.to_dict()) if self.config else None,
            "status": self.status.value,
            "created_at": self.created_at,  # Keep as float
            "updated_at": self.updated_at,  # Keep as float
            "output_path": self.output_path,
            "error": self.error,
            "extra_json": json.dumps(self.extra) if self.extra else json.dumps({}),
        }
