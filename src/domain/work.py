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
    WorkItem - Domain entity representing a work execution unit.

    Represents a single work item in the CSPBench system, containing
    configuration, status, timing information, and execution metadata.

    This entity handles status transitions, configuration management,
    and serialization to/from dictionary format for persistence.

    Attributes:
        id: Unique work identifier
        config: CSPBench configuration for execution
        status: Current execution status
        created_at: Creation timestamp (Unix time)
        updated_at: Last update timestamp (Unix time)
        output_path: Path for work output files
        error: Error message if execution failed
        extra: Additional metadata dictionary
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
        instance, converting from dict if necessary.

        Raises:
            TypeError: If config is not CSPBenchConfig or dict
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
        Create WorkItem from dictionary, ensuring proper config type.

        Factory method that creates a WorkItem instance from a dictionary
        representation, typically used when loading from persistence storage.

        Args:
            data: Dictionary containing work item data

        Returns:
            WorkItem: Created work item instance

        Raises:
            TypeError: If config data is not in expected format
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
        """Update the last modified timestamp."""
        self.updated_at = time.time()

    def can_transition(self, new: BaseStatus) -> bool:
        """
        Check if transition to new status is allowed.

        Args:
            new: Target status for transition

        Returns:
            bool: True if transition is allowed
        """
        return new in ALLOWEDSTATUS.get(self.status, set())

    def _set_status(self, new: BaseStatus) -> bool:
        """
        Attempt to set new status with validation.

        Args:
            new: Target status

        Returns:
            bool: True if status change was successful
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

        Returns:
            bool: True if transition successful
        """
        return self._set_status(BaseStatus.RUNNING)

    def mark_finished(self) -> bool:
        """
        Mark work item as finished/completed.

        Returns:
            bool: True if transition successful
        """
        ok = self._set_status(BaseStatus.COMPLETED)
        return ok

    def mark_error(self, error: str) -> bool:
        """
        Mark work item as failed with error message.

        Args:
            error: Error message to store

        Returns:
            bool: True if transition successful
        """
        ok = self._set_status(BaseStatus.FAILED)
        if ok:
            self.error = error
        return ok

    def pause(self) -> bool:
        """
        Pause work item execution.

        Returns:
            bool: True if transition successful
        """
        return self._set_status(BaseStatus.PAUSED)

    def resume(self) -> bool:
        """
        Resume paused work item.

        Returns:
            bool: True if transition successful
        """
        return self._set_status(BaseStatus.RUNNING)

    def cancel(self) -> bool:
        """
        Cancel work item execution.

        Returns:
            bool: True if transition successful
        """
        return self._set_status(BaseStatus.CANCELED)

    def restart(self) -> bool:
        """
        Restart work item (reset to QUEUED status).

        Only allows restart from terminal states (COMPLETED, FAILED).

        Returns:
            bool: True if restart successful
        """
        if self.status in (
            BaseStatus.COMPLETED,
            BaseStatus.FAILED,
        ):
            return False
        self.status = BaseStatus.QUEUED
        self._touch()
        return True

    # --- Representation ---
    def to_dict(self) -> dict[str, Any]:
        """
        Convert WorkItem to dictionary with raw timestamp values.

        Creates a dictionary representation suitable for persistence,
        with proper JSON serialization of complex fields.

        Returns:
            dict: Dictionary representation of work item
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
