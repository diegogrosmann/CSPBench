from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import json
import time
from typing import Any, Dict, Optional

from src.domain.status import ALLOWEDSTATUS, BaseStatus

from .config import CSPBenchConfig

# Alias para manter compatibilidade
WorkStatus = BaseStatus


@dataclass
class WorkItem:
    id: str
    config: CSPBenchConfig
    status: BaseStatus
    created_at: float
    updated_at: float
    output_path: str
    error: Optional[str] = None
    extra: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        """Validate and convert config to CSPBenchConfig if needed."""
        if isinstance(self.config, dict):
            self.config = CSPBenchConfig.from_dict(self.config)
        elif not isinstance(self.config, CSPBenchConfig):
            raise TypeError(
                f"config must be CSPBenchConfig or dict, got {type(self.config)}"
            )

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkItem":
        """Create WorkItem from dictionary, ensuring proper config type."""
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

    # --- internal helpers ---
    def _touch(self) -> None:
        self.updated_at = time.time()

    def can_transition(self, new: BaseStatus) -> bool:
        return new in ALLOWEDSTATUS.get(self.status, set())

    def _set_status(self, new: BaseStatus) -> bool:
        if not self.can_transition(new):
            return False
        self.status = new
        self._touch()
        return True

    # --- transitions ---
    def mark_running(self) -> bool:
        return self._set_status(BaseStatus.RUNNING)

    def mark_finished(self) -> bool:
        ok = self._set_status(BaseStatus.COMPLETED)
        return ok

    def mark_error(self, error: str) -> bool:
        ok = self._set_status(BaseStatus.FAILED)
        if ok:
            self.error = error
        return ok

    def pause(self) -> bool:
        return self._set_status(BaseStatus.PAUSED)

    def resume(self) -> bool:
        return self._set_status(BaseStatus.RUNNING)

    def cancel(self) -> bool:
        return self._set_status(BaseStatus.CANCELED)

    def restart(self) -> bool:
        if self.status not in (
            BaseStatus.COMPLETED,
            BaseStatus.FAILED,
            BaseStatus.CANCELED,
        ):
            return False
        self.status = BaseStatus.QUEUED
        self._touch()
        return True

    # --- representation ---
    def to_dict(self) -> dict[str, Any]:
        """Convert WorkItem to dictionary with raw timestamp values."""

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
