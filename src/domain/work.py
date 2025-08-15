from __future__ import annotations

from dataclasses import dataclass, asdict
from enum import Enum
import time
from typing import Any, Dict, Optional

from .config import CSPBenchConfig


class WorkStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    PAUSED = "paused"
    CANCELED = "canceled"
    FINISHED = "finished"
    ERROR = "error"


_ALLOWED: dict[WorkStatus, set[WorkStatus]] = {
    WorkStatus.QUEUED: {WorkStatus.RUNNING, WorkStatus.CANCELED},
    WorkStatus.RUNNING: {
        WorkStatus.PAUSED,
        WorkStatus.FINISHED,
        WorkStatus.ERROR,
        WorkStatus.CANCELED,
    },
    WorkStatus.PAUSED: {WorkStatus.RUNNING, WorkStatus.CANCELED},
    WorkStatus.FINISHED: {WorkStatus.QUEUED},
    WorkStatus.ERROR: {WorkStatus.QUEUED},
    WorkStatus.CANCELED: {WorkStatus.QUEUED},
}


@dataclass
class WorkItem:
    id: str
    config: CSPBenchConfig
    status: WorkStatus
    created_at: float
    updated_at: float
    output_path: str
    error: Optional[str] = None
    extra: Optional[Dict[str, Any]] = None

    # --- internal helpers ---
    def _touch(self) -> None:
        self.updated_at = time.time()

    def can_transition(self, new: WorkStatus) -> bool:
        return new in _ALLOWED.get(self.status, set())

    def _set_status(self, new: WorkStatus) -> bool:
        if not self.can_transition(new):
            return False
        self.status = new
        self._touch()
        return True

    # --- transitions ---
    def mark_running(self) -> bool:
        return self._set_status(WorkStatus.RUNNING)

    def mark_finished(self) -> bool:
        ok = self._set_status(WorkStatus.FINISHED)
        return ok

    def mark_error(self, error: str) -> bool:
        ok = self._set_status(WorkStatus.ERROR)
        if ok:
            self.error = error
        return ok

    def pause(self) -> bool:
        return self._set_status(WorkStatus.PAUSED)

    def resume(self) -> bool:
        return self._set_status(WorkStatus.RUNNING)

    def cancel(self) -> bool:
        return self._set_status(WorkStatus.CANCELED)

    def restart(self) -> bool:
        if self.status not in (
            WorkStatus.FINISHED,
            WorkStatus.ERROR,
            WorkStatus.CANCELED,
        ):
            return False
        self.status = WorkStatus.QUEUED
        self._touch()
        return True

    # --- representation ---
    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["status"] = self.status.value
        return d
