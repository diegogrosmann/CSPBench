from __future__ import annotations

import os
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Optional

from src.infrastructure.persistence.work_service_persistence import (
    WorkServicePersistence,
)
from src.domain.config import CSPBenchConfig
from src.domain.work import WorkItem, WorkStatus
from src.domain.status import BaseStatus
from src.infrastructure.utils.path_utils import get_output_base_directory


class WorkManager:
    def __init__(self, repository: WorkServicePersistence):
        self._repo = repository
        self._lock = threading.RLock()

    # util
    def _now(self) -> float:
        return time.time()

    def _new_id(self) -> str:
        return uuid.uuid4().hex[:12]

    # --- public API ---
    def submit(
        self, *, config: CSPBenchConfig, extra: dict[str, Any] | None = None
    ) -> str:
        with self._lock:
            wid = self._new_id()
            work_dir = get_output_base_directory() / wid
            item = WorkItem(
                id=wid,
                config=config,
                status=BaseStatus.QUEUED,
                created_at=self._now(),
                updated_at=self._now(),
                output_path=str(work_dir),
                extra=extra or {},
            )
            self._repo.add(item)

            return wid

    def get(self, work_id: str) -> dict[str, Any] | None:
        with self._lock:
            item = self._repo.get(work_id)
            return item.to_dict() if item else None

    def list(self) -> list[dict[str, Any]]:
        with self._lock:
            return [w.to_dict() for w in self._repo.list()]

    def get_status(self, work_id: str) -> Optional[str]:
        with self._lock:
            item = self._repo.get(work_id)
            return item.status.value if item else None

    def mark_running(self, work_id: str) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            ok = item.mark_running()
            if ok:
                self._repo.update(item)
            return ok

    def mark_finished(self, work_id: str) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            ok = item.mark_finished()
            if ok:
                self._repo.update(item)
            return ok

    def mark_error(self, work_id: str, error: str) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            ok = item.mark_error(error)
            if ok:
                self._repo.update(item)
            return ok

    def pause(self, work_id: str) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            ok = item.pause()
            if ok:
                self._repo.update(item)
            return ok

    def resume(self, work_id: str) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            ok = item.resume()
            if ok:
                self._repo.update(item)
            return ok

    def cancel(self, work_id: str) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            ok = item.cancel()
            if ok:
                self._repo.update(item)
            return ok

    def restart(self, work_id: str) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            ok = item.restart()
            if ok:
                self._repo.update(item)
            return ok

    def wait_until_terminal(
        self, work_id: str, timeout: float | None = None, poll_interval: float = 0.5
    ) -> Optional[str]:
        """
        Bloqueia até o work atingir um estado terminal (FINISHED, ERROR, CANCELED) ou até 'timeout' (segundos).
        Retorna:
          - status final (str) se terminal
          - None se o work não existir
          - último status observado (não-terminal) se estourar o timeout
        """
        start = self._now()
        last_status: Optional[str] = None
        terminal_statuses = [
            BaseStatus.COMPLETED.value,
            BaseStatus.FAILED.value,
            BaseStatus.CANCELED.value,
        ]

        while True:
            status = self.get_status(work_id)
            if status is None:
                return None
            last_status = status
            if status in terminal_statuses:
                return status

            if timeout is not None and (self._now() - start) >= timeout:
                return last_status

            time.sleep(poll_interval)
