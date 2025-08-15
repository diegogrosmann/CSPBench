from __future__ import annotations

import os
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Optional


from src.infrastructure.monitoring.monitor_interface import Monitor
from src.domain.config import CSPBenchConfig
from src.domain.work import WorkItem, WorkStatus
from .repository import InMemoryWorkRepository, WorkRepository

# Singleton simples
_singleton_lock = threading.Lock()
_singleton_instance: "WorkManager | None" = None


def get_work_manager() -> "WorkManager":
    global _singleton_instance
    if _singleton_instance is None:
        with _singleton_lock:
            if _singleton_instance is None:
                _singleton_instance = WorkManager()
    return _singleton_instance


class WorkManager:
    def __init__(self, repository: WorkRepository | None = None):
        self._repo = repository or InMemoryWorkRepository()
        self._lock = threading.RLock()

    # util
    def _now(self) -> float:
        return time.time()

    def _new_id(self) -> str:
        return uuid.uuid4().hex[:12]

    # --- public API ---
    def submit(
        self, *, config: CSPBenchConfig, extra: dict[str, Any] | None = None
    ) -> tuple[str, dict[str, Any]]:
        with self._lock:
            wid = self._new_id()
            work_dir = Path(os.environ.get("OUTPUT_BASE_DIRECTORY", "results")) / wid
            item = WorkItem(
                id=wid,
                config=config,
                status=WorkStatus.QUEUED,
                created_at=self._now(),
                updated_at=self._now(),
                output_path=str(work_dir),
                extra=extra or {},
            )
            self._repo.add(item)
            return wid, item.to_dict()

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

    def attach_directory(self, work_id: str, directory: Path) -> bool:
        with self._lock:
            item = self._repo.get(work_id)
            if not item:
                return False
            item.output_path = str(directory)
            item.updated_at = self._now()
            self._repo.update(item)
            return True

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
        terminal = {
            WorkStatus.FINISHED.value,
            WorkStatus.ERROR.value,
            WorkStatus.CANCELED.value,
        }

        while True:
            status = self.get_status(work_id)
            if status is None:
                return None
            last_status = status
            if status in terminal:
                return status

            if timeout is not None and (self._now() - start) >= timeout:
                return last_status

            time.sleep(poll_interval)
