"""Monitoring interface (lean) for pipeline orchestration.

Provides a very small surface so we can plug different sinks (logging, web, memory).
"""

from __future__ import annotations
from typing import Protocol, Optional, Dict, Any
import time


class Monitor(Protocol):  # pragma: no cover - structural type
    def pipeline_started(self, name: str, total_tasks: int) -> None: ...
    def pipeline_finished(
        self, success: bool, summary: Optional[Dict[str, Any]] = None
    ) -> None: ...
    def task_started(self, task_id: str, meta: Dict[str, Any]) -> None: ...
    def task_finished(self, task_id: str, summary: Dict[str, Any]) -> None: ...
    def unit_started(self, unit_id: str, meta: Dict[str, Any]) -> None: ...
    def unit_progress(
        self, unit_id: str, progress: float, message: Optional[str] = None
    ) -> None: ...
    def unit_finished(self, unit_id: str, result: Dict[str, Any]) -> None: ...
    def log(
        self, level: str, message: str, ctx: Optional[Dict[str, Any]] = None
    ) -> None: ...
    def error(
        self,
        unit_id: Optional[str],
        exc: Exception,
        ctx: Optional[Dict[str, Any]] = None,
    ) -> None: ...


class NoOpMonitor:
    """Default monitor that does nothing (safe fallback)."""

    def pipeline_started(self, name: str, total_tasks: int) -> None:
        pass

    def pipeline_finished(
        self, success: bool, summary: Optional[Dict[str, Any]] = None
    ) -> None:
        pass

    def task_started(self, task_id: str, meta: Dict[str, Any]) -> None:
        pass

    def task_finished(self, task_id: str, summary: Dict[str, Any]) -> None:
        pass

    def unit_started(self, unit_id: str, meta: Dict[str, Any]) -> None:
        pass

    def unit_progress(
        self, unit_id: str, progress: float, message: Optional[str] = None
    ) -> None:
        pass

    def unit_finished(self, unit_id: str, result: Dict[str, Any]) -> None:
        pass

    def log(
        self, level: str, message: str, ctx: Optional[Dict[str, Any]] = None
    ) -> None:
        pass

    def error(
        self,
        unit_id: Optional[str],
        exc: Exception,
        ctx: Optional[Dict[str, Any]] = None,
    ) -> None:
        pass


class LoggingMonitor(NoOpMonitor):
    """Simple logging based monitor."""

    def __init__(self, logger):
        self._logger = logger
        self._t0 = time.time()

    def pipeline_started(self, name: str, total_tasks: int) -> None:
        self._logger.info(f"[PIPELINE] started name={name} total_tasks={total_tasks}")

    def pipeline_finished(
        self, success: bool, summary: Optional[Dict[str, Any]] = None
    ) -> None:
        dt = time.time() - self._t0
        self._logger.info(
            f"[PIPELINE] finished success={success} elapsed={dt:.2f}s summary={summary}"
        )

    def task_started(self, task_id: str, meta: Dict[str, Any]) -> None:
        self._logger.info(f"[TASK] start id={task_id} meta={meta}")

    def task_finished(self, task_id: str, summary: Dict[str, Any]) -> None:
        self._logger.info(f"[TASK] finish id={task_id} summary={summary}")

    def unit_started(self, unit_id: str, meta: Dict[str, Any]) -> None:
        self._logger.debug(f"[UNIT] start id={unit_id} meta={meta}")

    def unit_progress(
        self, unit_id: str, progress: float, message: Optional[str] = None
    ) -> None:
        self._logger.debug(
            f"[UNIT] progress id={unit_id} p={progress:.2%} msg={message}"
        )

    def unit_finished(self, unit_id: str, result: Dict[str, Any]) -> None:
        self._logger.debug(f"[UNIT] finish id={unit_id} result={result}")

    def log(
        self, level: str, message: str, ctx: Optional[Dict[str, Any]] = None
    ) -> None:
        log_fn = getattr(self._logger, level, self._logger.info)
        log_fn(f"[LOG] {message} ctx={ctx}")

    def error(
        self,
        unit_id: Optional[str],
        exc: Exception,
        ctx: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._logger.error(f"[ERROR] unit={unit_id} exc={exc} ctx={ctx}")
