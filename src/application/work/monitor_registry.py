"""Registro simples de monitores por work_id.

Fornece API thread-safe para registrar/recuperar/desregistrar monitores associados a um work_id.
"""
from __future__ import annotations

import threading
from typing import Any, Dict, Optional

from src.infrastructure.monitoring.monitor_interface import NoOpMonitor


class MonitorRegistry:
    """Singleton-like registry (módulo) para monitores em memória."""

    def __init__(self):
        self._lock = threading.RLock()
        self._registry: Dict[str, Any] = {}

    def register(self, work_id: str, monitor: Any) -> None:
        with self._lock:
            self._registry[work_id] = monitor

    def get(self, work_id: str) -> Optional[Any]:
        with self._lock:
            return self._registry.get(work_id)

    def unregister(self, work_id: str) -> None:
        with self._lock:
            self._registry.pop(work_id, None)

    def list_ids(self) -> list[str]:
        with self._lock:
            return list(self._registry.keys())


# Módulo export: instância compartilhada
registry = MonitorRegistry()

__all__ = ["registry"]
