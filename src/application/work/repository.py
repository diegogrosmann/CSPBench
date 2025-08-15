from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterable, Optional

from src.domain.work import WorkItem


class WorkRepository(ABC):
    """Porta (abstração) para armazenar WorkItems.

    Implementações podem ser in-memory, file-system, banco, etc.
    Cada método deve ser thread-safe se for compartilhado entre threads.
    """

    @abstractmethod
    def add(self, item: WorkItem) -> None:
        """Persiste o novo WorkItem.
        Deve falhar se id já existir (opcional: sobrescrever).
        """
        raise NotImplementedError

    @abstractmethod
    def get(self, work_id: str) -> Optional[WorkItem]:
        """Retorna WorkItem ou None se não encontrado."""
        raise NotImplementedError

    @abstractmethod
    def list(self) -> Iterable[WorkItem]:
        """Lista todos os WorkItems."""
        raise NotImplementedError

    @abstractmethod
    def update(self, item: WorkItem) -> None:
        """Atualiza WorkItem existente (no-op se não existir)."""
        raise NotImplementedError

    @abstractmethod
    def remove(self, work_id: str) -> None:
        """Remove WorkItem se existir (idempotente)."""
        raise NotImplementedError


class InMemoryWorkRepository(WorkRepository):
    def __init__(self):
        self._items: dict[str, WorkItem] = {}

    def add(self, item: WorkItem) -> None:
        self._items[item.id] = item

    def get(self, work_id: str) -> Optional[WorkItem]:
        return self._items.get(work_id)

    def list(self):
        return list(self._items.values())

    def update(self, item: WorkItem) -> None:
        # para in-memory nada adicional
        self._items[item.id] = item

    def remove(self, work_id: str) -> None:
        self._items.pop(work_id, None)
