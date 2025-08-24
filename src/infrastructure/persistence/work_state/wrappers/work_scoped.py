"""Work-scoped persistence wrapper."""

from typing import Any, Optional
from src.infrastructure.persistence.work_state.core import WorkStatePersistence


class WorkScopedPersistence:
    """
    Wrapper para WorkStatePersistence com work_id preconfigurado.

    Evita repetir work_id em cada chamada, facilitando operações
    focadas em um trabalho específico.
    """

    def __init__(self, store: WorkStatePersistence, work_id: str):
        self._store = store
        self._work_id = work_id

    @property
    def work_id(self) -> str:
        return self._work_id

    @property
    def store(self) -> WorkStatePersistence:
        return self._store

    def __repr__(self) -> str:
        return f"WorkScopedPersistence(work_id={self._work_id})"

    # Permite acessar métodos/atributos não sobrescritos diretamente no store
    def __getattr__(self, name: str):
        return getattr(self._store, name)

    # === WORK ===
    def update_work_status(self, status: str, **fields: Any) -> None:
        """Atualiza status do work atual."""
        self._store.update_work_status(self._work_id, status, **fields)

    def get_work_status(self) -> str | None:
        """Obtém status do work atual do WorkService."""
        return self._store.get_work_status(self._work_id)

    def get_work_statistics(self) -> dict[str, Any]:
        """Obtém estatísticas do work atual."""
        return self._store.get_work_statistics(self._work_id)

    def algorithm_error(self, algorithm_id: str, error: Exception) -> None:
        """Log algorithm error for the current work."""
        from src.infrastructure.persistence.work_state.events import EventsMixin
        if isinstance(self._store, EventsMixin):
            # Use generic event logging for algorithm errors
            data = {
                "algorithm_id": algorithm_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            }
            self._store._log_event(self._work_id, "error", "other", data)

    # === COMBINATIONS ===
    def submit_combinations(self, tasks_combinations: list[dict[str, Any]]) -> int:
        """Submete combinações para o work atual."""
        return self._store.submit_combinations(self._work_id, tasks_combinations)

    def update_combination_status(
        self,
        task_id: str,
        dataset_id: str,
        preset_id: str,
        algorithm_id: str,
        status: str,
    ) -> None:
        """Atualiza o status de uma combinação do work atual."""
        self._store.update_combination_status(
            self._work_id, task_id, dataset_id, preset_id, algorithm_id, status
        )

    def get_next_queued_combination(self) -> Optional[dict[str, Any]]:
        """Obtém a próxima combinação em fila do work atual."""
        return self._store.get_next_queued_combination(self._work_id)

    def get_next_pending_combination(self) -> Optional[dict[str, Any]]:
        """Compat: próxima combinação pendente do work atual."""
        return self._store.get_next_pending_combination(self._work_id)

    def get_combinations(
        self,
        *,
        status: Optional[str] = None,
        task_id: Optional[str] = None,
        dataset_id: Optional[str] = None,
        preset_id: Optional[str] = None,
        algorithm_id: Optional[str] = None,
        mode: Optional[str] = None,
        order_by: str = "created_at",
        order_direction: str = "ASC",
        limit: Optional[int] = None,
        offset: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """Lista combinações do work atual com filtros opcionais."""
        return self._store.get_combinations(
            work_id=self._work_id,
            status=status,
            task_id=task_id,
            dataset_id=dataset_id,
            preset_id=preset_id,
            algorithm_id=algorithm_id,
            mode=mode,
            order_by=order_by,
            order_direction=order_direction,
            limit=limit,
            offset=offset,
        )

    def init_combination(self) -> bool:
        """Reinicia combinações em andamento/pausadas/canceladas do work atual para 'queued'."""
        return self._store.init_combination(self._work_id)

    # === EVENTS ===
    # Work
    def work_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        self._store.work_warning(self._work_id, message, context)

    def work_error(self, error: Exception) -> None:
        self._store.work_error(self._work_id, error)

    # Task
    def task_warning(
        self, task_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.task_warning(self._work_id, task_id, message, context)

    def task_error(self, task_id: str, error: Exception) -> None:
        self._store.task_error(self._work_id, task_id, error)

    # Dataset
    def dataset_warning(
        self, dataset_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.dataset_warning(self._work_id, dataset_id, message, context)

    def dataset_error(self, dataset_id: str, error: Exception) -> None:
        self._store.dataset_error(self._work_id, dataset_id, error)

    # Preset
    def preset_warning(
        self, preset_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.preset_warning(self._work_id, preset_id, message, context)

    def preset_error(self, preset_id: str, error: Exception) -> None:
        self._store.preset_error(self._work_id, preset_id, error)

    # Combination
    def combination_warning(
        self, combination_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.combination_warning(self._work_id, combination_id, message, context)

    def combination_error(self, combination_id: str, error: Exception) -> None:
        self._store.combination_error(self._work_id, combination_id, error)

    # Unit
    def unit_warning(
        self, unit_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.unit_warning(self._work_id, unit_id, message, context)

    def unit_error(self, unit_id: str, error: Exception) -> None:
        self._store.unit_error(self._work_id, unit_id, error)

    # Genérico
    def generic_event(
        self,
        unit_id: str,
        event_type: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        self._store.generic_event(self._work_id, unit_id, event_type, message, context)

    # Consultas de eventos
    def get_events(
        self,
        *,
        event_category: str | None = None,
        event_type: str | None = None,
        unit_id: str | None = None,
        limit: int = 100,
    ) -> list[dict[str, Any]]:
        return self._store.get_events(
            self._work_id,
            event_category=event_category,
            event_type=event_type,
            unit_id=unit_id,
            limit=limit,
        )

    def get_events_by_category(
        self, event_category: str, limit: int = 100
    ) -> list[dict[str, Any]]:
        return self._store.get_events_by_category(self._work_id, event_category, limit)

    def get_warnings(
        self, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        return self._store.get_warnings(
            self._work_id, event_category=event_category, limit=limit
        )

    def get_errors(
        self, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        return self._store.get_errors(
            self._work_id, event_category=event_category, limit=limit
        )

    def get_event_summary_by_category(self) -> dict[str, dict[str, int]]:
        return self._store.get_event_summary_by_category(self._work_id)
