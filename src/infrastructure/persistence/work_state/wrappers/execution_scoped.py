"""Unit-scoped persistence wrapper."""

from typing import Any, Optional
from src.infrastructure.persistence.work_state.core import WorkStatePersistence


class ExecutionScopedPersistence:
    """
    Wrapper para WorkStatePersistence com unit_id preconfigurado.

    Carrega automaticamente work_id e combination_id a partir do unit_id.
    Útil para operações focadas em uma unidade de execução específica.
    """

    def __init__(self, store: WorkStatePersistence, unit_id: str):
        self._store = store
        self._unit_id = unit_id
        # Campos derivados (serão carregados do banco)
        self._work_id: str | None = None
        self._combination_id: int | None = None
        self._load_from_unit_id()

    def _load_from_unit_id(self) -> None:
        """Carrega work_id e combination_id a partir do unit_id."""
        row = self._store._fetch_one(
            """
            SELECT c.work_id, e.combination_id
            FROM executions e
            JOIN combinations c ON e.combination_id = c.id
            WHERE e.unit_id = ?
            """,
            (self._unit_id,),
        )
        if not row:
            raise ValueError(
                f"Unidade de execução unit_id={self._unit_id} não encontrada"
            )
        self._work_id, self._combination_id = row

    @property
    def unit_id(self) -> str:
        return self._unit_id

    @property
    def work_id(self) -> str:
        return self._work_id  # type: ignore[return-value]

    @property
    def combination_id(self) -> int:
        return self._combination_id  # type: ignore[return-value]

    @property
    def store(self) -> WorkStatePersistence:
        return self._store

    def __repr__(self) -> str:
        return (
            f"ExecutionScopedPersistence(unit_id={self._unit_id}, "
            f"work_id={self._work_id}, combination_id={self._combination_id})"
        )

    def __getattr__(self, name: str):
        """Delegar ao store para métodos/atributos não cobertos pelo wrapper."""
        return getattr(self._store, name)

    # === Execução ===
    def update_execution_status(
        self,
        status: str,
        result: dict[str, Any] | None = None,
        objective: float | None = None,
        params: dict[str, Any] | None = None,
    ) -> None:
        """Atualiza status da execução atual."""
        self._store.update_execution_status(
            self._unit_id, status, result=result, objective=objective, params=params
        )

    def get_execution_info(self) -> dict[str, Any] | None:
        """Retorna informações da execução atual."""
        executions = self._store.get_executions(unit_id=self._unit_id)
        return executions[0] if executions else None

    # === Eventos ===
    def unit_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        """Registra warning para a unidade atual."""
        self._store.unit_warning(self.work_id, self._unit_id, message, context)

    def unit_error(self, error: Exception) -> None:
        """Registra erro para a unidade atual."""
        self._store.unit_error(self.work_id, self._unit_id, error)

    def generic_event(
        self, event_type: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Registra evento genérico para a unidade atual."""
        self._store.generic_event(
            self.work_id, self._unit_id, event_type, message, context
        )

    # === Consultas relacionadas ===
    def get_combination_info(self) -> dict[str, Any] | None:
        """Retorna informações da combinação associada."""
        combinations = self._store.get_combinations(work_id=self.work_id)
        for combo in combinations:
            if combo["id"] == self.combination_id:
                return combo
        return None

    def get_related_executions(self) -> list[dict[str, Any]]:
        """Retorna todas as execuções da mesma combinação."""
        return self._store.get_executions(combination_id=self.combination_id)

    def get_unit_events(
        self,
        event_type: str | None = None,
        limit: int = 100,
    ) -> list[dict[str, Any]]:
        """Retorna eventos relacionados à unidade atual."""
        return self._store.get_events(
            self.work_id,
            event_category="unit",
            event_type=event_type,
            unit_id=self._unit_id,
            limit=limit,
        )

    # === Progresso da Execução ===
    def add_progress(
        self,
        progress: float,
        message: str | None = None,
    ) -> None:
        """Adiciona entrada de progresso para a execução atual."""
        self._store.add_execution_progress(self._unit_id, progress, message)

    def get_progress(self, limit: int | None = None) -> list[dict[str, Any]]:
        """Retorna entradas de progresso da execução atual."""
        return self._store.get_execution_progress(self._unit_id, limit)
