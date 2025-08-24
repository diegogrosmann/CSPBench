"""Pipeline combination-scoped persistence wrapper."""

from typing import Any, Optional
from src.infrastructure.persistence.work_state.core import WorkStatePersistence


class CombinationScopedPersistence:
    """
    Wrapper for WorkStatePersistence with work_id and combination identifiers pre-configured.

    Eliminates the need to pass work_id and combination identifiers repeatedly.
    Useful for operations focused on a specific pipeline combination.
    """

    def __init__(self, store: WorkStatePersistence, combination_id: int):
        self._store = store
        self._combination_id = combination_id
        # Campos derivados (serão carregados do banco)
        self._work_id: str | None = None
        self._task_id: str | None = None
        self._dataset_id: str | None = None
        self._preset_id: str | None = None
        self._algorithm_id: str | None = None
        self._mode: str | None = None
        self._load_from_id()

    def _load_from_id(self) -> None:
        """Carrega os campos da combinação a partir do ID."""
        row = self._store._fetch_one(  # uso interno controlado
            """
            SELECT work_id, task_id, dataset_id, preset_id, algorithm_id, mode
            FROM combinations
            WHERE id = ?
            """,
            (self._combination_id,),
        )
        if not row:
            raise ValueError(f"Combinação id={self._combination_id} não encontrada")
        (
            self._work_id,
            self._task_id,
            self._dataset_id,
            self._preset_id,
            self._algorithm_id,
            self._mode,
        ) = row

    @property
    def combination_id(self) -> int:
        return self._combination_id

    @property
    def work_id(self) -> str:
        return self._work_id  # type: ignore[return-value]

    @property
    def task_id(self) -> str:
        return self._task_id  # type: ignore[return-value]

    @property
    def dataset_id(self) -> str:
        return self._dataset_id  # type: ignore[return-value]

    @property
    def preset_id(self) -> str:
        return self._preset_id  # type: ignore[return-value]

    @property
    def algorithm_id(self) -> str:
        return self._algorithm_id  # type: ignore[return-value]

    @property
    def mode(self) -> str:
        return self._mode  # type: ignore[return-value]

    @property
    def store(self) -> WorkStatePersistence:
        return self._store

    def __repr__(self) -> str:
        return (
            "CombinationScopedPersistence("
            f"combination_id={self._combination_id}, work_id={self._work_id}, "
            f"task_id={self._task_id}, dataset_id={self._dataset_id}, "
            f"preset_id={self._preset_id}, algorithm_id={self._algorithm_id}, mode={self._mode}"
            ")"
        )

    def __getattr__(self, name: str):
        """Delegar ao store para métodos/atributos não cobertos pelo wrapper."""
        return getattr(self._store, name)

    # === Helpers internos ===
    def _get_combination(self, required: bool = True) -> Optional[dict[str, Any]]:
        # Busca direta pelo ID (mais eficiente)
        row = self._store._fetch_one(
            """
            SELECT id, work_id, task_id, dataset_id, preset_id, algorithm_id, mode, status, total_sequences,
                   created_at, started_at, finished_at
            FROM combinations
            WHERE id = ?
            """,
            (self._combination_id,),
        )
        if row:
            return {
                "id": row[0],
                "work_id": row[1],
                "task_id": row[2],
                "dataset_id": row[3],
                "preset_id": row[4],
                "algorithm_id": row[5],
                "mode": row[6],
                "status": row[7],
                "total_sequences": row[8],
                "created_at": row[9],
                "started_at": row[10],
                "finished_at": row[11],
            }
        if required:
            raise ValueError("Combinação não encontrada")
        return None

    def _get_combination_id(self, required: bool = True) -> Optional[int]:
        return self._combination_id

    # === Combinação ===
    def update_combination_status(self, status: str) -> None:
        """Atualiza status da combinação atual."""
        self._store.update_combination_status(
            self.work_id,
            self.task_id,
            self.dataset_id,
            self.preset_id,
            self.algorithm_id,
            status,
        )

    def _get_combination(self) -> dict[str, Any]:
        """Retorna os dados da combinação atual."""
        return self._get_combination(required=True)  # type: ignore[return-value]

    def get_combination_progress(self) -> Optional[float]:
        """Retorna progresso da combinação atual."""
        return self._store.get_combination_progress(self._combination_id)

    # === Execuções da combinação ===
    def submit_execution(
        self, *, unit_id: str, sequencia: Optional[int] = None
    ) -> None:
        """Cria uma execução associada à combinação atual."""
        self._store.submit_execution(
            unit_id=unit_id,
            combination_id=self._combination_id,
            sequencia=sequencia,
        )

    def get_executions(
        self,
        *,
        unit_id: Optional[str] = None,
        status: Optional[str] = None,
        sequencia: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """Lista execuções associadas à combinação atual."""
        return self._store.get_executions(
            unit_id=unit_id,
            status=status,
            combination_id=self._combination_id,
            sequencia=sequencia,
        )

    # === Eventos ===
    def combination_warning(
        self, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.combination_warning(
            self.work_id, str(self._combination_id), message, context
        )

    def combination_error(self, error: Exception) -> None:
        self._store.combination_error(self.work_id, str(self._combination_id), error)

    def record_error(self, error: Exception) -> None:
        """Alias for combination_error for backward compatibility."""
        self.combination_error(error)

    def unit_warning(
        self, unit_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.unit_warning(self.work_id, unit_id, message, context)

    def unit_error(self, unit_id: str, error: Exception) -> None:
        self._store.unit_error(self.work_id, unit_id, error)

    def generic_event(
        self,
        unit_id: str,
        event_type: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        self._store.generic_event(self.work_id, unit_id, event_type, message, context)
