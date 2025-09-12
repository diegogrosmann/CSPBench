"""Pipeline combination-scoped persistence wrapper."""

from typing import Any, Optional, TYPE_CHECKING

from src.infrastructure.persistence.work_state.core import WorkPersistence
from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence

if TYPE_CHECKING:
    from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence


class CombinationScopedPersistence(WorkScopedPersistence):
    """
    Wrapper for WorkPersistence with work_id and combination identifiers pre-configured.

    Inherits from WorkScopedPersistence to provide work-level operations
    while adding combination-specific functionality.
    """

    def __init__(self, combination_id: int, store: Optional[WorkPersistence] = None):
        self._combination_id = combination_id
        
        # Load work_id and other fields from combination_id
        temp_store = store if store is not None else WorkPersistence()
        work_id, task_id, dataset_id, preset_id, algorithm_id, mode = self._load_combination_data(temp_store)
        
        # Initialize parent with work_id
        super().__init__(work_id, store)
        
        # Store combination-specific fields
        self._task_id = task_id
        self._dataset_id = dataset_id
        self._preset_id = preset_id
        self._algorithm_id = algorithm_id
        self._mode = mode

    def _load_combination_data(self, store: WorkPersistence) -> tuple[str, str, str, str, str, str]:
        """Load combination data from database."""
        combination = store.combination_get(self._combination_id)
        if not combination:
            raise ValueError(f"Combinação id={self._combination_id} não encontrada")
        
        return (
            combination['work_id'],
            combination['task_id'],
            combination['dataset_id'], 
            combination['preset_id'],
            combination['algorithm_id'],
            combination['mode']
        )

    @property
    def combination_id(self) -> int:
        return self._combination_id

    @property
    def task_id(self) -> str:
        return self._task_id

    @property
    def dataset_id(self) -> str:
        return self._dataset_id

    @property
    def preset_id(self) -> str:
        return self._preset_id

    @property
    def algorithm_id(self) -> str:
        return self._algorithm_id

    @property
    def mode(self) -> str:
        return self._mode

    def __repr__(self) -> str:
        return (
            "CombinationScopedPersistence("
            f"combination_id={self._combination_id}, work_id={self.work_id}, "
            f"task_id={self._task_id}, dataset_id={self._dataset_id}, "
            f"preset_id={self._preset_id}, algorithm_id={self._algorithm_id}, mode={self._mode}"
            ")"
        )

    # === Helpers internos ===
    def _get_combination(self, required: bool = True) -> Optional[dict[str, Any]]:
        """Get combination data by ID."""
        combination = self.store.combination_get(self._combination_id)
        if not combination and required:
            raise ValueError("Combinação não encontrada")
        return combination

    def _get_combination_id(self, required: bool = True) -> Optional[int]:
        return self._combination_id

    # === Combinação ===
    def update_combination_status(self, status: str) -> None:
        """Atualiza status da combinação atual."""
        import time
        
        fields = {'status': status}
        
        # Automatically set timestamps based on status
        current_time = time.time()
        if status == "running":
            fields['started_at'] = current_time
        elif status in ["completed", "failed", "error", "canceled"]:
            fields['finished_at'] = current_time
            
        self._store.combination_update(self._combination_id, **fields)

    def get_combination_data(self) -> dict[str, Any]:
        """Retorna os dados da combinação atual."""
        return self._get_combination(required=True)  # type: ignore[return-value]

    def get_combination_progress(self) -> Optional[float]:
        """Retorna progresso da combinação atual."""
        # Get all executions for this combination
        executions, _ = self.store.execution_list(
            filters={'combination_id': self._combination_id}
        )
        
        if not executions:
            return 0.0
            
        completed = sum(1 for ex in executions if ex['status'] in ['completed', 'failed', 'error'])
        total = len(executions)
        
        return completed / total if total > 0 else 0.0

    # === Execuções da combinação ===
    def submit_execution(
        self, *, unit_id: str, sequencia: int
    ) -> None:
        """Cria uma execução associada à combinação atual."""
        self.store.execution_create(
            unit_id=unit_id,
            combination_id=self._combination_id,
            sequencia=sequencia,
        )

    def submit_executions(self, executions: list[dict[str, Any]]) -> int:
        """Submete execuções para a combinação atual usando inserção em lote."""
        if not executions:
            return 0
        
        # Adiciona combination_id para todas as execuções
        executions_with_combination_id = []
        for exec_data in executions:
            exec_with_combination_id = exec_data.copy()
            exec_with_combination_id['combination_id'] = self._combination_id
            executions_with_combination_id.append(exec_with_combination_id)
        
        # Usa inserção em lote para melhor performance
        return self.store.execution_bulk_create(executions_with_combination_id)

    def get_executions(
        self,
        *,
        unit_id: Optional[str] = None,
        status: Optional[str] = None,
        sequencia: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """Lista execuções associadas à combinação atual."""
        filters = {'combination_id': self._combination_id}
        if unit_id is not None:
            filters['unit_id'] = unit_id
        if status is not None:
            filters['status'] = status
        if sequencia is not None:
            filters['sequencia'] = sequencia
            
        executions, _ = self.store.execution_list(filters=filters)
        return executions

    # === Eventos ===
    # Override parent methods to use combination-specific context
    def log_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        """Log warning for the current combination."""
        self.combination_warning(message, context)

    def log_error(self, error: Exception) -> None:
        """Log error for the current combination."""
        self.combination_error(error)

    def combination_warning(
        self, message: str, context: dict[str, Any] | None = None
    ) -> None:
        # Use parent class's method for combination warnings
        super().combination_warning(str(self._combination_id), message, context)

    def combination_error(self, error: Exception) -> None:
        # Use parent class's method for combination errors  
        super().combination_error(str(self._combination_id), error)

    def record_error(self, error: Exception) -> None:
        """Alias for combination_error for backward compatibility."""
        self.combination_error(error)

    def unit_warning(
        self, unit_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        # Use parent class's method for unit warnings
        super().unit_warning(unit_id, message, context)

    def unit_error(self, unit_id: str, error: Exception) -> None:
        # Use parent class's method for unit errors
        super().unit_error(unit_id, error)

    def generic_event(
        self,
        unit_id: str,
        event_type: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        # Use the parent class's generic_event method from WorkScopedPersistence
        super().generic_event(unit_id, event_type, message, context)

    # === FACTORY METHODS ===
    def for_execution(self, unit_id: str) -> "ExecutionScopedPersistence":
        """
        Create an ExecutionScopedPersistence for the specified execution unit within this combination.
        
        Args:
            unit_id: ID of the execution unit
            
        Returns:
            ExecutionScopedPersistence instance
            
        Raises:
            ValueError: If execution doesn't belong to this combination
        """
        from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence
        
        # Verify execution belongs to this combination
        executions = self.get_executions()
        unit_ids = [exec["unit_id"] for exec in executions]
        
        if unit_id not in unit_ids:
            raise ValueError(
                f"Execution unit {unit_id} not found in combination {self._combination_id}. "
                f"Available units: {unit_ids}"
            )
        
        return ExecutionScopedPersistence(unit_id, self.work_id, self.store)

    def get_all_executions(self) -> list["ExecutionScopedPersistence"]:
        """
        Get all executions in this combination as ExecutionScopedPersistence instances.
        
        Returns:
            List of ExecutionScopedPersistence instances
        """
        from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence
        
        executions = self.get_executions()
        return [ExecutionScopedPersistence(exec["unit_id"], self.work_id, self.store) for exec in executions]

    # === Limpeza de Progresso ===
    def clear_progress_for_non_finalized_executions(self) -> int:
        """
        Limpa entradas de progresso de todas as execuções não finalizadas desta combinação.
        
        Returns:
            Número de entradas de progresso removidas
        """
        return self.store.execution_progress_clear_for_non_finalized(
            combination_id=self._combination_id
        )
