"""Unit-scoped persistence wrapper."""

from typing import Any, Optional

from src.infrastructure.persistence.work_state.core import WorkPersistence
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import CombinationScopedPersistence


class ExecutionScopedPersistence(CombinationScopedPersistence):
    """
    Wrapper para WorkPersistence com unit_id preconfigurado.

    Inherits from CombinationScopedPersistence to provide work and combination-level operations
    while adding execution-specific functionality.
    """

    def __init__(self, unit_id: str, work_id: str, store: Optional[WorkPersistence] = None):
        self._unit_id = unit_id
        self._expected_work_id = work_id
        
        # Load combination_id and execution_id from unit_id
        temp_store = store if store is not None else WorkPersistence()
        combination_id, execution_id = self._load_execution_data(temp_store, unit_id, work_id)
        
        # Initialize parent with combination_id
        super().__init__(combination_id, store)
        
        # Verify work_id matches expected
        if self.work_id != work_id:
            raise ValueError(
                f"Work ID mismatch: expected {work_id}, got {self.work_id} for unit {unit_id}"
            )
        
        # Store execution-specific fields
        self._execution_id = execution_id

    def _load_execution_data(self, store: WorkPersistence, unit_id: str, expected_work_id: str) -> tuple[int, int]:
        """Load execution data from database."""
        from src.infrastructure.persistence.work_state.models import Execution, Combination
        
        # Search for execution with unit_id and filter by work_id through combination
        with store.session_scope() as session:
            query = session.query(Execution).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Execution.unit_id == unit_id,
                Combination.work_id == expected_work_id
            )
            
            execution_model = query.first()
            
            if not execution_model:
                raise ValueError(
                    f"Unidade de execução unit_id={unit_id} e work_id={expected_work_id} não encontrada"
                )
            
            execution = execution_model.to_dict()
        
        # Get combination to verify work_id
        combination = store.combination_get(execution['combination_id'])
        if not combination:
            raise ValueError(f"Combination {execution['combination_id']} not found")
            
        if combination['work_id'] != expected_work_id:
            raise ValueError(
                f"Work ID mismatch: expected {expected_work_id}, got {combination['work_id']} for unit {unit_id}"
            )
        
        return execution['combination_id'], execution['id']

    @property
    def unit_id(self) -> str:
        return self._unit_id

    @property
    def execution_id(self) -> int:
        return self._execution_id

    def __repr__(self) -> str:
        return (
            f"ExecutionScopedPersistence(unit_id={self._unit_id}, "
            f"work_id={self.work_id}, combination_id={self.combination_id}, "
            f"execution_id={self._execution_id})"
        )

    # === Execução ===
    def update_execution_status(
        self,
        status: str,
        result: dict[str, Any] | None = None,
        objective: float | None = None,
        params: dict[str, Any] | None = None,
    ) -> None:
        """Atualiza status da execução atual."""
        import time
        
        fields = {'status': status}
        if result is not None:
            fields['result'] = result
        if objective is not None:
            fields['objective'] = objective
        if params is not None:
            fields['params'] = params
        
        # Automatically set timestamps based on status
        current_time = time.time()
        if status == "running":
            fields['started_at'] = current_time
        elif status in ["completed", "failed", "error", "canceled"]:
            fields['finished_at'] = current_time
            
        self.store.execution_update(self._execution_id, **fields)

    def get_execution_info(self) -> dict[str, Any] | None:
        """Retorna informações da execução atual."""
        return self.store.execution_get(self._execution_id)

    # === Eventos ===
    # Override parent methods to use unit-specific context
    def log_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        """Log warning for the current execution unit."""
        self.unit_warning(self._unit_id, message, context)

    def log_error(self, error: Exception) -> None:
        """Log error for the current execution unit."""
        self.unit_error(self._unit_id, error)

    def unit_warning(self, unit_id: str, message: str, context: dict[str, Any] | None = None) -> None:
        """Registra warning para a unidade especificada (ou atual se não especificada)."""
        actual_unit_id = unit_id if unit_id != self._unit_id else self._unit_id
        # Use event_create directly since it's the actual method that exists
        self.store.event_create(
            work_id=self.work_id,
            event_type="warning",
            event_category="unit",
            entity_data={"unit_id": actual_unit_id, "message": message, **(context or {})}
        )

    def unit_error(self, unit_id: str, error: Exception) -> None:
        """Registra erro para a unidade especificada (ou atual se não especificada)."""
        actual_unit_id = unit_id if unit_id != self._unit_id else self._unit_id
        # Use event_create directly since it's the actual method that exists
        self.store.event_create(
            work_id=self.work_id,
            event_type="error",
            event_category="unit",
            entity_data={
                "unit_id": actual_unit_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error)
            }
        )

    def generic_event(
        self, event_type: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Registra evento genérico para a unidade atual."""
        # Use event_create directly since it's the actual method that exists
        self.store.event_create(
            work_id=self.work_id,
            event_type=event_type,
            event_category="unit",
            entity_data={"unit_id": self._unit_id, "message": message, **(context or {})}
        )

    # === Consultas relacionadas ===
    def get_combination_info(self) -> dict[str, Any] | None:
        """Retorna informações da combinação associada."""
        return self.get_combination_data()

    def get_related_executions(self) -> list[dict[str, Any]]:
        """Retorna todas as execuções da mesma combinação."""
        executions, _ = self.store.execution_list(
            filters={'combination_id': self.combination_id}
        )
        return executions

    def get_unit_events(
        self,
        event_type: str | None = None,
        limit: int = 100,
    ) -> list[dict[str, Any]]:
        """Retorna eventos relacionados à unidade atual."""
        # Build filters for the event_list method
        filters = {
            'work_id': self.work_id,
            'event_category': 'unit'
        }
        
        if event_type is not None:
            filters['event_type'] = event_type
            
        # Use event_list which is the actual method that exists
        events, _ = self.store.event_list(
            filters=filters,
            limit=limit,
            order_by='timestamp',
            order_desc=True
        )
        
        # Filter by unit_id from entity_data since it's stored as JSON
        unit_events = []
        for event in events:
            entity_data = event.get('entity_data_json', {})
            if entity_data.get('unit_id') == self._unit_id:
                unit_events.append(event)
                
        return unit_events

    # === Progresso da Execução ===
    def add_progress(
        self,
        progress: float,
        message: str | None = None,
    ) -> None:
        """Adiciona entrada de progresso para a execução atual."""
        self.store.execution_progress_create(execution_id=self._execution_id, progress=progress, message=message)

    def get_progress(self, limit: int | None = None) -> list[dict[str, Any]]:
        """Retorna entradas de progresso da execução atual."""
        filters = {'execution_id': self._execution_id}
        results, _ = self.store.execution_progress_list(
            filters=filters,
            order_by='timestamp',
            order_desc=True,
            limit=limit
        )
        return results
