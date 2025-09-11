"""Work-scoped persistence wrapper."""

from typing import Any, Optional, TYPE_CHECKING

from src.infrastructure.persistence.work_state.core import WorkPersistence

if TYPE_CHECKING:
    from src.infrastructure.persistence.work_state.wrappers.combination_scoped import CombinationScopedPersistence
    from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence


class WorkScopedPersistence:
    """
    Wrapper para WorkPersistence com work_id preconfigurado.

    Evita repetir work_id em cada chamada, facilitando operações
    focadas em um trabalho específico.
    """

    def __init__(self, work_id: str, store: Optional[WorkPersistence] = None):
        self._store = store if store is not None else WorkPersistence()
        self._work_id = work_id

    @staticmethod
    def submit(work_id: str) -> "WorkScopedPersistence":
        """
        Método estático para criar um WorkScopedPersistence com store padrão.
        
        Args:
            work_id: ID do trabalho
            
        Returns:
            Nova instância de WorkScopedPersistence
        """
        return WorkScopedPersistence(work_id)

    @property
    def work_id(self) -> str:
        return self._work_id

    @property
    def store(self) -> WorkPersistence:
        return self._store

    def __repr__(self) -> str:
        return f"WorkScopedPersistence(work_id={self._work_id})"

    # Permite acessar métodos/atributos não sobrescritos diretamente no store
    def __getattr__(self, name: str):
        return getattr(self._store, name)

    # === WORK ===
    def update_work_status(self, status: str, **fields: Any) -> None:
        """Atualiza status do work atual."""
        self._store.work_update(self._work_id, status=status, **fields)

    def get_work_status(self) -> str | None:
        """Obtém status do work atual."""
        work = self._store.work_get(self._work_id)
        return work['status'] if work else None

    def get_work_data(self) -> dict[str, Any] | None:
        """Obtém dados completos do work atual."""
        return self._store.work_get(self._work_id)

    def algorithm_error(self, algorithm_id: str, error: Exception) -> None:
        """Log algorithm error for the current work."""
        # Use event logging for algorithm errors
        self._store.event_create(
            work_id=self._work_id,
            event_type="error", 
            event_category="algorithm",
            entity_data={
                "algorithm_id": algorithm_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            }
        )

    # === COMBINATIONS ===
    def submit_combinations(self, tasks_combinations: list[dict[str, Any]]) -> int:
        """Submete combinações para o work atual."""
        count = 0
        for combo in tasks_combinations:
            self._store.combination_create(
                work_id=self._work_id,
                **combo
            )
            count += 1
        return count

    def update_combination_status(
        self,
        combination_id: int,
        status: str,
    ) -> None:
        """Atualiza o status de uma combinação do work atual."""
        import time
        
        fields = {'status': status}
        
        # Automatically set timestamps based on status
        current_time = time.time()
        if status == "running":
            fields['started_at'] = current_time
        elif status in ["completed", "failed", "error", "canceled"]:
            fields['finished_at'] = current_time
            
        self._store.combination_update(combination_id, **fields)

    def get_next_queued_combination(self) -> Optional[dict[str, Any]]:
        """Obtém a próxima combinação em fila do work atual."""
        combinations, _ = self._store.combination_list(
            filters={'work_id': self._work_id, 'status': 'queued'},
            order_by='created_at',
            limit=1
        )
        return combinations[0] if combinations else None

    def get_next_pending_combination(self) -> Optional[dict[str, Any]]:
        """Compat: próxima combinação pendente do work atual."""
        return self.get_next_queued_combination()

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
        # Build filters dict
        filters = {'work_id': self._work_id}
        if status is not None:
            filters['status'] = status
        if task_id is not None:
            filters['task_id'] = task_id
        if dataset_id is not None:
            filters['dataset_id'] = dataset_id
        if preset_id is not None:
            filters['preset_id'] = preset_id
        if algorithm_id is not None:
            filters['algorithm_id'] = algorithm_id
        if mode is not None:
            filters['mode'] = mode
            
        order_desc = (order_direction.upper() == "DESC")
        
        results, _ = self._store.combination_list(
            filters=filters,
            order_by=order_by,
            order_desc=order_desc,
            limit=limit,
            offset=offset or 0
        )
        return results

    def init_combination(self) -> bool:
        """Reinicia combinações em andamento/pausadas/canceladas do work atual para 'queued'."""
        # Get combinations to reset
        combinations, _ = self._store.combination_list(
            filters={
                'work_id': self._work_id,
                'status': ['running', 'paused', 'canceled']
            }
        )
        
        # Reset their status to queued
        for combo in combinations:
            self._store.combination_update(combo['id'], status='queued')
            
        return len(combinations) > 0

    # === EVENTS ===
    # Work-level events
    def log_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        """Log warning for the current work."""
        self.work_warning(message, context)

    def log_error(self, error: Exception) -> None:
        """Log error for the current work."""
        self.work_error(error)

    # Work
    def work_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="work",
            entity_data={"message": message, **(context or {})}
        )

    def work_error(self, error: Exception) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="work",
            entity_data={
                "error_type": error.__class__.__name__,
                "error_message": str(error)
            }
        )

    # Task
    def task_warning(
        self, task_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="task",
            entity_data={"task_id": task_id, "message": message, **(context or {})}
        )

    def task_error(self, task_id: str, error: Exception) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="task",
            entity_data={
                "task_id": task_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error)
            }
        )

    # Dataset
    def dataset_warning(
        self, dataset_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="dataset",
            entity_data={"dataset_id": dataset_id, "message": message, **(context or {})}
        )

    def dataset_error(self, dataset_id: str, error: Exception) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="dataset",
            entity_data={
                "dataset_id": dataset_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error)
            }
        )

    def submit_dataset(self, id: str, dataset_obj, meta: dict[str, Any] | None = None) -> None:
        """Submit a dataset for the current work."""
        # Create dataset entry with auto-generated ID
        dataset_auto_id = self._store.dataset_create(
            dataset_id=id,
            work_id=self._work_id,
            name=getattr(dataset_obj, 'name', None),
            meta=meta
        )
        
        # Store sequences if the dataset has them
        if hasattr(dataset_obj, 'sequences'):
            for seq_index, sequence in enumerate(dataset_obj.sequences):
                self._store.dataset_sequence_create(
                    dataset_id=dataset_auto_id,
                    seq_index=seq_index,
                    sequence=sequence
                )

    def get_dataset(self, dataset_id: str):
        """Get a dataset for the current work."""
        # Find dataset by dataset_id and work_id
        datasets, _ = self._store.dataset_list(
            filters={'dataset_id': dataset_id, 'work_id': self._work_id},
            limit=1
        )
        if not datasets:
            return None
            
        dataset = datasets[0]
        auto_id = dataset['id']
        
        # Get sequences
        sequences, _ = self._store.dataset_sequence_list(
            filters={'dataset_id': auto_id},
            order_by='seq_index'
        )
        
        # Return a simple object with the data
        class DatasetResult:
            def __init__(self, data, sequences):
                self.name = data.get('name')
                self.meta = data.get('meta', {})
                self.sequences = [seq['sequence'] for seq in sequences]
                
        return DatasetResult(dataset, sequences)

    # Preset
    def preset_warning(
        self, preset_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="preset",
            entity_data={"preset_id": preset_id, "message": message, **(context or {})}
        )

    def preset_error(self, preset_id: str, error: Exception) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="preset",
            entity_data={
                "preset_id": preset_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error)
            }
        )

    # Combination
    def combination_warning(
        self, combination_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="combination",
            entity_data={"combination_id": combination_id, "message": message, **(context or {})}
        )

    def combination_error(self, combination_id: str, error: Exception) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="combination",
            entity_data={
                "combination_id": combination_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error)
            }
        )

    # Unit
    def unit_warning(
        self, unit_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="unit",
            entity_data={"unit_id": unit_id, "message": message, **(context or {})}
        )

    def unit_error(self, unit_id: str, error: Exception) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="unit",
            entity_data={
                "unit_id": unit_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error)
            }
        )

    # Genérico
    def generic_event(
        self,
        unit_id: str,
        event_type: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        self._store.event_create(
            work_id=self._work_id,
            event_type=event_type,
            event_category="unit",
            entity_data={"unit_id": unit_id, "message": message, **(context or {})}
        )

    # Consultas de eventos
    def get_events(
        self,
        *,
        event_category: str | None = None,
        event_type: str | None = None,
        unit_id: str | None = None,
        limit: int = 100,
    ) -> list[dict[str, Any]]:
        filters = {'work_id': self._work_id}
        if event_category:
            filters['event_category'] = event_category
        if event_type:
            filters['event_type'] = event_type
        if unit_id:
            # For unit_id, we need to check the entity_data JSON field
            pass  # We'll need to implement custom filtering for JSON fields
        
        events, _ = self._store.event_list(
            filters=filters,
            order_by='timestamp',
            order_desc=True,
            limit=limit
        )
        return events

    def get_events_by_category(
        self, event_category: str, limit: int = 100
    ) -> list[dict[str, Any]]:
        return self.get_events(event_category=event_category, limit=limit)

    def get_warnings(
        self, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        return self.get_events(event_type="warning", event_category=event_category, limit=limit)

    def get_errors(
        self, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        return self.get_events(event_type="error", event_category=event_category, limit=limit)

    def get_event_summary_by_category(self) -> dict[str, dict[str, int]]:
        # Get all events for this work
        events, _ = self._store.event_list(
            filters={'work_id': self._work_id},
            limit=None  # Get all events for summary
        )
        
        # Group by category and type
        summary = {}
        for event in events:
            category = event['event_category']
            event_type = event['event_type']
            
            if category not in summary:
                summary[category] = {}
            if event_type not in summary[category]:
                summary[category][event_type] = 0
            summary[category][event_type] += 1
            
        return summary

    # === FACTORY METHODS ===
    def for_combination(self, combination_id: int) -> "CombinationScopedPersistence":
        """
        Create a CombinationScopedPersistence for the specified combination within this work.
        
        Args:
            combination_id: ID of the combination
            
        Returns:
            CombinationScopedPersistence instance
            
        Raises:
            ValueError: If combination doesn't belong to this work
        """
        from src.infrastructure.persistence.work_state.wrappers.combination_scoped import CombinationScopedPersistence
        
        # Verify combination belongs to this work
        combinations = self.get_combinations()
        combination_ids = [combo["id"] for combo in combinations]
        
        if combination_id not in combination_ids:
            raise ValueError(
                f"Combination {combination_id} not found in work {self._work_id}. "
                f"Available combinations: {combination_ids}"
            )
        
        return CombinationScopedPersistence(combination_id, self._store)

    def for_execution(self, unit_id: str) -> "ExecutionScopedPersistence":
        """
        Create an ExecutionScopedPersistence for the specified execution unit within this work.
        
        Args:
            unit_id: ID of the execution unit
            
        Returns:
            ExecutionScopedPersistence instance
            
        Raises:
            ValueError: If execution doesn't belong to this work
        """
        from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence
        
        return ExecutionScopedPersistence(unit_id, self._work_id, self._store)

    def get_all_combinations(self) -> list["CombinationScopedPersistence"]:
        """
        Get all combinations in this work as CombinationScopedPersistence instances.
        
        Returns:
            List of CombinationScopedPersistence instances
        """
        combinations = self.get_combinations()
        return [self.for_combination(combo["id"]) for combo in combinations]

    def get_all_executions(self) -> list["ExecutionScopedPersistence"]:
        """
        Get all executions in this work as ExecutionScopedPersistence instances.
        
        Returns:
            List of ExecutionScopedPersistence instances
        """
        from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence
        
        # Get all executions across all combinations in this work
        all_executions = []
        for combination in self.get_combinations():
            combination_id = combination["id"]
            combo_scoped = self.for_combination(combination_id)
            executions = combo_scoped.get_executions()
            for execution in executions:
                unit_id = execution["unit_id"]
                all_executions.append(ExecutionScopedPersistence(unit_id, self._work_id, self._store))
        
        return all_executions

    def get_running_executions(self) -> list[dict[str, Any]]:
        """
        Get all executions that are currently in 'running' status.
        
        Returns:
            List of execution dictionaries for executions in running status
        """
        running_executions = []
        
        # Get all executions across all combinations in this work
        for combination in self.get_combinations():
            combination_id = combination["id"]
            combo_scoped = self.for_combination(combination_id)
            executions = combo_scoped.get_executions()
            
            for execution in executions:
                if execution.get("status") == "running":
                    running_executions.append(execution)
        
        return running_executions
