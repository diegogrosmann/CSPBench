"""Work-scoped persistence wrapper."""

from typing import TYPE_CHECKING, Any, Optional

from src.domain.dataset import Dataset
from src.infrastructure.persistence.work_state.core import WorkPersistence

if TYPE_CHECKING:
    from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
        CombinationScopedPersistence,
    )
    from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
        ExecutionScopedPersistence,
    )


class WorkScopedPersistence:
    """Wrapper for WorkPersistence with pre-configured work_id.

    Avoids repeating work_id in each call, facilitating operations
    focused on a specific work.

    Args:
        work_id: The work identifier
        store: Optional WorkPersistence instance. If None, creates a new one.
    """

    def __init__(self, work_id: str, store: Optional[WorkPersistence] = None):
        """Initialize WorkScopedPersistence.

        Args:
            work_id: The work identifier
            store: Optional WorkPersistence instance
        """
        self._store = store if store is not None else WorkPersistence()
        self._work_id = work_id

    @staticmethod
    def submit(work_id: str) -> "WorkScopedPersistence":
        """Static method to create a WorkScopedPersistence with default store.

        Args:
            work_id: The work identifier

        Returns:
            New WorkScopedPersistence instance
        """
        return WorkScopedPersistence(work_id)

    @property
    def work_id(self) -> str:
        """Get the work ID."""
        return self._work_id

    @property
    def store(self) -> WorkPersistence:
        """Get the underlying WorkPersistence store."""
        return self._store

    def __repr__(self) -> str:
        """Return string representation."""
        return f"WorkScopedPersistence(work_id={self._work_id})"

    # Allows accessing non-overridden methods/attributes directly on the store
    def __getattr__(self, name: str):
        """Delegate attribute access to the underlying store."""
        return getattr(self._store, name)

    # === WORK ===
    def update_work_status(self, status: str, **fields: Any) -> None:
        """Update the status of the current work.

        Args:
            status: New status value
            **fields: Additional fields to update
        """
        self._store.work_update(self._work_id, status=status, **fields)

    def get_work_status(self) -> str | None:
        """Get the status of the current work.

        Returns:
            Work status or None if work not found
        """
        work = self._store.work_get(self._work_id)
        return work["status"] if work else None

    def get_work_data(self) -> dict[str, Any] | None:
        """Get complete data for the current work.

        Returns:
            Work data dictionary or None if work not found
        """
        return self._store.work_get(self._work_id)

    def get_work_output_path(self, work_id: str | None = None):
        """Get output_path for a work, defaulting to current work.

        Args:
            work_id: Optional explicit work id. If None, uses self.work_id.
        Returns:
            Path or None if not set.
        """
        from pathlib import Path

        target_id = work_id or self._work_id
        data = self._store.work_get(target_id)
        if not data:
            return None
        out = data.get("output_path")
        return Path(out) if out else None

    def algorithm_error(self, algorithm_id: str, error: Exception) -> None:
        """Log algorithm error for the current work.

        Args:
            algorithm_id: Algorithm identifier
            error: Exception that occurred
        """
        # Use event logging for algorithm errors
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="algorithm",
            entity_data={
                "algorithm_id": algorithm_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    # === COMBINATIONS ===
    def submit_combinations(self, tasks_combinations: list[dict[str, Any]]) -> int:
        """Submit combinations for the current work using batch insertion.

        Args:
            tasks_combinations: List of combination dictionaries

        Returns:
            Number of combinations created
        """
        if not tasks_combinations:
            return 0

        # Add work_id to all combinations
        combinations_with_work_id = []
        for combo in tasks_combinations:
            combo_with_work_id = combo.copy()
            combo_with_work_id["work_id"] = self._work_id
            combinations_with_work_id.append(combo_with_work_id)

        # Use batch insertion for better performance
        return self._store.combination_bulk_create(combinations_with_work_id)

    def update_combination_status(self, *args):
        """Update the status of a combination in the current work.

        Supports two call styles for backward compatibility:
        - update_combination_status(combination_id: int, status: str)
        - update_combination_status(task_id, dataset_id, preset_id, algorithm_id, status)
        """
        import time

        # Parse arguments
        if len(args) == 2 and isinstance(args[0], int):
            combination_id, status = args
        elif len(args) == 5:
            task_id, dataset_id, preset_id, algorithm_id, status = args
            # Find combination id by composite keys
            combos, _ = self._store.combination_list(
                filters={
                    "work_id": self._work_id,
                    "task_id": task_id,
                    "dataset_id": dataset_id,
                    "preset_id": preset_id,
                    "algorithm_id": algorithm_id,
                },
                limit=1,
            )
            if not combos:
                return
            combination_id = combos[0]["id"]
        else:
            raise TypeError(
                "update_combination_status expects (id, status) or (task_id, dataset_id, preset_id, algorithm_id, status)"
            )

        fields = {"status": status}
        current_time = time.time()
        if status == "running":
            fields["started_at"] = current_time
        elif status in ["completed", "failed", "error", "canceled"]:
            fields["finished_at"] = current_time

        self._store.combination_update(combination_id, **fields)

    def get_next_queued_combination(self) -> Optional[dict[str, Any]]:
        """Get the next queued combination for the current work.

        Returns:
            Next queued combination or None if none available
        """
        combinations, _ = self._store.combination_list(
            filters={"work_id": self._work_id, "status": "queued"}, limit=1
        )
        return combinations[0] if combinations else None

    def get_next_pending_combination(self) -> Optional[dict[str, Any]]:
        """Compatibility: get next pending combination for the current work.

        Returns:
            Next pending combination or None if none available
        """
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
        """List combinations for the current work with optional filters.

        Args:
            status: Filter by status
            task_id: Filter by task ID
            dataset_id: Filter by dataset ID
            preset_id: Filter by preset ID
            algorithm_id: Filter by algorithm ID
            mode: Filter by mode
            order_by: Field to order by
            order_direction: Order direction (ASC/DESC)
            limit: Maximum number of results
            offset: Number of results to skip

        Returns:
            List of combination dictionaries
        """
        # Build filters dict
        filters = {"work_id": self._work_id}
        if status is not None:
            filters["status"] = status
        if task_id is not None:
            filters["task_id"] = task_id
        if dataset_id is not None:
            filters["dataset_id"] = dataset_id
        if preset_id is not None:
            filters["preset_id"] = preset_id
        if algorithm_id is not None:
            filters["algorithm_id"] = algorithm_id
        if mode is not None:
            filters["mode"] = mode

        order_desc = order_direction.upper() == "DESC"

        results, _ = self._store.combination_list(
            filters=filters,
            order_by=order_by,
            order_desc=order_desc,
            limit=limit,
            offset=offset or 0,
        )
        return results

    def init_combination(self) -> bool:
        """Reset running/paused/canceled combinations for the current work to 'queued'.

        Returns:
            True if any combinations were reset, False otherwise
        """
        # Get combinations to reset
        combinations, _ = self._store.combination_list(
            filters={
                "work_id": self._work_id,
                "status": ["running", "paused", "canceled"],
            }
        )

        # Reset their status to queued
        for combo in combinations:
            self._store.combination_update(combo["id"], status="queued")

        return len(combinations) > 0

    # === EVENTS ===
    # Work-level events
    def log_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        """Log warning for the current work.

        Args:
            message: Warning message
            context: Additional context data
        """
        self.work_warning(message, context)

    def log_error(self, error: Exception) -> None:
        """Log error for the current work.

        Args:
            error: Exception that occurred
        """
        self.work_error(error)

    # Work
    def work_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        """Log a work-level warning.

        Args:
            message: Warning message
            context: Additional context data
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="work",
            entity_data={"message": message, **(context or {})},
        )

    def work_error(self, error: Exception) -> None:
        """Log a work-level error.

        Args:
            error: Exception that occurred
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="work",
            entity_data={
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    # Task
    def task_warning(
        self, task_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Log a task-level warning.

        Args:
            task_id: Task identifier
            message: Warning message
            context: Additional context data
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="task",
            entity_data={"task_id": task_id, "message": message, **(context or {})},
        )

    def task_error(self, task_id: str, error: Exception) -> None:
        """Log a task-level error.

        Args:
            task_id: Task identifier
            error: Exception that occurred
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="task",
            entity_data={
                "task_id": task_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    # Dataset
    def dataset_warning(
        self, dataset_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Log a dataset-level warning.

        Args:
            dataset_id: Dataset identifier
            message: Warning message
            context: Additional context data
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="dataset",
            entity_data={
                "dataset_id": dataset_id,
                "message": message,
                **(context or {}),
            },
        )

    def dataset_error(self, dataset_id: str, error: Exception) -> None:
        """Log a dataset-level error.

        Args:
            dataset_id: Dataset identifier
            error: Exception that occurred
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="dataset",
            entity_data={
                "dataset_id": dataset_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    def submit_dataset(
        self, id: str, dataset_obj, meta: dict[str, Any] | None = None
    ) -> None:
        """Submit a dataset for the current work.

        Args:
            id: Dataset identifier
            dataset_obj: Dataset object containing sequences
            meta: Additional metadata
        """
        # Create dataset entry with auto-generated ID
        dataset_auto_id = self._store.dataset_create(
            dataset_id=id,
            work_id=self._work_id,
            name=getattr(dataset_obj, "name", None),
            meta=meta,
        )

        # Store sequences if the dataset has them
        if hasattr(dataset_obj, "sequences"):
            for seq_index, sequence in enumerate(dataset_obj.sequences):
                self._store.dataset_sequence_create(
                    dataset_id=dataset_auto_id, seq_index=seq_index, sequence=sequence
                )

    def get_dataset(self, dataset_id: str):
        """Get a dataset for the current work.

        Args:
            dataset_id: Dataset identifier

        Returns:
            Dataset object or None if not found
        """
        # Find dataset by dataset_id and work_id
        datasets, _ = self._store.dataset_list(
            filters={"dataset_id": dataset_id, "work_id": self._work_id}, limit=1
        )
        if not datasets:
            return None

        dataset = datasets[0]
        auto_id = dataset["id"]

        # Get sequences
        sequences, _ = self._store.dataset_sequence_list(
            filters={"dataset_id": auto_id}, order_by="seq_index"
        )

        # Create Dataset instance with the retrieved data
        sequence_list = [seq["sequence"] for seq in sequences]
        return Dataset(
            id=dataset_id,
            name=dataset.get("name", "unnamed_dataset"),
            sequences=sequence_list,
        )

    # Preset
    def preset_warning(
        self, preset_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Log a preset-level warning.

        Args:
            preset_id: Preset identifier
            message: Warning message
            context: Additional context data
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="preset",
            entity_data={"preset_id": preset_id, "message": message, **(context or {})},
        )

    def preset_error(self, preset_id: str, error: Exception) -> None:
        """Log a preset-level error.

        Args:
            preset_id: Preset identifier
            error: Exception that occurred
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="preset",
            entity_data={
                "preset_id": preset_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    # Combination
    def combination_warning(
        self, combination_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Log a combination-level warning.

        Args:
            combination_id: Combination identifier
            message: Warning message
            context: Additional context data
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="combination",
            entity_data={
                "combination_id": combination_id,
                "message": message,
                **(context or {}),
            },
        )

    def combination_error(self, combination_id: str, error: Exception) -> None:
        """Log a combination-level error.

        Args:
            combination_id: Combination identifier
            error: Exception that occurred
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="combination",
            entity_data={
                "combination_id": combination_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    # Unit
    def unit_warning(
        self, unit_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Log a unit-level warning.

        Args:
            unit_id: Unit identifier
            message: Warning message
            context: Additional context data
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="warning",
            event_category="unit",
            entity_data={"unit_id": unit_id, "message": message, **(context or {})},
        )

    def unit_error(self, unit_id: str, error: Exception) -> None:
        """Log a unit-level error.

        Args:
            unit_id: Unit identifier
            error: Exception that occurred
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type="error",
            event_category="unit",
            entity_data={
                "unit_id": unit_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    # Generic
    def generic_event(
        self,
        unit_id: str,
        event_type: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log a generic unit-level event.

        Args:
            unit_id: Unit identifier
            event_type: Type of event
            message: Event message
            context: Additional context data
        """
        self._store.event_create(
            work_id=self._work_id,
            event_type=event_type,
            event_category="unit",
            entity_data={"unit_id": unit_id, "message": message, **(context or {})},
        )

    # Event queries
    def get_events(
        self,
        *,
        event_category: str | None = None,
        event_type: str | None = None,
        unit_id: str | None = None,
        limit: int = 100,
    ) -> list[dict[str, Any]]:
        """Get events for the current work with optional filters.

        Args:
            event_category: Filter by event category
            event_type: Filter by event type
            unit_id: Filter by unit ID
            limit: Maximum number of events to return

        Returns:
            List of event dictionaries
        """
        filters = {"work_id": self._work_id}
        if event_category:
            filters["event_category"] = event_category
        if event_type:
            filters["event_type"] = event_type
        if unit_id:
            # For unit_id, we need to check the entity_data JSON field
            pass  # We'll need to implement custom filtering for JSON fields

        events, _ = self._store.event_list(
            filters=filters, order_by="timestamp", order_desc=True, limit=limit
        )
        return events

    def get_events_by_category(
        self, event_category: str, limit: int = 100
    ) -> list[dict[str, Any]]:
        """Get events for the current work filtered by category.

        Args:
            event_category: Event category to filter by
            limit: Maximum number of events to return

        Returns:
            List of event dictionaries
        """
        return self.get_events(event_category=event_category, limit=limit)

    def get_warnings(
        self, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        """Get warning events for the current work.

        Args:
            event_category: Optional category filter
            limit: Maximum number of warnings to return

        Returns:
            List of warning event dictionaries
        """
        return self.get_events(
            event_type="warning", event_category=event_category, limit=limit
        )

    def get_errors(
        self, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        """Get error events for the current work.

        Args:
            event_category: Optional category filter
            limit: Maximum number of errors to return

        Returns:
            List of error event dictionaries
        """
        return self.get_events(
            event_type="error", event_category=event_category, limit=limit
        )

    def get_event_summary_by_category(self) -> dict[str, dict[str, int]]:
        """Get event summary grouped by category and type.

        Returns:
            Dictionary with event counts by category and type
        """
        # Get all events for this work
        events, _ = self._store.event_list(
            filters={"work_id": self._work_id},
            limit=None,  # Get all events for summary
        )

        # Group by category and type
        summary = {}
        for event in events:
            category = event["event_category"]
            event_type = event["event_type"]

            if category not in summary:
                summary[category] = {}
            if event_type not in summary[category]:
                summary[category][event_type] = 0
            summary[category][event_type] += 1

        return summary

    # === FACTORY METHODS ===
    def for_combination(self, combination_id: int) -> "CombinationScopedPersistence":
        """Create a CombinationScopedPersistence for the specified combination within this work.

        Args:
            combination_id: ID of the combination

        Returns:
            CombinationScopedPersistence instance

        Raises:
            ValueError: If combination doesn't belong to this work
        """
        from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
            CombinationScopedPersistence,
        )

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
        """Create an ExecutionScopedPersistence for the specified execution unit within this work.

        Args:
            unit_id: ID of the execution unit

        Returns:
            ExecutionScopedPersistence instance

        Raises:
            ValueError: If execution doesn't belong to this work
        """
        from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
            ExecutionScopedPersistence,
        )

        return ExecutionScopedPersistence(unit_id, self._work_id, self._store)

    def get_all_combinations(self) -> list["CombinationScopedPersistence"]:
        """Get all combinations in this work as CombinationScopedPersistence instances.

        Returns:
            List of CombinationScopedPersistence instances
        """
        combinations = self.get_combinations()
        return [self.for_combination(combo["id"]) for combo in combinations]

    def get_all_executions(self) -> list["ExecutionScopedPersistence"]:
        """Get all executions in this work as ExecutionScopedPersistence instances.

        Returns:
            List of ExecutionScopedPersistence instances
        """
        from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
            ExecutionScopedPersistence,
        )

        # Get all executions across all combinations in this work
        all_executions = []
        for combination in self.get_combinations():
            combination_id = combination["id"]
            combo_scoped = self.for_combination(combination_id)
            executions = combo_scoped.get_executions()
            for execution in executions:
                unit_id = execution["unit_id"]
                all_executions.append(
                    ExecutionScopedPersistence(unit_id, self._work_id, self._store)
                )

        return all_executions

    def get_running_executions(self) -> list[dict[str, Any]]:
        """Get all executions that are currently in 'running' status.

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

    # === Progress Cleanup ===
    def clear_progress_for_non_finalized_executions(self) -> int:
        """Clear progress entries for all non-finalized executions in this work.

        Returns:
            Number of progress entries removed
        """
        return self._store.execution_progress_clear_for_non_finalized(
            work_id=self._work_id
        )
