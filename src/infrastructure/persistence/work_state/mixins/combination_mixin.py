"""Combination CRUD operations mixin.

This module provides CRUD operations for the Combination table, which represents
unique combinations of task, dataset, preset, and algorithm configurations
within a work execution context.
"""

import time
from typing import Any, Dict, List, Optional, Tuple

from sqlalchemy import and_


class CombinationCRUDMixin:
    """Mixin providing CRUD operations for Combination table.

    This mixin handles combination entities which represent unique configurations
    of task, dataset, preset, and algorithm combinations within a work context.
    Each combination can contain multiple executions (sequences) and tracks
    its own execution state and progress.
    """

    def combination_create(
        self,
        *,
        work_id: str,
        task_id: Optional[str],
        dataset_id: Optional[str],
        preset_id: Optional[str],
        algorithm_id: Optional[str],
        mode: Optional[str],
        status: str = "queued",
        total_sequences: Optional[int] = None,
        created_at: Optional[float] = None,
        started_at: Optional[float] = None,
        finished_at: Optional[float] = None,
    ) -> None:
        """Create a new combination entry.

        Args:
            work_id: ID of the work this combination belongs to
            task_id: Task identifier for this combination
            dataset_id: Dataset identifier for this combination
            preset_id: Configuration preset identifier
            algorithm_id: Algorithm identifier
            mode: Execution mode (e.g., 'optimization', 'evaluation')
            status: Initial status, defaults to "queued"
            total_sequences: Total number of sequences/executions planned
            created_at: Creation timestamp, defaults to current time
            started_at: Start timestamp, if applicable
            finished_at: Completion timestamp, if applicable
        """
        from ..models import Combination

        # Set created_at timestamp if not provided
        created_at = created_at if created_at is not None else time.time()

        with self.session_scope() as session:
            combination = Combination(
                work_id=work_id,
                task_id=task_id,
                dataset_id=dataset_id,
                preset_id=preset_id,
                algorithm_id=algorithm_id,
                mode=mode,
                status=status,
                total_sequences=total_sequences,
                created_at=created_at,
                started_at=started_at,
                finished_at=finished_at,
            )
            session.add(combination)

    def combination_get(self, id: int) -> Optional[Dict[str, Any]]:
        """Retrieve a combination by its ID.

        Args:
            id: Primary key of the combination

        Returns:
            Dictionary containing combination data if found, None otherwise.
        """
        from ..models import Combination

        with self.session_scope() as session:
            combination = (
                session.query(Combination).filter(Combination.id == id).first()
            )
            return combination.to_dict() if combination else None

    def combination_update(self, id: int, **fields: Any) -> None:
        """Update a combination entry.

        Args:
            id: Primary key of the combination to update
            **fields: Key-value pairs of fields to update

        Note:
            This operation is idempotent - no error if combination doesn't exist.
        """
        from ..models import Combination

        if not fields:
            return

        with self.session_scope() as session:
            combination = (
                session.query(Combination).filter(Combination.id == id).first()
            )
            if combination:
                for key, value in fields.items():
                    setattr(combination, key, value)

    def combination_delete(self, id: int) -> None:
        """Delete a combination entry.

        Args:
            id: Primary key of the combination to delete

        Warning:
            This will also delete all associated executions due to foreign
            key constraints. This operation is idempotent.
        """
        from ..models import Combination

        with self.session_scope() as session:
            combination = (
                session.query(Combination).filter(Combination.id == id).first()
            )
            if combination:
                session.delete(combination)

    def combination_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False,
    ) -> Tuple[List[Dict[str, Any]], int]:
        """List combinations with filtering and pagination.

        Args:
            filters: Dictionary of field:value pairs for filtering. Supports:
                - Simple values: {'work_id': 'abc', 'status': 'running'}
                - List values: {'status': ['running', 'queued']}
                - Advanced operators: {'created_at': {'operator': 'gte', 'value': 123456}}
            offset: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            order_by: Field name to order results by
            order_desc: Whether to order in descending order

        Returns:
            Tuple containing:
                - List of combination dictionaries matching the criteria
                - Total count of matching records (before pagination)

        Note:
            Advanced operators supported: 'like', 'ilike', 'gt', 'gte',
            'lt', 'lte', 'ne' for flexible querying.
        """
        from ..models import Combination

        with self.session_scope() as session:
            query = session.query(Combination)

            # Apply filters
            if filters:
                for field, value in filters.items():
                    if hasattr(Combination, field):
                        if isinstance(value, list):
                            query = query.filter(getattr(Combination, field).in_(value))
                        elif isinstance(value, dict) and "operator" in value:
                            # Support for advanced operators
                            op = value["operator"]
                            val = value["value"]
                            column = getattr(Combination, field)

                            if op == "like":
                                query = query.filter(column.like(f"%{val}%"))
                            elif op == "ilike":
                                query = query.filter(column.ilike(f"%{val}%"))
                            elif op == "gt":
                                query = query.filter(column > val)
                            elif op == "gte":
                                query = query.filter(column >= val)
                            elif op == "lt":
                                query = query.filter(column < val)
                            elif op == "lte":
                                query = query.filter(column <= val)
                            elif op == "ne":
                                query = query.filter(column != val)
                        else:
                            query = query.filter(getattr(Combination, field) == value)

            # Get total count before pagination
            total_count = query.count()

            # Apply ordering
            if order_by and hasattr(Combination, order_by):
                column = getattr(Combination, order_by)
                if order_desc:
                    query = query.order_by(column.desc())
                else:
                    query = query.order_by(column)

            # Apply pagination
            if offset > 0:
                query = query.offset(offset)
            if limit is not None:
                query = query.limit(limit)

            # Execute query and convert to dicts
            results = [combination.to_dict() for combination in query.all()]

            return results, total_count

    def combination_bulk_create(self, combinations: List[Dict[str, Any]]) -> int:
        """Create multiple combinations efficiently in a single transaction.

        Uses bulk insert operations for better performance when creating
        large numbers of combinations. Falls back to individual inserts
        if bulk operation fails due to constraint violations.

        Args:
            combinations: List of combination dictionaries with required fields.
                         Each should contain at least work_id and identifiers.

        Returns:
            Number of combinations successfully inserted. May be less than
            the input list length if some insertions fail due to constraints.

        Note:
            - Automatically sets created_at timestamp if not provided
            - Automatically sets default status to 'queued' if not provided
            - Ignores duplicate key constraint violations gracefully
        """
        from ..models import Combination

        if not combinations:
            return 0

        # Set created_at timestamp for all combinations if not provided
        current_time = time.time()
        for combo in combinations:
            if "created_at" not in combo or combo["created_at"] is None:
                combo["created_at"] = current_time
            # Set default status if not provided
            if "status" not in combo:
                combo["status"] = "queued"

        with self.session_scope() as session:
            try:
                # Use SQLAlchemy's bulk_insert_mappings for better performance
                session.bulk_insert_mappings(Combination, combinations)
                return len(combinations)
            except Exception:
                # If bulk insert fails (e.g., duplicate key constraint),
                # fallback to individual inserts with ignore duplicates
                session.rollback()
                inserted_count = 0

                for combo in combinations:
                    try:
                        # Try to insert individually
                        combination = Combination(**combo)
                        session.add(combination)
                        session.flush()  # Force execution to catch constraint violations
                        inserted_count += 1
                    except Exception:
                        # Ignore duplicates and other constraint violations
                        session.rollback()
                        continue

                return inserted_count

    def combination_bulk_update(
        self, filters: Dict[str, Any], update_fields: Dict[str, Any]
    ) -> int:
        """Bulk update combinations matching the specified filters.

        Efficiently updates multiple combinations in a single database operation
        based on filter criteria. Useful for status transitions, batch operations,
        and cleanup tasks.

        Args:
            filters: Dictionary with filter conditions. Supports:
                - Simple values: {'work_id': 'abc', 'status': 'running'}
                - List values: {'work_id': ['abc', 'def'], 'status': ['running', 'paused']}
            update_fields: Dictionary with fields to update:
                - Any combination field can be updated
                - Common use cases: status transitions, timestamp updates

        Returns:
            Number of combinations that were actually updated.

        Example:
            # Reset all running combinations in work 'abc' to queued
            count = combination_bulk_update(
                filters={'work_id': 'abc', 'status': 'running'},
                update_fields={'status': 'queued', 'started_at': None}
            )
        """
        from ..models import Combination

        if not filters or not update_fields:
            return 0

        with self.session_scope() as session:
            query = session.query(Combination)

            # Apply filters
            filter_conditions = []
            for key, value in filters.items():
                if isinstance(value, list):
                    filter_conditions.append(getattr(Combination, key).in_(value))
                else:
                    filter_conditions.append(getattr(Combination, key) == value)

            if filter_conditions:
                query = query.filter(and_(*filter_conditions))

            # Perform bulk update
            result = query.update(update_fields, synchronize_session=False)
            return result
