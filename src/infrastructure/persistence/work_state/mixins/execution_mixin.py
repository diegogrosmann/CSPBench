"""Execution CRUD operations mixin.

This module provides CRUD operations for the Execution table, including
individual operations and bulk operations for efficient execution management.
"""

from typing import Any, Dict, List, Optional, Tuple

from sqlalchemy import and_


class ExecutionCRUDMixin:
    """Mixin providing CRUD operations for Execution table.

    This mixin handles individual execution records within combinations,
    supporting both single operations and bulk operations for performance
    when dealing with large numbers of executions.
    """

    def execution_create(
        self,
        *,
        unit_id: str,
        combination_id: int,
        sequencia: int,
        status: str = "queued",
        started_at: Optional[float] = None,
        finished_at: Optional[float] = None,
        params: Dict[str, Any] | None = None,
        result: Dict[str, Any] | None = None,
        objective: Optional[float] = None,
    ) -> None:
        """Create a new execution entry.

        Args:
            unit_id: Unique identifier for the execution unit
            combination_id: ID of the combination this execution belongs to
            sequencia: Sequence number within the combination
            status: Initial execution status, defaults to "queued"
            started_at: Timestamp when execution started
            finished_at: Timestamp when execution finished
            params: Execution parameters (stored as JSON)
            result: Execution results (stored as JSON)
            objective: Objective function value, if applicable
        """
        from ..models import Execution

        with self.session_scope() as session:
            execution = Execution(
                unit_id=unit_id,
                combination_id=combination_id,
                sequencia=sequencia,
                status=status,
                started_at=started_at,
                finished_at=finished_at,
                params_json=params or {},
                result_json=result or {},
                objective=objective,
            )
            session.add(execution)

    def execution_get(self, id: int) -> Optional[Dict[str, Any]]:
        """Retrieve an execution by its ID.

        Args:
            id: Primary key of the execution

        Returns:
            Dictionary containing execution data if found, None otherwise.
        """
        from ..models import Execution

        with self.session_scope() as session:
            execution = session.query(Execution).filter(Execution.id == id).first()
            return execution.to_dict() if execution else None

    def execution_update(self, id: int, **fields: Any) -> None:
        """Update an execution entry.

        Args:
            id: Primary key of the execution to update
            **fields: Key-value pairs of fields to update. Special handling
                     for 'params' and 'result' keys which map to JSON fields.

        Note:
            - 'params' maps to params_json field
            - 'result' maps to result_json field
            - Other fields are set directly on the model
        """
        from ..models import Execution

        if not fields:
            return

        with self.session_scope() as session:
            execution = session.query(Execution).filter(Execution.id == id).first()
            if execution:
                for key, value in fields.items():
                    # Handle JSON fields specifically
                    if key == "params":
                        execution.params_json = value
                    elif key == "result":
                        execution.result_json = value
                    else:
                        setattr(execution, key, value)

    def execution_delete(self, id: int) -> None:
        """Delete an execution entry.

        Args:
            id: Primary key of the execution to delete

        Note:
            This operation is idempotent - no error if execution doesn't exist.
        """
        from ..models import Execution

        with self.session_scope() as session:
            execution = session.query(Execution).filter(Execution.id == id).first()
            if execution:
                session.delete(execution)

    def execution_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False,
    ) -> Tuple[List[Dict[str, Any]], int]:
        """List executions with advanced filtering and pagination.

        Args:
            filters: Dictionary of field:value pairs for filtering. Supports:
                - Simple values: {'status': 'running'}
                - List values: {'status': ['running', 'queued']}
                - Advanced operators: {'started_at': {'operator': 'gte', 'value': 123456}}
            offset: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            order_by: Field name to order results by
            order_desc: Whether to order in descending order

        Returns:
            Tuple containing:
                - List of execution dictionaries matching the criteria
                - Total count of matching records (before pagination)

        Note:
            Advanced operators supported: 'like', 'ilike', 'gt', 'gte',
            'lt', 'lte', 'ne' for flexible querying.
        """
        from ..models import Execution

        with self.session_scope() as session:
            query = session.query(Execution)

            # Apply filters
            if filters:
                for field, value in filters.items():
                    if hasattr(Execution, field):
                        if isinstance(value, list):
                            query = query.filter(getattr(Execution, field).in_(value))
                        elif isinstance(value, dict) and "operator" in value:
                            # Support for advanced operators
                            op = value["operator"]
                            val = value["value"]
                            column = getattr(Execution, field)

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
                            query = query.filter(getattr(Execution, field) == value)

            # Get total count before pagination
            total_count = query.count()

            # Apply ordering
            if order_by and hasattr(Execution, order_by):
                column = getattr(Execution, order_by)
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
            results = [execution.to_dict() for execution in query.all()]

            return results, total_count

    def execution_bulk_create(self, executions: List[Dict[str, Any]]) -> int:
        """Create multiple executions efficiently in a single transaction.

        Uses bulk insert operations for better performance when creating
        large numbers of executions. Falls back to individual inserts
        if bulk operation fails due to constraint violations.

        Args:
            executions: List of execution dictionaries with required fields.
                       Each dictionary should contain at least unit_id,
                       combination_id, and sequencia.

        Returns:
            Number of executions successfully inserted. May be less than
            the input list length if some insertions fail due to constraints.

        Note:
            - Automatically sets default status to 'queued' if not provided
            - Automatically initializes empty JSON fields if not provided
            - Ignores duplicate key constraint violations gracefully
        """
        from ..models import Execution

        if not executions:
            return 0

        # Set default values for all executions if not provided
        for exec_data in executions:
            if "status" not in exec_data:
                exec_data["status"] = "queued"
            if "params_json" not in exec_data:
                exec_data["params_json"] = {}
            if "result_json" not in exec_data:
                exec_data["result_json"] = {}

        with self.session_scope() as session:
            try:
                # Use SQLAlchemy's bulk_insert_mappings for better performance
                session.bulk_insert_mappings(Execution, executions)
                return len(executions)
            except Exception:
                # If bulk insert fails (e.g., duplicate key constraint),
                # fallback to individual inserts with ignore duplicates
                session.rollback()
                inserted_count = 0

                for exec_data in executions:
                    try:
                        # Try to insert individually
                        execution = Execution(**exec_data)
                        session.add(execution)
                        session.flush()  # Force execution to catch constraint violations
                        inserted_count += 1
                    except Exception:
                        # Ignore duplicates and other constraint violations
                        session.rollback()
                        continue

                return inserted_count

    def execution_bulk_update(
        self, filters: Dict[str, Any], update_fields: Dict[str, Any]
    ) -> int:
        """Bulk update executions matching the specified filters.

        Efficiently updates multiple executions in a single database operation
        based on filter criteria.

        Args:
            filters: Dictionary with filter conditions. Supports:
                - Simple values: {'combination_id': 123, 'status': 'running'}
                - List values: {'combination_id': [1,2,3], 'status': ['running', 'paused']}
            update_fields: Dictionary with fields to update:
                - 'result' maps to result_json field
                - 'params' maps to params_json field
                - Other fields are updated directly

        Returns:
            Number of executions that were actually updated.

        Example:
            # Reset all running executions in combinations 1,2,3 to queued
            count = execution_bulk_update(
                filters={'combination_id': [1,2,3], 'status': 'running'},
                update_fields={'status': 'queued', 'started_at': None}
            )
        """
        from ..models import Execution

        if not filters or not update_fields:
            return 0

        # Handle JSON fields properly
        processed_fields = {}
        for key, value in update_fields.items():
            if key == "result":
                processed_fields["result_json"] = value
            elif key == "params":
                processed_fields["params_json"] = value
            else:
                processed_fields[key] = value

        with self.session_scope() as session:
            query = session.query(Execution)

            # Apply filters
            filter_conditions = []
            for key, value in filters.items():
                if isinstance(value, list):
                    filter_conditions.append(getattr(Execution, key).in_(value))
                else:
                    filter_conditions.append(getattr(Execution, key) == value)

            if filter_conditions:
                query = query.filter(and_(*filter_conditions))

            # Perform bulk update
            result = query.update(processed_fields, synchronize_session=False)
            return result
