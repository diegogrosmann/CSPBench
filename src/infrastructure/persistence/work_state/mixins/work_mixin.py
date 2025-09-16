"""Work CRUD operations mixin.

This module provides CRUD operations for the Work table, including creation,
retrieval, updates with status transition validation, deletion, and advanced
listing with filtering and pagination capabilities.
"""

import time
from typing import Any, Dict, List, Optional, Tuple


class WorkCRUDMixin:
    """Mixin providing CRUD operations for Work table.

    This mixin implements all basic database operations for Work entities,
    including status transition validation according to business rules.
    It should be used as a base class alongside other mixins to build
    a complete persistence layer.
    """

    def work_create(
        self,
        *,
        id: str,
        config: Dict[str, Any] | None = None,
        status: str = "queued",
        output_path: Optional[str] = None,
        error: Optional[str] = None,
        extra: Dict[str, Any] | None = None,
        created_at: Optional[float] = None,
        updated_at: Optional[float] = None,
    ) -> None:
        """Create a new work entry in the database.

        Args:
            id: Unique identifier for the work entry
            config: Configuration data for the work (stored as JSON)
            status: Initial status of the work, defaults to "queued"
            output_path: Path to output files, if any
            error: Error message, if applicable
            extra: Additional metadata (stored as JSON)
            created_at: Creation timestamp, defaults to current time
            updated_at: Last update timestamp, defaults to current time

        Note:
            If created_at or updated_at are not provided, they will be set
            to the current timestamp.
        """
        from ..models import Work

        ts = time.time()
        created_at = created_at if created_at is not None else ts
        updated_at = updated_at if updated_at is not None else ts

        with self.session_scope() as session:
            work = Work(
                id=id,
                config_json=config or {},
                status=status,
                created_at=created_at,
                updated_at=updated_at,
                output_path=output_path,
                error=error,
                extra_json=extra or {},
            )
            session.add(work)

    def work_get(self, id: str) -> Optional[Dict[str, Any]]:
        """Retrieve a work entry by its ID.

        Args:
            id: The unique identifier of the work entry

        Returns:
            Dictionary containing work data if found, None otherwise.
            The dictionary includes all work fields converted to appropriate types.
        """
        from ..models import Work

        with self.session_scope() as session:
            work = session.query(Work).filter(Work.id == id).first()
            return work.to_dict() if work else None

    def work_update(self, id: str, **fields: Any) -> bool:
        """Update a work entry with status transition validation.

        This method updates work fields while enforcing status transition rules
        defined in the domain layer. Invalid status transitions will be rejected.

        Args:
            id: The unique identifier of the work entry to update
            **fields: Key-value pairs of fields to update

        Returns:
            True if the update was successful, False if the work was not found
            or if an invalid status transition was attempted.

        Note:
            The updated_at field is automatically set to the current timestamp
            unless explicitly provided in the fields.
        """
        from src.domain.status import ALLOWEDSTATUS, BaseStatus

        from ..models import Work

        if not fields:
            return True

        # Auto-update updated_at if not provided
        fields.setdefault("updated_at", time.time())

        with self.session_scope() as session:
            work = session.query(Work).filter(Work.id == id).first()
            if not work:
                return False

            # Validate status transition if status is being updated
            if "status" in fields:
                new_status = fields["status"]
                current_status = work.status

                try:
                    current_status_enum = BaseStatus(current_status)
                    new_status_enum = BaseStatus(new_status)

                    # Check if transition is allowed using only the defined rules
                    allowed_transitions = ALLOWEDSTATUS.get(current_status_enum, set())
                    if new_status_enum not in allowed_transitions:
                        return False  # Invalid transition

                except ValueError:
                    # Invalid status value
                    return False

            # Apply updates
            for key, value in fields.items():
                if key.endswith("_json") or key in {"config_json", "extra_json"}:
                    setattr(work, key, value)
                else:
                    setattr(work, key, value)
            return True

    def work_delete(self, id: str) -> None:
        """Delete a work entry from the database.

        Args:
            id: The unique identifier of the work entry to delete

        Note:
            This operation is idempotent - no error is raised if the
            work entry doesn't exist.
        """
        from ..models import Work

        with self.session_scope() as session:
            work = session.query(Work).filter(Work.id == id).first()
            if work:
                session.delete(work)

    def work_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False,
    ) -> Tuple[List[Dict[str, Any]], int]:
        """List work entries with advanced filtering and pagination.

        Supports various filtering operators and pagination for efficient
        data retrieval. Filters can use simple equality, list membership,
        or advanced operators like LIKE, comparison operators, etc.

        Args:
            filters: Dictionary of field:value pairs for filtering. Supports:
                - Simple values: {'status': 'running'}
                - List values: {'status': ['running', 'queued']}
                - Advanced operators: {'created_at': {'operator': 'gte', 'value': 123456}}
            offset: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            order_by: Field name to order results by
            order_desc: Whether to order in descending order

        Returns:
            Tuple containing:
                - List of work dictionaries matching the criteria
                - Total count of matching records (before pagination)

        Note:
            Advanced operators supported: 'like', 'ilike', 'gt', 'gte',
            'lt', 'lte', 'ne' for flexible querying.
        """
        from ..models import Work

        with self.session_scope() as session:
            query = session.query(Work)

            # Apply filters
            if filters:
                for field, value in filters.items():
                    if hasattr(Work, field):
                        if isinstance(value, list):
                            query = query.filter(getattr(Work, field).in_(value))
                        elif isinstance(value, dict) and "operator" in value:
                            # Support for advanced operators
                            op = value["operator"]
                            val = value["value"]
                            column = getattr(Work, field)

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
                            query = query.filter(getattr(Work, field) == value)

            # Get total count before pagination
            total_count = query.count()

            # Apply ordering
            if order_by and hasattr(Work, order_by):
                column = getattr(Work, order_by)
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
            results = [work.to_dict() for work in query.all()]

            return results, total_count
