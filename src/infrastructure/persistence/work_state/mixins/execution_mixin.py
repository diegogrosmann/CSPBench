"""Execution CRUD mixin."""

import time
from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import and_, or_


class ExecutionCRUDMixin:
    """CRUD operations for Execution table."""

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
        """Create execution entry."""
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
        """Get execution by ID."""
        from ..models import Execution
        
        with self.session_scope() as session:
            execution = session.query(Execution).filter(Execution.id == id).first()
            return execution.to_dict() if execution else None

    def execution_update(self, id: int, **fields: Any) -> None:
        """Update execution entry."""
        from ..models import Execution
        
        if not fields:
            return
        
        with self.session_scope() as session:
            execution = session.query(Execution).filter(Execution.id == id).first()
            if execution:
                for key, value in fields.items():
                    # Handle JSON fields specifically
                    if key == 'params':
                        execution.params_json = value
                    elif key == 'result':
                        execution.result_json = value
                    else:
                        setattr(execution, key, value)

    def execution_delete(self, id: int) -> None:
        """Delete execution entry."""
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
        order_desc: bool = False
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List executions with generic filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering
            offset: Number of records to skip
            limit: Maximum number of records to return
            order_by: Field to order by
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple of (records, total_count)
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
                        elif isinstance(value, dict) and 'operator' in value:
                            # Support for advanced operators
                            op = value['operator']
                            val = value['value']
                            column = getattr(Execution, field)
                            
                            if op == 'like':
                                query = query.filter(column.like(f"%{val}%"))
                            elif op == 'ilike':
                                query = query.filter(column.ilike(f"%{val}%"))
                            elif op == 'gt':
                                query = query.filter(column > val)
                            elif op == 'gte':
                                query = query.filter(column >= val)
                            elif op == 'lt':
                                query = query.filter(column < val)
                            elif op == 'lte':
                                query = query.filter(column <= val)
                            elif op == 'ne':
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