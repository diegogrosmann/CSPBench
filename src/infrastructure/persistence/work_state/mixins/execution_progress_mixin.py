"""Execution Progress CRUD mixin."""

import time
from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import and_, or_


class ExecutionProgressCRUDMixin:
    """CRUD operations for ExecutionProgress table."""

    def execution_progress_create(
        self,
        *,
        execution_id: int,
        progress: float,
        message: Optional[str] = None,
        timestamp: Optional[float] = None,
    ) -> None:
        """Create execution progress entry."""
        from ..models import ExecutionProgress
        
        with self.session_scope() as session:
            progress_entry = ExecutionProgress(
                execution_id=execution_id,
                progress=progress,
                message=message,
                timestamp=timestamp if timestamp is not None else time.time(),
            )
            session.add(progress_entry)

    def execution_progress_get(self, id: int) -> Optional[Dict[str, Any]]:
        """Get execution progress by ID."""
        from ..models import ExecutionProgress
        
        with self.session_scope() as session:
            progress = session.query(ExecutionProgress).filter(ExecutionProgress.id == id).first()
            return progress.to_dict() if progress else None

    def execution_progress_update(self, id: int, **fields: Any) -> None:
        """Update execution progress entry."""
        from ..models import ExecutionProgress
        
        if not fields:
            return
        
        with self.session_scope() as session:
            progress = session.query(ExecutionProgress).filter(ExecutionProgress.id == id).first()
            if progress:
                for key, value in fields.items():
                    setattr(progress, key, value)

    def execution_progress_delete(self, id: int) -> None:
        """Delete execution progress entry."""
        from ..models import ExecutionProgress
        
        with self.session_scope() as session:
            progress = session.query(ExecutionProgress).filter(ExecutionProgress.id == id).first()
            if progress:
                session.delete(progress)

    def execution_progress_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List execution progress entries with generic filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering
            offset: Number of records to skip
            limit: Maximum number of records to return
            order_by: Field to order by
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple of (records, total_count)
        """
        from ..models import ExecutionProgress
        
        with self.session_scope() as session:
            query = session.query(ExecutionProgress)
            
            # Apply filters
            if filters:
                for field, value in filters.items():
                    if hasattr(ExecutionProgress, field):
                        if isinstance(value, list):
                            query = query.filter(getattr(ExecutionProgress, field).in_(value))
                        elif isinstance(value, dict) and 'operator' in value:
                            # Support for advanced operators
                            op = value['operator']
                            val = value['value']
                            column = getattr(ExecutionProgress, field)
                            
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
                            query = query.filter(getattr(ExecutionProgress, field) == value)
            
            # Get total count before pagination
            total_count = query.count()
            
            # Apply ordering
            if order_by and hasattr(ExecutionProgress, order_by):
                column = getattr(ExecutionProgress, order_by)
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
            results = [progress.to_dict() for progress in query.all()]
            
            return results, total_count