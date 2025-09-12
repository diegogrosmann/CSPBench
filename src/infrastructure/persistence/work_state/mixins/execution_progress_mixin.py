"""Execution Progress CRUD operations mixin.

This module provides CRUD operations for tracking execution progress,
including creation, retrieval, updates, deletion, and bulk operations
for managing progress data efficiently.
"""

import time
from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import and_, or_


class ExecutionProgressCRUDMixin:
    """Mixin providing CRUD operations for ExecutionProgress table.
    
    This mixin handles progress tracking data for executions, allowing
    detailed monitoring of execution progress over time. It supports
    both individual operations and bulk operations for performance.
    """

    def execution_progress_create(
        self,
        *,
        execution_id: int,
        progress: float,
        message: Optional[str] = None,
        timestamp: Optional[float] = None,
    ) -> None:
        """Create a new execution progress entry.
        
        Args:
            execution_id: ID of the execution this progress belongs to
            progress: Progress value (typically 0.0 to 1.0)
            message: Optional progress message or description
            timestamp: Progress timestamp, defaults to current time
        """
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
        """Retrieve an execution progress entry by ID.
        
        Args:
            id: Primary key of the progress entry
            
        Returns:
            Dictionary containing progress data if found, None otherwise.
        """
        from ..models import ExecutionProgress
        
        with self.session_scope() as session:
            progress = session.query(ExecutionProgress).filter(ExecutionProgress.id == id).first()
            return progress.to_dict() if progress else None

    def execution_progress_update(self, id: int, **fields: Any) -> None:
        """Update an execution progress entry.
        
        Args:
            id: Primary key of the progress entry to update
            **fields: Key-value pairs of fields to update
            
        Note:
            This operation is idempotent - no error if entry doesn't exist.
        """
        from ..models import ExecutionProgress
        
        if not fields:
            return
        
        with self.session_scope() as session:
            progress = session.query(ExecutionProgress).filter(ExecutionProgress.id == id).first()
            if progress:
                for key, value in fields.items():
                    setattr(progress, key, value)

    def execution_progress_delete(self, id: int) -> None:
        """Delete an execution progress entry.
        
        Args:
            id: Primary key of the progress entry to delete
            
        Note:
            This operation is idempotent - no error if entry doesn't exist.
        """
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
        """List execution progress entries with filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering. Supports:
                - Simple values: {'execution_id': 123}
                - List values: {'execution_id': [123, 456]}
                - Advanced operators: {'progress': {'operator': 'gte', 'value': 0.5}}
            offset: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            order_by: Field name to order results by
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple containing:
                - List of progress entry dictionaries
                - Total count of matching records (before pagination)
                
        Note:
            Advanced operators supported: 'like', 'ilike', 'gt', 'gte', 
            'lt', 'lte', 'ne' for flexible querying.
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

    def execution_progress_clear_for_non_finalized(
        self,
        *,
        work_id: Optional[str] = None,
        combination_id: Optional[int] = None,
        execution_id: Optional[int] = None,
    ) -> int:
        """Clear progress entries for executions that are not finalized.
        
        Removes progress tracking data for executions in non-finalized states
        (queued, running, paused). This is useful for cleanup operations when
        restarting or resetting work execution.
        
        Args:
            work_id: If provided, only clear progress for executions in this work
            combination_id: If provided, only clear progress for executions in this combination
            execution_id: If provided, only clear progress for this specific execution
            
        Returns:
            Number of progress entries that were deleted.
            
        Note:
            Non-finalized statuses are: 'queued', 'running', 'paused'
            Finalized statuses are: 'completed', 'failed', 'error', 'canceled'
        """
        from ..models import ExecutionProgress, Execution, Combination
        
        with self.session_scope() as session:
            # Build base query to find progress entries for non-finalized executions
            query = session.query(ExecutionProgress).join(
                Execution, ExecutionProgress.execution_id == Execution.id
            )
            
            # Filter by non-finalized execution statuses
            non_finalized_statuses = ['queued', 'running', 'paused']
            query = query.filter(Execution.status.in_(non_finalized_statuses))
            
            # Apply additional filters if provided
            if work_id is not None:
                # Join with combination to filter by work_id
                query = query.join(Combination, Execution.combination_id == Combination.id)
                query = query.filter(Combination.work_id == work_id)
            
            if combination_id is not None:
                query = query.filter(Execution.combination_id == combination_id)
                
            if execution_id is not None:
                query = query.filter(Execution.id == execution_id)
            
            # Get the progress entries to delete
            progress_entries = query.all()
            deleted_count = len(progress_entries)
            
            # Delete the progress entries
            for progress in progress_entries:
                session.delete(progress)
            
            return deleted_count