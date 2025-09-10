"""Combination CRUD mixin."""

import time
from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import and_, or_


class CombinationCRUDMixin:
    """CRUD operations for Combination table."""

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
        """Create combination entry."""
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
        """Get combination by ID."""
        from ..models import Combination
        
        with self.session_scope() as session:
            combination = session.query(Combination).filter(Combination.id == id).first()
            return combination.to_dict() if combination else None

    def combination_update(self, id: int, **fields: Any) -> None:
        """Update combination entry."""
        from ..models import Combination
        
        if not fields:
            return
        
        with self.session_scope() as session:
            combination = session.query(Combination).filter(Combination.id == id).first()
            if combination:
                for key, value in fields.items():
                    setattr(combination, key, value)

    def combination_delete(self, id: int) -> None:
        """Delete combination entry."""
        from ..models import Combination
        
        with self.session_scope() as session:
            combination = session.query(Combination).filter(Combination.id == id).first()
            if combination:
                session.delete(combination)

    def combination_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List combinations with generic filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering
            offset: Number of records to skip
            limit: Maximum number of records to return
            order_by: Field to order by
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple of (records, total_count)
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
                        elif isinstance(value, dict) and 'operator' in value:
                            # Support for advanced operators
                            op = value['operator']
                            val = value['value']
                            column = getattr(Combination, field)
                            
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