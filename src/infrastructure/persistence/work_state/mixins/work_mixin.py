"""Work CRUD mixin."""

import time
from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import and_, or_


class WorkCRUDMixin:
    """CRUD operations for Work table."""

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
        """Create a new work entry."""
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
        """Get work by ID."""
        from ..models import Work
        
        with self.session_scope() as session:
            work = session.query(Work).filter(Work.id == id).first()
            return work.to_dict() if work else None

    def work_update(self, id: str, **fields: Any) -> bool:
        """Update work entry."""
        from ..models import Work
        
        if not fields:
            return True
        
        # Auto-update updated_at if not provided
        fields.setdefault("updated_at", time.time())
        
        with self.session_scope() as session:
            work = session.query(Work).filter(Work.id == id).first()
            if work:
                for key, value in fields.items():
                    if key.endswith('_json') or key in {'config_json', 'extra_json'}:
                        setattr(work, key, value)
                    else:
                        setattr(work, key, value)
                return True
            return False

    def work_delete(self, id: str) -> None:
        """Delete work entry."""
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
        order_desc: bool = False
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List works with generic filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering
            offset: Number of records to skip
            limit: Maximum number of records to return
            order_by: Field to order by
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple of (records, total_count)
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
                        elif isinstance(value, dict) and 'operator' in value:
                            # Support for advanced operators
                            op = value['operator']
                            val = value['value']
                            column = getattr(Work, field)
                            
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