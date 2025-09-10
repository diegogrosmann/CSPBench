"""Event CRUD mixin."""

import time
from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import desc


class EventCRUDMixin:
    """CRUD operations for Event table."""

    def event_create(
        self,
        *,
        work_id: str,
        event_type: str,
        event_category: str,
        entity_data: Dict[str, Any] | None = None,
        timestamp: Optional[float] = None,
    ) -> None:
        """Create event entry."""
        from ..models import Event
        
        with self.session_scope() as session:
            event = Event(
                work_id=work_id,
                event_type=event_type,
                event_category=event_category,
                entity_data_json=entity_data or {},
                timestamp=timestamp if timestamp is not None else time.time(),
            )
            session.add(event)

    def event_get(self, id: int) -> Optional[Dict[str, Any]]:
        """Get event by ID."""
        from ..models import Event
        
        with self.session_scope() as session:
            event = session.query(Event).filter(Event.id == id).first()
            return event.to_dict() if event else None

    def event_update(self, id: int, **fields: Any) -> None:
        """Update event entry."""
        from ..models import Event
        
        if not fields:
            return
        
        with self.session_scope() as session:
            event = session.query(Event).filter(Event.id == id).first()
            if event:
                for key, value in fields.items():
                    if key == 'entity_data':
                        event.entity_data_json = value
                    else:
                        setattr(event, key, value)

    def event_delete(self, id: int) -> None:
        """Delete event entry."""
        from ..models import Event
        
        with self.session_scope() as session:
            event = session.query(Event).filter(Event.id == id).first()
            if event:
                session.delete(event)

    def event_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List events with generic filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering
            offset: Number of records to skip
            limit: Maximum number of records to return
            order_by: Field to order by
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple of (records, total_count)
        """
        from ..models import Event
        
        with self.session_scope() as session:
            query = session.query(Event)
            
            # Apply filters
            if filters:
                for field, value in filters.items():
                    if hasattr(Event, field):
                        if isinstance(value, list):
                            query = query.filter(getattr(Event, field).in_(value))
                        elif isinstance(value, dict) and 'operator' in value:
                            # Support for advanced operators
                            op = value['operator']
                            val = value['value']
                            column = getattr(Event, field)
                            
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
                            query = query.filter(getattr(Event, field) == value)
            
            # Get total count before pagination
            total_count = query.count()
            
            # Apply ordering
            if order_by and hasattr(Event, order_by):
                column = getattr(Event, order_by)
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
            results = [event.to_dict() for event in query.all()]
            
            return results, total_count