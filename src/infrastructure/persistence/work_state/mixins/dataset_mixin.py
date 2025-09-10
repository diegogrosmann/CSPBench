"""Dataset CRUD mixin."""

from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import and_, or_


class DatasetCRUDMixin:
    """CRUD operations for Dataset table."""

    def dataset_create(
        self,
        *,
        dataset_id: str,
        work_id: Optional[str],
        name: Optional[str] = None,
        meta: Dict[str, Any] | None = None,
    ) -> int:
        """Create a new dataset entry. Returns the auto-generated ID."""
        from ..models import Dataset
        
        with self.session_scope() as session:
            dataset = Dataset(
                dataset_id=dataset_id,
                work_id=work_id,
                name=name,
                meta_json=meta or {},
            )
            session.add(dataset)
            session.flush()  # Force ID generation
            return dataset.id

    def dataset_get(self, id: int) -> Optional[Dict[str, Any]]:
        """Get dataset by auto-generated ID."""
        from ..models import Dataset
        
        with self.session_scope() as session:
            dataset = session.query(Dataset).filter(Dataset.id == id).first()
            return dataset.to_dict() if dataset else None

    def dataset_update(self, id: int, **fields: Any) -> None:
        """Update dataset entry by auto-generated ID."""
        from ..models import Dataset
        
        if not fields:
            return
        
        with self.session_scope() as session:
            dataset = session.query(Dataset).filter(Dataset.id == id).first()
            if dataset:
                for field, value in fields.items():
                    if field == "meta":
                        dataset.meta_json = value
                    elif hasattr(dataset, field):
                        setattr(dataset, field, value)

    def dataset_delete(self, id: int) -> None:
        """Delete dataset entry by auto-generated ID."""
        from ..models import Dataset
        
        with self.session_scope() as session:
            dataset = session.query(Dataset).filter(Dataset.id == id).first()
            if dataset:
                session.delete(dataset)

    def dataset_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List datasets with generic filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering
            offset: Number of records to skip
            limit: Maximum number of records to return
            order_by: Field to order by
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple of (records, total_count)
        """
        from ..models import Dataset
        
        with self.session_scope() as session:
            query = session.query(Dataset)
            
            # Apply filters
            if filters:
                for field, value in filters.items():
                    if hasattr(Dataset, field):
                        if isinstance(value, list):
                            query = query.filter(getattr(Dataset, field).in_(value))
                        elif isinstance(value, dict) and 'operator' in value:
                            # Support for advanced operators
                            op = value['operator']
                            val = value['value']
                            column = getattr(Dataset, field)
                            
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
                            query = query.filter(getattr(Dataset, field) == value)
            
            # Get total count before pagination
            total_count = query.count()
            
            # Apply ordering
            if order_by and hasattr(Dataset, order_by):
                column = getattr(Dataset, order_by)
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
            results = [dataset.to_dict() for dataset in query.all()]
            
            return results, total_count