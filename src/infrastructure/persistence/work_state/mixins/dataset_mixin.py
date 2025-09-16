"""Dataset CRUD operations mixin.

This module provides CRUD operations for the Dataset table, which stores
dataset metadata and information for CSP benchmark problems.
"""

from typing import Any, Dict, List, Optional, Tuple


class DatasetCRUDMixin:
    """Mixin providing CRUD operations for Dataset table.

    This mixin handles dataset metadata and information, including dataset
    identifiers, names, and associated metadata. Datasets contain collections
    of sequences used for CSP benchmarking and algorithm evaluation.
    """

    def dataset_create(
        self,
        *,
        dataset_id: str,
        work_id: Optional[str],
        name: Optional[str] = None,
        meta: Dict[str, Any] | None = None,
    ) -> int:
        """Create a new dataset entry and return its auto-generated ID.

        Args:
            dataset_id: Logical identifier for the dataset (not the primary key)
            work_id: ID of the work this dataset belongs to (can be None)
            name: Human-readable name for the dataset
            meta: Additional metadata about the dataset (stored as JSON)

        Returns:
            Auto-generated primary key ID of the created dataset.

        Note:
            The returned ID should be used for creating associated dataset sequences.
        """
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
        """Retrieve a dataset by its auto-generated primary key ID.

        Args:
            id: Auto-generated primary key of the dataset

        Returns:
            Dictionary containing dataset data if found, None otherwise.
            Includes dataset_id, work_id, name, and meta fields.
        """
        from ..models import Dataset

        with self.session_scope() as session:
            dataset = session.query(Dataset).filter(Dataset.id == id).first()
            return dataset.to_dict() if dataset else None

    def dataset_update(self, id: int, **fields: Any) -> None:
        """Update a dataset entry by its auto-generated primary key ID.

        Args:
            id: Auto-generated primary key of the dataset to update
            **fields: Key-value pairs of fields to update. The 'meta'
                     field maps to meta_json in the database.

        Note:
            This operation is idempotent - no error if dataset doesn't exist.
        """
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
        """Delete a dataset entry by its auto-generated primary key ID.

        Args:
            id: Auto-generated primary key of the dataset to delete

        Warning:
            This will also delete all associated dataset sequences due to
            foreign key constraints. This operation is idempotent.
        """
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
        order_desc: bool = False,
    ) -> Tuple[List[Dict[str, Any]], int]:
        """List datasets with filtering and pagination.

        Args:
            filters: Dictionary of field:value pairs for filtering. Supports:
                - Simple values: {'work_id': 'abc', 'dataset_id': 'dataset1'}
                - List values: {'work_id': ['abc', 'def']}
                - Advanced operators: {'name': {'operator': 'like', 'value': 'test'}}
            offset: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            order_by: Field name to order results by
            order_desc: Whether to order in descending order

        Returns:
            Tuple containing:
                - List of dataset dictionaries matching the criteria
                - Total count of matching records (before pagination)

        Note:
            Advanced operators supported: 'like', 'ilike', 'gt', 'gte',
            'lt', 'lte', 'ne' for flexible querying. Useful for searching
            datasets by name patterns or filtering by work association.
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
                        elif isinstance(value, dict) and "operator" in value:
                            # Support for advanced operators
                            op = value["operator"]
                            val = value["value"]
                            column = getattr(Dataset, field)

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
