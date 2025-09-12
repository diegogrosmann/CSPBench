"""Dataset Sequence CRUD operations mixin.

This module provides CRUD operations for the DatasetSequence table,
which stores individual sequences within datasets for CSP benchmark problems.
"""

from typing import Any, Dict, List, Optional, Tuple
from sqlalchemy import and_, or_


class DatasetSequenceCRUDMixin:
    """Mixin providing CRUD operations for DatasetSequence table.
    
    This mixin handles individual sequence data within datasets. Each sequence
    represents a specific problem instance or data point within a larger dataset
    used for CSP benchmarking and algorithm evaluation.
    """

    def dataset_sequence_create(self, *, dataset_id: int, seq_index: int, sequence: str) -> None:
        """Create a new dataset sequence entry.
        
        Args:
            dataset_id: Auto-generated ID of the dataset this sequence belongs to
            seq_index: Index/position of this sequence within the dataset
            sequence: The actual sequence data/content
        """
        from ..models import DatasetSequence
        
        with self.session_scope() as session:
            seq = DatasetSequence(
                dataset_id=dataset_id,
                seq_index=seq_index,
                sequence=sequence
            )
            session.add(seq)

    def dataset_sequence_get(self, dataset_id: int, seq_index: int) -> Optional[Dict[str, Any]]:
        """Retrieve a dataset sequence by dataset ID and sequence index.
        
        Args:
            dataset_id: Auto-generated ID of the dataset
            seq_index: Index of the sequence within the dataset
            
        Returns:
            Dictionary containing sequence data if found, None otherwise.
        """
        from ..models import DatasetSequence
        
        with self.session_scope() as session:
            seq = session.query(DatasetSequence).filter(
                and_(DatasetSequence.dataset_id == dataset_id, DatasetSequence.seq_index == seq_index)
            ).first()
            return seq.to_dict() if seq else None

    def dataset_sequence_update(self, dataset_id: int, seq_index: int, *, sequence: str) -> None:
        """Update a dataset sequence's content.
        
        Args:
            dataset_id: Auto-generated ID of the dataset
            seq_index: Index of the sequence within the dataset
            sequence: New sequence content to update
            
        Note:
            This operation is idempotent - no error if sequence doesn't exist.
        """
        from ..models import DatasetSequence
        
        with self.session_scope() as session:
            seq = session.query(DatasetSequence).filter(
                and_(DatasetSequence.dataset_id == dataset_id, DatasetSequence.seq_index == seq_index)
            ).first()
            if seq:
                seq.sequence = sequence

    def dataset_sequence_delete(self, dataset_id: int, seq_index: int) -> None:
        """Delete a dataset sequence.
        
        Args:
            dataset_id: Auto-generated ID of the dataset
            seq_index: Index of the sequence within the dataset
            
        Note:
            This operation is idempotent - no error if sequence doesn't exist.
        """
        from ..models import DatasetSequence
        
        with self.session_scope() as session:
            seq = session.query(DatasetSequence).filter(
                and_(DatasetSequence.dataset_id == dataset_id, DatasetSequence.seq_index == seq_index)
            ).first()
            if seq:
                session.delete(seq)

    def dataset_sequence_list(
        self,
        *,
        filters: Optional[Dict[str, Any]] = None,
        offset: int = 0,
        limit: Optional[int] = None,
        order_by: Optional[str] = None,
        order_desc: bool = False
    ) -> Tuple[List[Dict[str, Any]], int]:
        """List dataset sequences with filtering and pagination.
        
        Args:
            filters: Dictionary of field:value pairs for filtering. Supports:
                - Simple values: {'dataset_id': 123}
                - List values: {'dataset_id': [123, 456]}
                - Advanced operators: {'seq_index': {'operator': 'gte', 'value': 10}}
            offset: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            order_by: Field name to order results by (typically 'seq_index')
            order_desc: Whether to order in descending order
            
        Returns:
            Tuple containing:
                - List of dataset sequence dictionaries
                - Total count of matching records (before pagination)
                
        Note:
            Advanced operators supported: 'like', 'ilike', 'gt', 'gte', 
            'lt', 'lte', 'ne' for flexible querying. Useful for filtering
            sequences by content patterns or index ranges.
        """
        from ..models import DatasetSequence
        
        with self.session_scope() as session:
            query = session.query(DatasetSequence)
            
            # Apply filters
            if filters:
                for field, value in filters.items():
                    if hasattr(DatasetSequence, field):
                        if isinstance(value, list):
                            query = query.filter(getattr(DatasetSequence, field).in_(value))
                        elif isinstance(value, dict) and 'operator' in value:
                            # Support for advanced operators
                            op = value['operator']
                            val = value['value']
                            column = getattr(DatasetSequence, field)
                            
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
                            query = query.filter(getattr(DatasetSequence, field) == value)
            
            # Get total count before pagination
            total_count = query.count()
            
            # Apply ordering
            if order_by and hasattr(DatasetSequence, order_by):
                column = getattr(DatasetSequence, order_by)
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
            results = [seq.to_dict() for seq in query.all()]
            
            return results, total_count