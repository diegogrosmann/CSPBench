"""CRUD mixins package for work state persistence.

This module provides a collection of mixin classes that implement CRUD operations
for different entities in the work state persistence layer. Each mixin focuses
on operations for a specific database table or domain concept.

The mixins can be combined to create a complete persistence layer by multiple
inheritance, providing a modular and maintainable approach to database operations.
"""

from .combination_mixin import CombinationCRUDMixin
from .dataset_mixin import DatasetCRUDMixin
from .dataset_sequence_mixin import DatasetSequenceCRUDMixin
from .event_mixin import EventCRUDMixin
from .execution_mixin import ExecutionCRUDMixin
from .execution_progress_mixin import ExecutionProgressCRUDMixin
from .queries_mixin import ErrorSummary, ExecutionDetail, ProgressSummary, QueriesMixin
from .work_mixin import WorkCRUDMixin

__all__ = [
    "WorkCRUDMixin",
    "CombinationCRUDMixin",
    "DatasetCRUDMixin",
    "DatasetSequenceCRUDMixin",
    "EventCRUDMixin",
    "ExecutionCRUDMixin",
    "ExecutionProgressCRUDMixin",
    "QueriesMixin",
    "ExecutionDetail",
    "ProgressSummary",
    "ErrorSummary",
]
