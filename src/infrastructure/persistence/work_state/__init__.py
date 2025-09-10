"""
Work state persistence module.

Provides comprehensive work state management using SQLAlchemy ORM:
- Core database operations with SQLAlchemy
- Event logging and retrieval
- Execution tracking
- Pipeline management
- Queries and reporting
- Convenient wrapper classes

Migration Note:
This module has been migrated from custom database drivers to SQLAlchemy ORM
for better maintainability, type safety, and database compatibility.
"""

from .core import WorkPersistence
from .models import (
    Work,
    Dataset,
    DatasetSequence,
    Combination,
    Execution,
    ExecutionProgress,
    Event,
)
from .wrappers import (
    ExecutionScopedPersistence,
    CombinationScopedPersistence,
    WorkScopedPersistence,
)
from .mixins import (
    ExecutionDetail,
    ProgressSummary,
    ErrorSummary,
)

__all__ = [
    "WorkPersistence",
    "SQLAlchemyWorkPersistence",
    "WorkScopedPersistence",
    "ExecutionScopedPersistence",
    "CombinationScopedPersistence",
    # SQLAlchemy models
    "Work",
    "Dataset",
    "DatasetSequence",
    "Combination",
    "Execution",
    "ExecutionProgress",
    "Event",
    # Query result classes
    "ExecutionDetail",
    "ProgressSummary",
    "ErrorSummary",
]
