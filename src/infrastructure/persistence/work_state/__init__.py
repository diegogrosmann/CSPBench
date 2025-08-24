"""
Work state persistence module.

Provides comprehensive work state management with organized, focused components:
- Core database operations
- Event logging and retrieval
- Execution tracking
- Pipeline management
- Queries and reporting
- Convenient wrapper classes
"""

from .core import WorkStatePersistence
from .wrappers import (
    WorkScopedPersistence,
    ExecutionScopedPersistence,  # nome novo
    UnitScopedPersistence,  # alias legado
    CombinationScopedPersistence,
)

__all__ = [
    "WorkStatePersistence",
    "WorkScopedPersistence",
    "ExecutionScopedPersistence",
    "UnitScopedPersistence",  # legado
    "CombinationScopedPersistence",
]
