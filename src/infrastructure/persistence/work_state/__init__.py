"""
Work state persistence module.

This module provides comprehensive work state management using SQLAlchemy ORM,
offering robust database operations for CSPBench's execution tracking and
result management needs.

Key Components:
    - Core database operations with SQLAlchemy ORM
    - Event logging and retrieval for audit trails
    - Execution tracking with progress monitoring
    - Pipeline management for complex workflows
    - Queries and reporting capabilities
    - Convenient wrapper classes for scoped operations

Features:
    - Thread-safe database operations
    - Automatic schema management
    - Support for multiple database backends
    - Comprehensive error handling and logging
    - Type-safe ORM models with validation
    - Migration utilities for easy upgrades

Migration Note:
    This module has been migrated from custom database drivers to SQLAlchemy ORM
    for better maintainability, type safety, and database compatibility. The public
    API remains largely unchanged, ensuring backward compatibility while providing
    improved reliability and performance.

Usage Example:
    ```python
    from src.infrastructure.persistence.work_state import WorkPersistence

    # Create persistence instance
    persistence = WorkPersistence()

    # Create a new work
    persistence.work_create(
        id="benchmark-001",
        config={"algorithms": ["greedy", "genetic"], "datasets": ["dna1"]}
    )

    # Query work status
    work = persistence.work_get("benchmark-001")
    print(f"Work status: {work['status']}")
    ```
"""

from .core import WorkPersistence
from .mixins import (
    ErrorSummary,
    ExecutionDetail,
    ProgressSummary,
)
from .models import (
    Combination,
    Dataset,
    DatasetSequence,
    Event,
    Execution,
    ExecutionProgress,
    Work,
)
from .wrappers import (
    CombinationScopedPersistence,
    ExecutionScopedPersistence,
    WorkScopedPersistence,
)

# Alias for backward compatibility
SQLAlchemyWorkPersistence = WorkPersistence

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
