"""Persistence wrapper classes.

This module exposes persistence wrapper classes with support for a backward-compatible
alias: ``UnitScopedPersistence`` which points to ``ExecutionScopedPersistence``.
The old name is still imported in various places in the code and tests; maintaining
the alias prevents breaking these points while the rest of the refactor progresses.
"""

from .combination_scoped import CombinationScopedPersistence
from .execution_scoped import ExecutionScopedPersistence
from .work_scoped import WorkScopedPersistence

# Backward-compatible alias
UnitScopedPersistence = ExecutionScopedPersistence  # type: ignore

__all__ = [
    "WorkScopedPersistence",
    "ExecutionScopedPersistence",
    "UnitScopedPersistence",  # legacy alias
    "CombinationScopedPersistence",
]
