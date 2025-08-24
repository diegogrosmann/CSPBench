"""Main work state persistence class composed from mixins."""

from pathlib import Path

from src.infrastructure.persistence.work_state.combinations import CombinationsMixin
from src.infrastructure.persistence.work_state.dataset import DatasetMixin
from src.infrastructure.persistence.work_state.executions import ExecutionsMixin
from src.infrastructure.persistence.work_state.events import EventsMixin
from src.infrastructure.persistence.work_state.work import WorkMixin
from src.infrastructure.persistence.work_state.base import WorkBase


class WorkStatePersistence(
    WorkBase, WorkMixin, EventsMixin, ExecutionsMixin, CombinationsMixin, DatasetMixin
):
    """
    Main work state persistence class composed from focused mixins.

    Provides comprehensive work state management including:
    - Core database operations (WorkStatePersistenceBase)
    - Event logging and retrieval (EventsMixin)
    - Execution tracking (ExecutionsMixin)
    - Combinations management (CombinationsMixin)
    - Dataset management (DatasetMixin)
    """

    def __init__(self, db_path: Path):
        """Initialize work state persistence with all functionality."""
        super().__init__(db_path)
        self._logger.info(f"WorkStatePersistence initialized with database: {db_path}")

    def __repr__(self) -> str:
        return f"WorkStatePersistence(db_path={self.db_path})"
