"""Base work state persistence class with core functionality."""

import json
import logging
import sqlite3
import threading
import time
from pathlib import Path
from typing import Any, Sequence

from .schema import SCHEMA_STATEMENTS
from .utils import validate_required_field


class WorkBase:
    """Core work state persistence with database connection and basic operations."""

    def __init__(self, db_path: Path):
        # Caminho e diretórios
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        # Conexão SQLite
        self._conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row  # rows como mapeamentos
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")

        # Infra de concorrência e estado
        self._lock = threading.RLock()
        self._init_done = False
        self._logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

        # Inicializa schema
        self.init_schema()

    def _execute(self, sql: str, params: Sequence[Any] | None = None) -> None:
        """Execute SQL statement with proper locking."""
        with self._lock:
            self._conn.execute(sql, params or [])
            self._conn.commit()

    def _executemany(self, sql: str, seq_of_params):
        """Execute SQL statement with multiple parameter sets."""
        with self._lock:
            self._conn.executemany(sql, seq_of_params)
            self._conn.commit()

    def _fetch_one(self, sql: str, params: Sequence[Any] | None = None) -> tuple | None:
        """Fetch one row with proper locking."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(sql, params or [])
            return cursor.fetchone()

    def _fetch_all(self, sql: str, params: Sequence[Any] | None = None) -> list[tuple]:
        """Fetch all rows with proper locking."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(sql, params or [])
            return cursor.fetchall()

    def init_schema(self) -> None:
        """Initialize database schema."""
        if self._init_done:
            return

        with self._lock:
            cur = self._conn.cursor()
            for stmt in SCHEMA_STATEMENTS:
                cur.execute(stmt)
            self._conn.commit()
            self._init_done = True

    def backup_database(self, backup_path: Path | None = None) -> Path:
        """Create database backup."""
        if backup_path is None:
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            backup_path = (
                self.db_path.parent / f"backup_{timestamp}_{self.db_path.name}"
            )

        backup_path = Path(backup_path)
        backup_path.parent.mkdir(parents=True, exist_ok=True)

        with self._lock:
            # Use SQLite backup API
            backup_conn = sqlite3.connect(backup_path)
            self._conn.backup(backup_conn)
            backup_conn.close()

        self._logger.info(f"Database backed up to: {backup_path}")
        return backup_path

    def vacuum_database(self) -> None:
        """Optimize database by removing unused space."""
        with self._lock:
            self._conn.execute("VACUUM")
            self._conn.commit()
        self._logger.info("Database vacuumed successfully")

    def get_database_info(self) -> dict[str, Any]:
        """Get database information and statistics."""
        with self._lock:
            cursor = self._conn.cursor()

            # Get database size
            cursor.execute("PRAGMA page_count")
            page_count = cursor.fetchone()[0]
            cursor.execute("PRAGMA page_size")
            page_size = cursor.fetchone()[0]

            # Get table information
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
            tables = [row[0] for row in cursor.fetchall()]

            table_stats = {}
            for table in tables:
                cursor.execute(f"SELECT COUNT(*) FROM {table}")
                count = cursor.fetchone()[0]
                table_stats[table] = count

            return {
                "db_path": str(self.db_path),
                "size_bytes": page_count * page_size,
                "page_count": page_count,
                "page_size": page_size,
                "tables": table_stats,
                "wal_mode": True,
            }
