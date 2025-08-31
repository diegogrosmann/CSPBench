"""
Persistent Work Repository Implementation
Provides SQLite-based storage for WorkManager state with global persistence.
"""

from __future__ import annotations

import logging
import os
import sqlite3
import threading
import time
from pathlib import Path
from typing import List, Optional

from src.domain.work import WorkItem, WorkStatus
from src.domain.status import BaseStatus
from src.infrastructure.utils.path_utils import get_work_db_path

logger = logging.getLogger(__name__)

_SCHEMA = [
    """
        CREATE TABLE IF NOT EXISTS work_items (
            id TEXT PRIMARY KEY,
            status TEXT CHECK(status IN ('queued', 'running', 'paused', 'canceled', 'completed', 'failed', 'error')),
            created_at REAL NOT NULL,
            updated_at REAL NOT NULL,
            output_path TEXT,
            config_json TEXT,
            extra_json TEXT DEFAULT '{}',     
            error TEXT                  
        )
    """
]


class WorkServicePersistence:
    """
    SQLite-based persistent work repository for global WorkManager state.
    Ensures data persists across FastAPI requests and application restarts.
    """

    def __init__(self, db_path: Optional[str] = None):
        """Initialize with database path from environment or parameter."""
        if db_path is None:
            self.db_path = get_work_db_path()
        else:
            self.db_path = Path(db_path)
            # Convert to absolute path if relative
            if not self.db_path.is_absolute():
                self.db_path = self.db_path.resolve()
            # Create parent directory if needed
            self.db_path.parent.mkdir(parents=True, exist_ok=True)

        self._conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        
        # Configurações para concorrência e performance
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.execute("PRAGMA busy_timeout=30000")  # 30 segundos
        self._conn.execute("PRAGMA wal_autocheckpoint=1000")
        self._conn.execute("PRAGMA cache_size=10000")
        
        self._lock = threading.RLock()
        self._init_done = False

        self._init_schema()

    def _init_schema(self) -> None:
        """Initialize database schema for work items."""

        if self._init_done:
            return
        with self._lock:
            cur = self._conn.cursor()
            for stmt in _SCHEMA:
                cur.execute(stmt)
            self._conn.commit()
            self._init_done = True

    def is_ready(self) -> bool:
        """
        Verifica se o banco de dados está pronto para uso.
        
        Returns:
            True se o banco está inicializado e funcionando
        """
        try:
            with self._lock:
                cursor = self._conn.cursor()
                cursor.execute("SELECT 1")
                return self._init_done
        except Exception:
            return False

    def _execute(self, sql: str, params: tuple = ()) -> sqlite3.Cursor:
        """Execute SQL with proper locking and retry logic."""
        max_retries = 3
        retry_delay = 0.1
        
        for attempt in range(max_retries):
            try:
                with self._lock:
                    cursor = self._conn.cursor()
                    cursor.execute(sql, params)
                    self._conn.commit()
                    return cursor
            except sqlite3.OperationalError as e:
                if "database is locked" in str(e).lower() and attempt < max_retries - 1:
                    logger.warning(f"Database locked, retrying in {retry_delay}s (attempt {attempt + 1}/{max_retries})")
                    time.sleep(retry_delay)
                    retry_delay *= 2
                    continue
                raise
            except Exception:
                raise

    def _fetch_one(self, sql: str, params: tuple = ()) -> Optional[sqlite3.Row]:
        """Fetch one row with proper locking."""
        with self._lock:
            self._conn.row_factory = sqlite3.Row
            cursor = self._conn.cursor()
            cursor.execute(sql, params)
            return cursor.fetchone()

    def _fetch_all(self, sql: str, params: tuple = ()) -> List[sqlite3.Row]:
        """Fetch all rows with proper locking."""
        with self._lock:
            self._conn.row_factory = sqlite3.Row
            cursor = self._conn.cursor()
            cursor.execute(sql, params)
            return cursor.fetchall()

    def _row_to_work_item(self, row: sqlite3.Row) -> WorkItem:
        """Convert database row to WorkItem."""

        return WorkItem.from_dict(
            {
                "id": row["id"],
                "config_json": row["config_json"],
                "status": row["status"],
                "created_at": row["created_at"],  # Already float from database
                "updated_at": row["updated_at"],  # Already float from database
                "output_path": row["output_path"],
                "error": row["error"],
                "extra_json": row["extra_json"],
            }
        )

    def add(self, item: WorkItem) -> None:
        """Add new work item to repository."""
        item_dict = item.to_dict()
        self._execute(
            """
            INSERT INTO work_items 
            (id, config_json, status, created_at, updated_at, output_path, extra_json, error)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
            (
                item_dict["id"],
                item_dict["config_json"],
                item_dict["status"],
                item_dict["created_at"],  # Keep as float for database storage
                item_dict["updated_at"],  # Keep as float for database storage
                item_dict["output_path"],
                item_dict["extra_json"],
                item_dict["error"] or "",
            ),
        )

    def get(self, work_id: str) -> Optional[WorkItem]:
        """Get work item by ID."""
        row = self._fetch_one("SELECT * FROM work_items WHERE id = ?", (work_id,))
        return self._row_to_work_item(row) if row else None

    def update(self, item: WorkItem) -> None:
        """Update existing work item."""
        item_dict = item.to_dict()
        self._execute(
            """
            UPDATE work_items 
            SET status=?, updated_at=?, extra_json=?, error=?
            WHERE id=?
        """,
            (
                item_dict["status"],
                item_dict["updated_at"],  # Keep as float for database storage
                item_dict["extra_json"],
                item_dict["error"],
                item_dict["id"],
            ),
        )

    def remove(self, work_id: str) -> bool:
        """Remove work item from repository."""
        cursor = self._execute("DELETE FROM work_items WHERE id = ?", (work_id,))
        return cursor.rowcount > 0

    def list(self) -> List[WorkItem]:
        """List all work items ordered by creation time."""
        rows = self._fetch_all(
            """
            SELECT * FROM work_items 
            ORDER BY created_at DESC
        """
        )
        return [self._row_to_work_item(row) for row in rows]

    def list_by_status(self, status: BaseStatus) -> List[WorkItem]:
        """List work items by status."""
        rows = self._fetch_all(
            """
            SELECT * FROM work_items 
            WHERE status = ?
            ORDER BY created_at DESC
        """,
            (status.value,),
        )
        return [self._row_to_work_item(row) for row in rows]

    def count(self) -> int:
        """Count total work items."""
        row = self._fetch_one("SELECT COUNT(*) as total FROM work_items")
        return row["total"] if row else 0

    def count_by_status(self, status: BaseStatus) -> int:
        """Count work items by status."""
        row = self._fetch_one(
            "SELECT COUNT(*) as total FROM work_items WHERE status = ?", (status.value,)
        )
        return row["total"] if row else 0

    def cleanup_old_items(self, older_than_days: int = 30) -> int:
        """
        Cleanup old completed/failed work items.

        Args:
            older_than_days: Remove items older than this many days

        Returns:
            Number of items removed
        """
        cutoff_time = time.time() - (older_than_days * 24 * 60 * 60)
        cursor = self._execute(
            """
            DELETE FROM work_items 
            WHERE updated_at < ? 
            AND status IN (?, ?, ?)
        """,
            (
                cutoff_time,
                BaseStatus.COMPLETED.value,
                BaseStatus.FAILED.value,
                BaseStatus.CANCELED.value,
            ),
        )
        return cursor.rowcount

    def get_status_counts(self) -> dict[str, int]:
        """Get count of work items by status."""
        rows = self._fetch_all(
            """
            SELECT status, COUNT(*) as count 
            FROM work_items 
            GROUP BY status
        """
        )
        return {row["status"]: row["count"] for row in rows}

    def close(self) -> None:
        """Close database connection."""
        if hasattr(self, "_conn") and self._conn:
            self._conn.close()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()