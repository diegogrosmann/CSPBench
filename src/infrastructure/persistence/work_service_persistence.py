"""
Persistent Work Repository Implementation
Provides SQLite-based storage for WorkManager state with global persistence.
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
import threading
import time
from pathlib import Path
from typing import List, Optional

from application.ports.repositories import WorkRepository
from src.domain.work import WorkItem, WorkStatus

logger = logging.getLogger(__name__)


class WorkServicePersistence(WorkRepository):
    """
    SQLite-based persistent work repository for global WorkManager state.
    Ensures data persists across FastAPI requests and application restarts.
    """

    def __init__(self):
        """
        Initialize persistent repository with SQLite database.
        
        Args:
            db_path: Path to SQLite database file. If None, uses WORK_DB_PATH environment variable
                    or defaults to "data/work_manager.db"
        """
        if db_path is None:
            db_path = os.getenv("WORK_DB_PATH", "data/work_manager.db")
        
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        self._lock = threading.RLock()
        self._conn = sqlite3.connect(
            str(self.db_path), 
            check_same_thread=False,
            timeout=30.0
        )
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.execute("PRAGMA foreign_keys=ON")
        
        self._init_schema()

    def _init_schema(self) -> None:
        """Initialize database schema for work items."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS work_items (
                    id TEXT PRIMARY KEY,
                    status TEXT NOT NULL,
                    created_at REAL NOT NULL,
                    updated_at REAL NOT NULL,
                    output_path TEXT,
                    extra_json TEXT DEFAULT '{}',     
                    error TEXT                  
                )
            """)
            
            self._conn.commit()

    def _execute(self, sql: str, params: tuple = ()) -> sqlite3.Cursor:
        """Execute SQL with proper locking."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(sql, params)
            self._conn.commit()
            return cursor

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
        extra_dict = json.loads(row["extra_json"] or "{}")
        
        # Keep batch_file in extra dict if it exists in the row
        if "batch_file" in row.keys() and row["batch_file"]:
            extra_dict["batch_file"] = row["batch_file"]
        
        return WorkItem(
            id=row["id"],
            config=None,  # No longer store config in database
            status=WorkStatus(row["status"]),
            created_at=row["created_at"],
            updated_at=row["updated_at"],
            output_path=row["output_path"],
            error=row["error"],
            extra=extra_dict
        )

    def _work_item_to_params(self, item: WorkItem) -> tuple:
        """Convert WorkItem to database parameters."""
        # Keep batch_file in extra dict, don't extract it
        extra_dict = item.extra.copy() if item.extra else {}
        batch_file = extra_dict.get("batch_file")  # Get but don't remove
        
        return (
            item.id,
            item.status.value,
            item.created_at,
            item.updated_at,
            item.output_path,
            batch_file,
            json.dumps(extra_dict, ensure_ascii=False),
            item.error
        )

    def add(self, item: WorkItem) -> None:
        """Add new work item to repository."""
        params = self._work_item_to_params(item)
        self._execute("""
            INSERT INTO work_items 
            (id, status, created_at, updated_at, output_path, extra_json, error)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, params)

    def get(self, work_id: str) -> Optional[WorkItem]:
        """Get work item by ID."""
        row = self._fetch_one(
            "SELECT * FROM work_items WHERE id = ?", 
            (work_id,)
        )
        return self._row_to_work_item(row) if row else None

    def update(self, item: WorkItem) -> None:
        """Update existing work item."""
        params = self._work_item_to_params(item)
        # Remove id from params for UPDATE (it's in WHERE clause)
        update_params = params[1:] + (params[0],)
        
        self._execute("""
            UPDATE work_items 
            SET status=?, created_at=?, updated_at=?, output_path=?, error=?, 
                extra_json=?
            WHERE id=?
        """, update_params)

    def remove(self, work_id: str) -> bool:
        """Remove work item from repository."""
        cursor = self._execute(
            "DELETE FROM work_items WHERE id = ?", 
            (work_id,)
        )
        return cursor.rowcount > 0

    def list(self) -> List[WorkItem]:
        """List all work items ordered by creation time."""
        rows = self._fetch_all("""
            SELECT * FROM work_items 
            ORDER BY created_at DESC
        """)
        return [self._row_to_work_item(row) for row in rows]

    def list_by_status(self, status: WorkStatus) -> List[WorkItem]:
        """List work items by status."""
        rows = self._fetch_all("""
            SELECT * FROM work_items 
            WHERE status = ?
            ORDER BY created_at DESC
        """, (status.value,))
        return [self._row_to_work_item(row) for row in rows]

    def count(self) -> int:
        """Count total work items."""
        row = self._fetch_one("SELECT COUNT(*) as total FROM work_items")
        return row["total"] if row else 0

    def count_by_status(self, status: WorkStatus) -> int:
        """Count work items by status."""
        row = self._fetch_one(
            "SELECT COUNT(*) as total FROM work_items WHERE status = ?",
            (status.value,)
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
        cursor = self._execute("""
            DELETE FROM work_items 
            WHERE updated_at < ? 
            AND status IN ('completed', 'failed', 'cancelled')
        """, (cutoff_time,))
        return cursor.rowcount

    def close(self) -> None:
        """Close database connection."""
        with self._lock:
            if self._conn:
                self._conn.close()

    def __del__(self):
        """Cleanup on destruction."""
        try:
            self.close()
        except:
            pass  # Ignore errors during cleanup
