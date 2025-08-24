"""Combinations mixin for work state persistence."""

import json
import time
from typing import Any, Optional

from src.infrastructure.persistence.work_state.events import EventsMixin

from .utils import validate_status


class CombinationsMixin:
    """Mixin providing combinations management functionality."""

    def init_combination(self, work_id: str) -> bool:
        """
        Reinicia todas as combinações em andamento, pausadas ou canceladas para 'queued'.
        Lista e loga as combinações afetadas antes de reiniciar.
        Retorna True se já existem combinações, False caso contrário.
        """
        event = EventsMixin()

        with self._lock:
            cursor = self._conn.cursor()
            # Buscar combinações em andamento, pausadas ou canceladas
            cursor.execute(
                """
                SELECT id, task_id, dataset_id, preset_id, algorithm_id, status
                FROM combinations
                WHERE work_id = ? AND status IN ('running', 'paused', 'canceled')
                """,
                (work_id,),
            )

            rows = cursor.fetchall()
            has_combinations = bool(rows)

            for cur in rows:
                # Logar cada combinação afetada
                self.combination_warning(
                    work_id,
                    combination_id=cur[0],
                    message=f"Reiniciando combinação {cur[0]} com status {cur[5]}",
                )

            if has_combinations:
                # Reiniciar todas para 'queued'
                cursor.execute(
                    """
                    UPDATE combinations
                    SET status = 'queued', started_at = NULL, finished_at = NULL
                    WHERE work_id = ? AND status IN ('running', 'paused', 'canceled')
                    """,
                    (work_id,),
                )
                self._conn.commit()

            return has_combinations

    def submit_combinations(
        self, work_id: str, tasks_combinations: list[dict[str, Any]]
    ) -> int:
        """Insere (ou ignora se já existir) TODAS as combinações fornecidas.

        Retorna quantas novas combinações foram de fato inseridas (count)."""
        if not tasks_combinations:
            return 0

        timestamp = time.time()
        rows_params = []
        for combo in tasks_combinations:
            rows_params.append(
                (
                    work_id,
                    combo["task_id"],
                    combo["dataset_id"],
                    combo["preset_id"],
                    combo["algorithm_id"],
                    combo.get("mode", "experiment"),
                    "queued",
                    combo.get("total_sequences", 0),
                    timestamp,
                    None,  # started_at
                    None,  # finished_at
                )
            )

        with self._lock:
            cursor = self._conn.cursor()
            cursor.executemany(
                """
                INSERT OR IGNORE INTO combinations(
                    work_id, task_id, dataset_id, preset_id, algorithm_id, mode,
                    status, total_sequences, created_at, started_at, finished_at
                ) VALUES(?,?,?,?,?,?,?,?,?,?,?)
                """,
                rows_params,
            )
            self._conn.commit()

            # Contar quantas foram efetivamente inseridas agora (status queued e created_at ~= timestamp)
            cursor.execute(
                """
                SELECT COUNT(*) FROM combinations
                WHERE work_id = ? AND created_at = ?
                """,
                (work_id, timestamp),
            )
            inserted = cursor.fetchone()[0]
            return inserted

    def update_combination_status(
        self,
        work_id: str,
        task_id: str,
        dataset_id: str,
        preset_id: str,
        algorithm_id: str,
        status: str,
    ) -> None:
        """Update status of specific combination with propagation."""
        validate_status(status)

        timestamp = time.time()

        if status == "running":
            self._execute(
                """
                UPDATE combinations 
                SET status=?, started_at=? 
                WHERE work_id=? AND task_id=? AND dataset_id=? AND preset_id=? AND algorithm_id=?
                """,
                (
                    status,
                    timestamp,
                    work_id,
                    task_id,
                    dataset_id,
                    preset_id,
                    algorithm_id,
                ),
            )
        elif status in ("completed", "failed"):
            self._execute(
                """
                UPDATE combinations 
                SET status=?, finished_at=? 
                WHERE work_id=? AND task_id=? AND dataset_id=? AND preset_id=? AND algorithm_id=?
                """,
                (
                    status,
                    timestamp,
                    work_id,
                    task_id,
                    dataset_id,
                    preset_id,
                    algorithm_id,
                ),
            )
        else:
            self._execute(
                """
                UPDATE combinations 
                SET status=? 
                WHERE work_id=? AND task_id=? AND dataset_id=? AND preset_id=? AND algorithm_id=?
                """,
                (status, work_id, task_id, dataset_id, preset_id, algorithm_id),
            )

    def get_next_queued_combination(self, work_id: str) -> dict[str, Any] | None:
        """Get next queued combination for execution following hierarchical order.
        
        Order: task_id -> dataset_id -> preset_id -> algorithm_id
        This ensures proper progression through the hierarchy.
        """
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                """
                SELECT id, task_id, dataset_id, preset_id, algorithm_id, mode, total_sequences
                FROM combinations 
                WHERE work_id = ? AND status = 'queued'
                ORDER BY task_id ASC, dataset_id ASC, preset_id ASC, algorithm_id ASC, created_at ASC
                LIMIT 1
                """,
                (work_id,),
            )

            row = cursor.fetchone()

            if not row:
                return None

            return {
                "id": row[0],
                "task_id": row[1],
                "dataset_id": row[2],
                "preset_id": row[3],
                "algorithm_id": row[4],
                "mode": row[5],
                "total_sequences": row[6],
            }

    def get_next_pending_combination(self, work_id: str) -> dict[str, Any] | None:
        """Get next pending combination (compatibility method)."""
        return self.get_next_queued_combination(work_id)

    def get_combination_progress(self, combination_id: int) -> float | None:
        """
        Get combination progress based on execution status.

        Calculates the average progress of all executions in the combination.
        Progress is calculated as: completed executions / total executions.

        Args:
            combination_id: ID of the combination

        Returns:
            Float between 0.0 and 1.0 representing progress, or None if combination not found
        """
        with self._lock:
            cursor = self._conn.cursor()

            # First check if combination exists
            cursor.execute(
                """
                SELECT id FROM combinations
                WHERE id = ?
                """,
                (combination_id,),
            )

            if not cursor.fetchone():
                return None

            # Get execution statistics
            cursor.execute(
                """
                SELECT 
                    COUNT(*) as total_executions,
                    COUNT(CASE WHEN status = 'completed' THEN 1 END) as completed_executions
                FROM executions
                WHERE combination_id = ?
                """,
                (combination_id,),
            )

            stats_row = cursor.fetchone()

            if not stats_row or stats_row[0] == 0:
                # No executions found, progress is 0
                return 0.0

            total_executions = stats_row[0]
            completed_executions = stats_row[1]

            # Calculate progress as ratio of completed executions
            progress = completed_executions / total_executions

            return round(progress, 4)

    def get_combinations(
        self,
        work_id: Optional[str] = None,
        status: Optional[str] = None,
        task_id: Optional[str] = None,
        dataset_id: Optional[str] = None,
        preset_id: Optional[str] = None,
        algorithm_id: Optional[str] = None,
        mode: Optional[str] = None,
        order_by: str = "created_at",
        order_direction: str = "ASC",
        limit: Optional[int] = None,
        offset: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """
        Get combinations with flexible filtering options.

        Args:
            work_id: Filter by work ID
            status: Filter by status
            task_id: Filter by task ID
            dataset_id: Filter by dataset ID
            preset_id: Filter by preset ID
            algorithm_id: Filter by algorithm ID
            mode: Filter by mode
            order_by: Column to order by (default: created_at)
            order_direction: ASC or DESC (default: ASC)
            limit: Maximum number of results
            offset: Number of results to skip

        Returns:
            List of combination dictionaries
        """
        with self._lock:
            cursor = self._conn.cursor()

            # Build WHERE clause dynamically
            where_conditions = []
            params = []

            if work_id is not None:
                where_conditions.append("work_id = ?")
                params.append(work_id)

            if status is not None:
                where_conditions.append("status = ?")
                params.append(status)

            if task_id is not None:
                where_conditions.append("task_id = ?")
                params.append(task_id)

            if dataset_id is not None:
                where_conditions.append("dataset_id = ?")
                params.append(dataset_id)

            if preset_id is not None:
                where_conditions.append("preset_id = ?")
                params.append(preset_id)

            if algorithm_id is not None:
                where_conditions.append("algorithm_id = ?")
                params.append(algorithm_id)

            if mode is not None:
                where_conditions.append("mode = ?")
                params.append(mode)

            # Build query
            query = """
                SELECT id, work_id, task_id, dataset_id, preset_id, algorithm_id, 
                       mode, status, total_sequences, created_at, started_at, finished_at
                FROM combinations
            """

            if where_conditions:
                query += " WHERE " + " AND ".join(where_conditions)

            # Add ordering
            valid_columns = {
                "id",
                "work_id",
                "task_id",
                "dataset_id",
                "preset_id",
                "algorithm_id",
                "mode",
                "status",
                "total_sequences",
                "created_at",
                "started_at",
                "finished_at",
            }

            if order_by in valid_columns:
                direction = "DESC" if order_direction.upper() == "DESC" else "ASC"
                query += f" ORDER BY {order_by} {direction}"

            # Add limit and offset
            if limit is not None:
                query += " LIMIT ?"
                params.append(limit)

                if offset is not None:
                    query += " OFFSET ?"
                    params.append(offset)

            cursor.execute(query, params)
            rows = cursor.fetchall()

            # Convert to dictionaries
            combinations = []
            for row in rows:
                combinations.append(
                    {
                        "id": row[0],
                        "work_id": row[1],
                        "task_id": row[2],
                        "dataset_id": row[3],
                        "preset_id": row[4],
                        "algorithm_id": row[5],
                        "mode": row[6],
                        "status": row[7],
                        "total_sequences": row[8],
                        "created_at": row[9],
                        "started_at": row[10],
                        "finished_at": row[11],
                    }
                )

            return combinations
