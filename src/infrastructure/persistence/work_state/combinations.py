"""Combinations mixin for work state persistence.

Adiciona (patch) logs detalhados para diagnosticar perda/skip de combinações:
 - submit_combinations: lista inseridas vs ignoradas
 - update_combination_status: transições (old -> new) com timestamps
 - get_next_queued_combination: seleção e fila remanescente
 - init_combination: combinações reiniciadas

Ativar via LOG_LEVEL=DEBUG para maior verbosidade.
"""

import json
import time
from typing import Any, Optional

from src.infrastructure.persistence.work_state.events import EventsMixin
from src.infrastructure.logging_config import get_logger

from .utils import validate_status

# Logger dedicado para operações de combinação
_combolog = get_logger("CSPBench.Persistence.Combinations")


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

            if rows:
                _combolog.debug(
                    "[init_combination] Reinicializando %d combinações work_id=%s: %s",
                    len(rows),
                    work_id,
                    [
                        {
                            "id": r[0],
                            "task": r[1],
                            "dataset": r[2],
                            "preset": r[3],
                            "algorithm": r[4],
                            "status": r[5],
                        }
                        for r in rows
                    ],
                )
                for cur in rows:
                    # Log estruturado também em events
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
        unique_keys: list[tuple[str, str, str, str, str, str]] = []
        for combo in tasks_combinations:
            mode = combo.get("mode", "experiment")
            rows_params.append(
                (
                    work_id,
                    combo["task_id"],
                    combo["dataset_id"],
                    combo["preset_id"],
                    combo["algorithm_id"],
                    mode,
                    "queued",
                    combo.get("total_sequences", 0),
                    timestamp,
                    None,  # started_at
                    None,  # finished_at
                )
            )
            unique_keys.append(
                (
                    combo["task_id"],
                    combo["dataset_id"],
                    combo["preset_id"],
                    combo["algorithm_id"],
                    mode,
                    str(combo.get("total_sequences", 0)),
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
                SELECT id, task_id, dataset_id, preset_id, algorithm_id, mode, total_sequences
                FROM combinations
                WHERE work_id = ? AND created_at = ?
                """,
                (work_id, timestamp),
            )
            inserted_rows = cursor.fetchall()
            inserted = len(inserted_rows)

            # Diagnóstico de duplicatas ignoradas
            if inserted != len(tasks_combinations):
                _combolog.debug(
                    "[submit_combinations] work_id=%s inseridas=%d ignoradas=%d total_enviadas=%d",
                    work_id,
                    inserted,
                    len(tasks_combinations) - inserted,
                    len(tasks_combinations),
                )
            else:
                _combolog.debug(
                    "[submit_combinations] work_id=%s todas %d combinações inseridas",
                    work_id,
                    inserted,
                )

            if inserted_rows:
                _combolog.debug(
                    "[submit_combinations] IDs inseridos: %s",
                    [
                        {
                            "id": r[0],
                            "task": r[1],
                            "dataset": r[2],
                            "preset": r[3],
                            "algorithm": r[4],
                            "mode": r[5],
                            "total_sequences": r[6],
                        }
                        for r in inserted_rows
                    ],
                )

            # Contagem atual de estados por status para visibilidade
            cursor.execute(
                """
                SELECT status, COUNT(*) FROM combinations WHERE work_id = ? GROUP BY status
                """,
                (work_id,),
            )
            status_counts = {row[0]: row[1] for row in cursor.fetchall()}
            _combolog.debug(
                "[submit_combinations] Distribuição de status após inserção: %s",
                status_counts,
            )

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
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                """
                SELECT id, status, started_at, finished_at FROM combinations
                WHERE work_id=? AND task_id=? AND dataset_id=? AND preset_id=? AND algorithm_id=?
                """,
                (work_id, task_id, dataset_id, preset_id, algorithm_id),
            )
            row = cursor.fetchone()
            orig = None
            if row:
                orig = {
                    "id": row[0],
                    "old_status": row[1],
                    "old_started_at": row[2],
                    "old_finished_at": row[3],
                }

            if status == "running":
                cursor.execute(
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
                cursor.execute(
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
                cursor.execute(
                    """
                    UPDATE combinations 
                    SET status=? 
                    WHERE work_id=? AND task_id=? AND dataset_id=? AND preset_id=? AND algorithm_id=?
                    """,
                    (status, work_id, task_id, dataset_id, preset_id, algorithm_id),
                )
            self._conn.commit()

            # Recarregar estado final para log
            cursor.execute(
                """
                SELECT id, status, started_at, finished_at FROM combinations
                WHERE work_id=? AND task_id=? AND dataset_id=? AND preset_id=? AND algorithm_id=?
                """,
                (work_id, task_id, dataset_id, preset_id, algorithm_id),
            )
            new_row = cursor.fetchone()
            if new_row:
                _combolog.debug(
                    "[update_combination_status] id=%s %s/%s/%s/%s: %s -> %s (started_at=%s finished_at=%s) orig_started_at=%s orig_finished_at=%s",
                    new_row[0],
                    task_id,
                    dataset_id,
                    preset_id,
                    algorithm_id,
                    orig["old_status"] if orig else None,
                    new_row[1],
                    new_row[2],
                    new_row[3],
                    orig["old_started_at"] if orig else None,
                    orig["old_finished_at"] if orig else None,
                )

    def get_next_queued_combination(self, work_id: str) -> dict[str, Any] | None:
        """Get next queued combination for execution following hierarchical order.
        
        Order: task_id -> dataset_id -> preset_id -> algorithm_id
        This ensures proper progression through the hierarchy.
        """
        with self._lock:
            cursor = self._conn.cursor()
            # Contagem antes
            cursor.execute(
                "SELECT COUNT(*) FROM combinations WHERE work_id = ? AND status='queued'",
                (work_id,),
            )
            queued_before = cursor.fetchone()[0]

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
                _combolog.debug(
                    "[get_next_queued_combination] Nenhuma combinação em fila (queued_before=%d) work_id=%s",
                    queued_before,
                    work_id,
                )
                return None

            combo = {
                "id": row[0],
                "task_id": row[1],
                "dataset_id": row[2],
                "preset_id": row[3],
                "algorithm_id": row[4],
                "mode": row[5],
                "total_sequences": row[6],
            }
            _combolog.debug(
                "[get_next_queued_combination] Selecionada id=%s task=%s dataset=%s preset=%s algorithm=%s queued_before=%d",
                combo["id"],
                combo["task_id"],
                combo["dataset_id"],
                combo["preset_id"],
                combo["algorithm_id"],
                queued_before,
            )
            return combo

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
