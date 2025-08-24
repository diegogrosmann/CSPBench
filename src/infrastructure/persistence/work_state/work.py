import time
from typing import Any


class WorkMixin:
    """Core work state persistence with database connection and basic operations."""

    def submit_work(self, work) -> None:
        """Insert or update work item."""
        work_dict = work.to_dict()
        self._execute(
            """
            INSERT INTO work(id, config_json, status, created_at, updated_at, output_path, error, extra_json)
            VALUES(?,?,?,?,?,?,?,?)
            ON CONFLICT(id) DO UPDATE SET
              config_json=excluded.config_json,
              status=excluded.status,
              updated_at=excluded.updated_at,
              output_path=excluded.output_path,
              error=excluded.error,
              extra_json=excluded.extra_json
            """,
            (
                work_dict["id"],
                work_dict["config_json"],
                work_dict["status"],
                work_dict["created_at"],
                work_dict["updated_at"],
                work_dict["output_path"],
                work_dict.get("error", ""),
                work_dict["extra_json"],
            ),
        )

    def get_work_status(self, work_id: str) -> str | None:
        """Get work status from WorkService, replicate to work table and return."""
        from .utils import validate_required_field
        
        validate_required_field(work_id, "work_id")
        
        try:
            # Import here to avoid circular dependencies
            from src.application.services.work_service import get_work_service
            
            # Get status from WorkService
            work_service = get_work_service()
            service_status = work_service.get_status(work_id)
            
            if service_status is None:
                return None
            
            # Update work table with the status from service
            self._execute(
                "UPDATE work SET status=?, updated_at=? WHERE id=?",
                (service_status, time.time(), work_id)
            )
            
            return service_status
            
        except Exception:
            # If WorkService fails, fallback to local database
            with self._lock:
                cursor = self._conn.cursor()
                cursor.execute("SELECT status FROM work WHERE id = ?", (work_id,))
                row = cursor.fetchone()
                return row[0] if row else None

    def update_work_status(self, work_id: str, status: str, **fields: Any) -> None:
        """Update work status with additional fields."""
        from .utils import validate_status, validate_required_field

        validate_required_field(work_id, "work_id")
        validate_status(status)

        # Update local work table
        cols = ["status=?", "updated_at=?"]
        params: list[Any] = [status, time.time()]

        for k, v in fields.items():
            cols.append(f"{k}=?")
            params.append(v)

        params.append(work_id)
        sql = f"UPDATE work SET {', '.join(cols)} WHERE id=?"
        self._execute(sql, params)
        
        # Reflect status change in WorkService
        try:
            from src.application.services.work_service import get_work_service
            from src.domain.status import BaseStatus
            
            work_service = get_work_service()
            
            # Map status to appropriate WorkService method
            if status == BaseStatus.RUNNING.value:
                work_service.mark_running(work_id)
            elif status == BaseStatus.COMPLETED.value:
                work_service.mark_finished(work_id)
            elif status == BaseStatus.FAILED.value:
                error_msg = fields.get("error", "Unknown error")
                work_service.mark_error(work_id, str(error_msg))
            elif status == BaseStatus.PAUSED.value:
                work_service.pause(work_id)
            elif status == BaseStatus.CANCELED.value:
                work_service.cancel(work_id)
            # Note: QUEUED status doesn't have a specific method in WorkService
            
        except Exception:
            # If WorkService update fails, log but don't raise
            # The local database update already succeeded
            pass

    def get_work_statistics(self, work_id: str) -> dict[str, Any]:
        """Get comprehensive work statistics."""
        with self._lock:
            cursor = self._conn.cursor()

            # Basic work info
            cursor.execute(
                "SELECT status, created_at, updated_at FROM work WHERE id = ?",
                (work_id,),
            )
            work_row = cursor.fetchone()

            if not work_row:
                return {"error": "Work not found"}

            # Dataset counts
            cursor.execute(
                """
                SELECT COUNT(DISTINCT dataset_id) as dataset_count
                FROM combinations 
                WHERE work_id = ?
                """,
                (work_id,),
            )
            dataset_count = cursor.fetchone()[0]

            # Algorithm counts
            cursor.execute(
                """
                SELECT COUNT(DISTINCT algorithm_id) as algorithm_count
                FROM combinations 
                WHERE work_id = ?
                """,
                (work_id,),
            )
            algorithm_count = cursor.fetchone()[0]

            # Combination statistics
            cursor.execute(
                """
                SELECT status, COUNT(*) as count
                FROM combinations 
                WHERE work_id = ?
                GROUP BY status
                """,
                (work_id,),
            )
            combination_stats = {row[0]: row[1] for row in cursor.fetchall()}

            # Execution statistics
            cursor.execute(
                """
                SELECT e.status, COUNT(*) as count
                FROM executions e
                JOIN combinations pc ON e.combination_id = pc.id
                WHERE pc.work_id = ?
                GROUP BY e.status
                """,
                (work_id,),
            )
            execution_stats = {row[0]: row[1] for row in cursor.fetchall()}

            # Event statistics
            cursor.execute(
                """
                SELECT event_type, COUNT(*) as count
                FROM events 
                WHERE work_id = ?
                GROUP BY event_type
                """,
                (work_id,),
            )
            event_stats = {row[0]: row[1] for row in cursor.fetchall()}

            # Calculate runtime
            current_time = time.time()
            runtime_seconds = current_time - work_row[1]  # created_at
            last_activity = current_time - work_row[2]  # updated_at

            return {
                "work_id": work_id,
                "status": work_row[0],
                "created_at": work_row[1],
                "updated_at": work_row[2],
                "runtime_seconds": runtime_seconds,
                "last_activity_seconds": last_activity,
                "dataset_count": dataset_count,
                "algorithm_count": algorithm_count,
                "combination_stats": combination_stats,
                "execution_stats": execution_stats,
                "event_stats": event_stats,
                "total_combinations": sum(combination_stats.values()),
                "total_executions": sum(execution_stats.values()),
                "total_events": sum(event_stats.values()),
            }
