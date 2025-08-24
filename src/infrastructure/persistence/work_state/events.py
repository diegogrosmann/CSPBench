"""Events mixin for work state persistence."""

import json
import time
import traceback
from typing import Any

from .utils import (
    validate_event_type,
    validate_event_category,
    add_timestamp_if_missing,
)


class EventsMixin:
    """Mixin providing event logging functionality."""

    def _log_event(
        self,
        work_id: str,
        event_type: str,
        event_category: str,
        entity_data: dict[str, Any],
    ) -> None:
        """Internal method to log events."""
        validate_event_type(event_type)
        validate_event_category(event_category)

        self._execute(
            "INSERT INTO events(work_id, event_type, event_category, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (
                work_id,
                event_type,
                event_category,
                json.dumps(entity_data, ensure_ascii=False),
                time.time(),
            ),
        )

    # === WORK EVENTS ===
    def work_warning(
        self, work_id: str, message: str, context: dict[str, Any] | None = None
    ) -> None:
        """Log work warning."""
        data = {"message": message, "context": context or {}}
        self._log_event(work_id, "warning", "work", data)

    def work_error(self, work_id: str, error: Exception) -> None:
        """Log work error."""
        error_data = {
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }
        self._log_event(work_id, "error", "work", error_data)

    # === TASK EVENTS ===
    def task_warning(
        self,
        work_id: str,
        task_id: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log task warning."""
        data = {"task_id": task_id, "message": message, "context": context or {}}
        self._log_event(work_id, "warning", "task", data)

    def task_error(self, work_id: str, task_id: str, error: Exception) -> None:
        """Log task error."""
        error_data = {
            "task_id": task_id,
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }
        self._log_event(work_id, "error", "task", error_data)

    # === DATASET EVENTS ===
    def dataset_warning(
        self,
        work_id: str,
        dataset_id: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log dataset warning."""
        data = {"dataset_id": dataset_id, "message": message, "context": context or {}}
        self._log_event(work_id, "warning", "dataset", data)

    def dataset_error(self, work_id: str, dataset_id: str, error: Exception) -> None:
        """Log dataset error."""
        error_data = {
            "dataset_id": dataset_id,
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }
        self._log_event(work_id, "error", "dataset", error_data)

    # === PRESET EVENTS ===
    def preset_warning(
        self,
        work_id: str,
        preset_id: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log preset warning."""
        data = {"preset_id": preset_id, "message": message, "context": context or {}}
        self._log_event(work_id, "warning", "preset", data)

    def preset_error(self, work_id: str, preset_id: str, error: Exception) -> None:
        """Log preset error."""
        error_data = {
            "preset_id": preset_id,
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }
        self._log_event(work_id, "error", "preset", error_data)

    # === COMBINATION EVENTS ===
    def combination_warning(
        self,
        work_id: str,
        combination_id: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log combination warning."""
        data = {
            "combination_id": combination_id,
            "message": message,
            "context": context or {},
        }
        self._log_event(work_id, "warning", "combination", data)

    def combination_error(
        self, work_id: str, combination_id: str, error: Exception
    ) -> None:
        """Log combination error."""
        error_data = {
            "combination_id": combination_id,
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }
        self._log_event(work_id, "error", "combination", error_data)

    # === UNIT EVENTS ===
    def unit_warning(
        self,
        work_id: str,
        unit_id: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log unit warning."""
        data = {"unit_id": unit_id, "message": message, "context": context or {}}
        self._log_event(work_id, "warning", "unit", data)

    def unit_error(self, work_id: str, unit_id: str, error: Exception) -> None:
        """Log unit error."""
        error_data = {
            "unit_id": unit_id,
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }
        self._log_event(work_id, "error", "unit", error_data)

    # === GENERIC EVENTS ===
    def generic_event(
        self,
        work_id: str,
        unit_id: str,
        event_type: str,
        message: str,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Log generic event."""
        data = {"unit_id": unit_id, "message": message, "context": context or {}}
        self._log_event(work_id, event_type, "other", data)

    # === EVENT QUERIES ===
    def get_events(
        self,
        work_id: str,
        event_category: str | None = None,
        event_type: str | None = None,
        unit_id: str | None = None,
        limit: int = 100,
    ) -> list[dict[str, Any]]:
        """Get events with flexible filters."""
        conditions = ["work_id = ?"]
        params = [work_id]

        if event_category:
            conditions.append("event_category = ?")
            params.append(event_category)

        if event_type:
            conditions.append("event_type = ?")
            params.append(event_type)

        if unit_id:
            conditions.append("JSON_EXTRACT(entity_data_json, '$.unit_id') = ?")
            params.append(unit_id)

        where_clause = " AND ".join(conditions)
        params.append(limit)

        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                f"""
                SELECT id, work_id, event_type, event_category, 
                       entity_data_json, timestamp
                FROM events 
                WHERE {where_clause}
                ORDER BY timestamp DESC 
                LIMIT ?
                """,
                params,
            )

            results = []
            for row in cursor.fetchall():
                entity_data = json.loads(row[4])
                results.append(
                    {
                        "id": row[0],
                        "work_id": row[1],
                        "unit_id": entity_data.get("unit_id"),
                        "event_type": row[2],
                        "event_category": row[3],
                        "entity_data": entity_data,
                        "timestamp": row[5],
                    }
                )

            return results

    def get_events_by_category(
        self, work_id: str, event_category: str, limit: int = 100
    ) -> list[dict[str, Any]]:
        """Get events by category."""
        return self.get_events(work_id, event_category=event_category, limit=limit)

    def get_warnings(
        self, work_id: str, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        """Get warning events."""
        return self.get_events(
            work_id, event_category=event_category, event_type="warning", limit=limit
        )

    def get_errors(
        self, work_id: str, event_category: str | None = None, limit: int = 100
    ) -> list[dict[str, Any]]:
        """Get error events."""
        return self.get_events(
            work_id, event_category=event_category, event_type="error", limit=limit
        )

    def get_event_summary_by_category(self, work_id: str) -> dict[str, dict[str, int]]:
        """Get event summary grouped by category and type."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                """
                SELECT event_category, event_type, COUNT(*) as count
                FROM events 
                WHERE work_id = ?
                GROUP BY event_category, event_type
                ORDER BY event_category, event_type
                """,
                (work_id,),
            )

            summary = {}
            for row in cursor.fetchall():
                category = row[0]
                event_type = row[1]
                count = row[2]

                if category not in summary:
                    summary[category] = {}
                summary[category][event_type] = count

            return summary
