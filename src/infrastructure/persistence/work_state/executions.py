"""Executions mixin for work state persistence."""

import json
import time
from typing import Any, Optional

from .utils import validate_status, validate_required_field


class ExecutionsMixin:
    """Mixin providing execution management functionality."""

    def submit_execution(
        self,
        *,
        unit_id: str,
        combination_id: int,
        sequencia: int | None = None,
    ) -> None:
        """Record execution start. Always creates as 'queued'."""

        self._execute(
            """
            INSERT OR REPLACE INTO executions(
              unit_id, combination_id, sequencia, status, started_at, 
              finished_at, params_json, result_json, objective
            ) VALUES(?,?,?,?,NULL,NULL,NULL,NULL,NULL)
            """,
            (
                unit_id,
                combination_id,
                sequencia,
                "queued",
            ),
        )

    def update_execution_status(
        self,
        unit_id: str,
        status: str,
        result: dict[str, Any] | None = None,
        objective: float | None = None,
        params: dict[str, Any] | None = None,
    ) -> None:
        """Update execution status with conditional field updates."""
        validate_required_field(unit_id, "unit_id")
        validate_status(status)

        timestamp = time.time()

        if status == "running":
            self._execute(
                """
                UPDATE executions 
                SET status=?, started_at=? 
                WHERE unit_id=?
                """,
                (status, timestamp, unit_id),
            )
        elif status in ("completed", "failed"):
            res_js = json.dumps(result or {}, ensure_ascii=False)
            params_js = json.dumps(params or {}, ensure_ascii=False)

            self._execute(
                """
                UPDATE executions 
                SET status=?, finished_at=?, result_json=?, objective=?, params_json=?
                WHERE unit_id=?
                """,
                (status, timestamp, res_js, objective, params_js, unit_id),
            )
        else:
            self._execute(
                """
                UPDATE executions 
                SET status=? 
                WHERE unit_id=?
                """,
                (status, unit_id),
            )

    def get_executions(
        self,
        *,
        unit_id: Optional[str] = None,
        status: Optional[str] = None,
        combination_id: Optional[int] = None,
        sequencia: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """Get executions with optional filters."""
        query = "SELECT * FROM executions WHERE 1=1"
        params = []

        if unit_id is not None:
            query += " AND unit_id = ?"
            params.append(unit_id)

        if status is not None:
            validate_status(status)
            query += " AND status = ?"
            params.append(status)

        if combination_id is not None:
            query += " AND combination_id = ?"
            params.append(combination_id)

        if sequencia is not None:
            query += " AND sequencia = ?"
            params.append(sequencia)

        query += " ORDER BY started_at DESC"

        rows = self._fetch_all(query, params)

        executions = []
        for row in rows:
            execution = dict(row)
            # Deserialize JSON fields
            if execution.get("params_json"):
                execution["params"] = json.loads(execution["params_json"])
            else:
                execution["params"] = {}

            if execution.get("result_json"):
                execution["result"] = json.loads(execution["result_json"])
            else:
                execution["result"] = {}

            # Remove raw JSON fields
            execution.pop("params_json", None)
            execution.pop("result_json", None)

            executions.append(execution)

        return executions

    def add_execution_progress(
        self,
        unit_id: str,
        progress: float,
        message: str | None = None,
    ) -> None:
        """Add progress entry for an execution."""
        validate_required_field(unit_id, "unit_id")

        if not (0.0 <= progress <= 1.0):
            raise ValueError("Progress must be between 0.0 and 1.0")

        # Get execution_id from unit_id
        execution_row = self._fetch_one(
            "SELECT id FROM executions WHERE unit_id = ?", (unit_id,)
        )

        if not execution_row:
            raise ValueError(f"No execution found for unit_id: {unit_id}")

        execution_id = execution_row["id"]

        self._execute(
            """
            INSERT INTO execution_progress(execution_id, progress, message, timestamp)
            VALUES(?, ?, ?, ?)
            """,
            (execution_id, progress, message, time.time()),
        )

    def get_execution_progress(
        self,
        unit_id: str,
        limit: int | None = None,
    ) -> list[dict[str, Any]]:
        """Get progress entries for an execution."""
        validate_required_field(unit_id, "unit_id")

        query = """
        SELECT ep.progress, ep.message, ep.timestamp
        FROM execution_progress ep
        JOIN executions e ON ep.execution_id = e.id
        WHERE e.unit_id = ?
        ORDER BY ep.timestamp DESC
        """

        params = [unit_id]

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = self._fetch_all(query, params)
        return [dict(row) for row in rows]
