"""
WebSocket message schemas and data structures.
"""

import json
from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

from src.infrastructure.persistence.work_state.core import (
    ErrorSummary,
    ExecutionDetail,
    ProgressSummary,
)


class MessageType(str, Enum):
    """WebSocket message types."""

    SNAPSHOT = "snapshot"
    UPDATE = "update"
    EVENT = "event"
    ERROR = "error"
    HEARTBEAT = "heartbeat"


class EventType(str, Enum):
    """Event types for real-time notifications."""

    WORK_STATUS_CHANGED = "work_status_changed"
    COMBINATION_CHANGED = "combination_changed"
    EXECUTION_STATUS_CHANGED = "execution_status_changed"
    ERROR_OCCURRED = "error_occurred"
    WARNING_OCCURRED = "warning_occurred"


@dataclass
class WebSocketMessage:
    """Base WebSocket message structure."""

    type: MessageType
    work_id: str
    timestamp: float
    payload: Dict[str, Any]

    def to_json(self) -> str:
        """Convert to JSON string."""
        return json.dumps(asdict(self), ensure_ascii=False)


@dataclass
class ProgressSnapshot:
    """Complete progress snapshot."""

    progress: Dict[str, Any]  # ProgressSummary serialized
    executions: List[Dict[str, Any]]  # ExecutionDetail list serialized
    logs: Dict[str, List[Dict[str, Any]]]  # errors and warnings
    events: Optional[List[Dict[str, Any]]] = None  # events from events table
    executions_full: Optional[List[Dict[str, Any]]] = (
        None  # lista completa por combinação
    )
    combinations: Optional[List[Dict[str, Any]]] = None  # metadados de combinações


@dataclass
class ExecutionChanges:
    """Changes in executions list."""

    updated: List[str] = None  # unit_ids with progress/status changes
    completed: List[str] = None  # unit_ids that completed
    new: List[str] = None  # new unit_ids started
    removed: List[str] = None  # unit_ids no longer in current combination


@dataclass
class ProgressUpdate:
    """Incremental progress update."""

    progress: Optional[Dict[str, Any]] = None  # Only changed fields
    executions_changed: Optional[ExecutionChanges] = None
    logs_appended: Optional[Dict[str, List[Dict[str, Any]]]] = None
    events_appended: Optional[List[Dict[str, Any]]] = None  # New events
    executions: Optional[List[Dict[str, Any]]] = (
        None  # Estado atual completo da combinação corrente (quando mudou)
    )
    executions_combination_id: Optional[int] = (
        None  # Combination id referente à lista completa enviada
    )


@dataclass
class ProgressMessage:
    """Structured progress message for WebSocket."""

    type: MessageType
    work_id: str
    timestamp: float
    snapshot: Optional[ProgressSnapshot] = None
    update: Optional[ProgressUpdate] = None
    event: Optional[Dict[str, Any]] = None
    error: Optional[Dict[str, Any]] = None

    def to_websocket_message(self) -> WebSocketMessage:
        """Convert to WebSocketMessage."""
        payload = {}
        if self.snapshot:
            payload = asdict(self.snapshot)
        elif self.update:
            payload = asdict(self.update)
        elif self.event:
            payload = self.event
        elif self.error:
            payload = self.error

        return WebSocketMessage(
            type=self.type,
            work_id=self.work_id,
            timestamp=self.timestamp,
            payload=payload,
        )


def serialize_progress_summary(progress: ProgressSummary) -> Dict[str, Any]:
    """Serialize ProgressSummary to dict."""
    return {
        "work_id": progress.work_id,
        "tasks": progress.tasks,
        "datasets": progress.datasets,
        "configs": progress.configs,
        "algorithms": progress.algorithms,
        "execution": progress.execution,
        "global_execution": progress.global_execution,
        "global_progress": progress.global_progress,
        "current_combination_details": progress.current_combination_details,
    }


def serialize_execution_detail(execution: ExecutionDetail) -> Dict[str, Any]:
    """Serialize ExecutionDetail to dict."""
    return {
        "unit_id": execution.unit_id,
        "combination_id": execution.combination_id,
        "sequencia": execution.sequencia,
        "status": execution.status,
        "progress": execution.progress,
        "progress_message": execution.progress_message,
        "started_at": execution.started_at,
        "finished_at": execution.finished_at,
        "objective": execution.objective,
        "task_id": execution.task_id,
        "dataset_id": execution.dataset_id,
        "preset_id": execution.preset_id,
        "algorithm_id": execution.algorithm_id,
        "mode": execution.mode,
        "total_sequences": execution.total_sequences,
    }


def serialize_error_summary(error: ErrorSummary) -> Dict[str, Any]:
    """Serialize ErrorSummary to dict."""
    return {
        "unit_id": error.unit_id,
        "error_type": error.error_type,
        "error_message": error.error_message,
        "timestamp": error.timestamp,
    }
