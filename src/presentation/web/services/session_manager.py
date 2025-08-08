"""
Execution Session Manager

This module manages execution sessions for tracking the status of
batch executions and other long-running tasks.
"""

import threading
import time
import uuid
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional


@dataclass
class ExecutionSession:
    """Represents an execution session."""

    session_id: str
    session_type: str
    status: str  # pending, running, completed, failed, cancelled
    created_at: float
    updated_at: float
    completed_at: Optional[float] = None
    config: Optional[Dict[str, Any]] = None
    message: Optional[str] = None
    error: Optional[str] = None
    progress: Optional[Dict[str, Any]] = None
    progress_state: Optional[Dict[str, Any]] = None  # For WebDisplay compatibility
    logs: Optional[List[Dict[str, Any]]] = None
    current_execution: Optional[Dict[str, Any]] = None
    last_updated: Optional[str] = None  # For WebDisplay compatibility

    def __post_init__(self):
        if self.logs is None:
            self.logs = []


class ExecutionSessionManager:
    """Manages execution sessions."""

    def __init__(self):
        self._sessions: Dict[str, ExecutionSession] = {}
        self._lock = threading.Lock()

    def create_session(
        self,
        session_type: str,
        config: Optional[Dict[str, Any]] = None,
        status: str = "pending",
    ) -> str:
        """Create a new execution session."""
        with self._lock:
            session_id = str(uuid.uuid4())
            current_time = time.time()

            session = ExecutionSession(
                session_id=session_id,
                session_type=session_type,
                status=status,
                created_at=current_time,
                updated_at=current_time,
                config=config or {},
                message=f"Session created: {session_type}",
            )

            self._sessions[session_id] = session
            return session_id

    def get_session(self, session_id: str) -> Optional[Dict[str, Any]]:
        """Get session details."""
        with self._lock:
            session = self._sessions.get(session_id)
            if session:
                return asdict(session)
            return None

    def update_session(self, session_id: str, updates: Dict[str, Any]) -> bool:
        """Update session with new data."""
        with self._lock:
            session = self._sessions.get(session_id)
            if not session:
                return False

            # Update allowed fields
            allowed_fields = {
                "status",
                "message",
                "error",
                "progress",
                "progress_state",  # Added for WebDisplay compatibility
                "completed_at",
                "config",
                "logs",
                "current_execution",
                "last_updated",  # Added for WebDisplay compatibility
            }

            for key, value in updates.items():
                if key in allowed_fields:
                    setattr(session, key, value)

            session.updated_at = time.time()
            return True

    def add_log(
        self, session_id: str, level: str, message: str, source: Optional[str] = None
    ) -> bool:
        """Add a log entry to a session."""
        with self._lock:
            session = self._sessions.get(session_id)
            if not session:
                return False

            if session.logs is None:
                session.logs = []

            # Gerar ID único para o log
            log_id = len(session.logs) + 1

            # Gerar timestamp legível
            timestamp = time.strftime("%H:%M:%S", time.localtime())

            log_entry = {
                "id": log_id,
                "timestamp": timestamp,
                "level": level.upper(),
                "message": message,
            }

            if source:
                log_entry["source"] = source

            session.logs.append(log_entry)
            session.updated_at = time.time()
            return True

    def get_session_logs(self, session_id: str) -> List[Dict[str, Any]]:
        """Get logs for a specific session."""
        with self._lock:
            session = self._sessions.get(session_id)
            if session and session.logs:
                return session.logs.copy()
            return []

    def clear_session_logs(self, session_id: str) -> bool:
        """Clear all logs for a session."""
        with self._lock:
            session = self._sessions.get(session_id)
            if not session:
                return False

            session.logs = []
            session.updated_at = time.time()
            return True

    def complete_session(
        self, session_id: str, status: str = "completed", message: Optional[str] = None
    ) -> bool:
        """Mark session as completed."""
        updates = {"status": status, "completed_at": time.time()}
        if message:
            updates["message"] = message

        return self.update_session(session_id, updates)

    def fail_session(
        self, session_id: str, error: str, message: Optional[str] = None
    ) -> bool:
        """Mark session as failed."""
        updates = {"status": "failed", "error": error, "completed_at": time.time()}
        if message:
            updates["message"] = message
        else:
            updates["message"] = f"Session failed: {error}"

        return self.update_session(session_id, updates)

    def cancel_session(self, session_id: str) -> bool:
        """Cancel a running session."""
        with self._lock:
            session = self._sessions.get(session_id)
            if not session:
                return False

            if session.status not in ["pending", "running"]:
                return False

            session.status = "cancelled"
            session.completed_at = time.time()
            session.updated_at = time.time()
            session.message = "Session cancelled by user"

            return True

    def list_sessions(
        self,
        session_type: Optional[str] = None,
        status: Optional[str] = None,
        limit: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """List sessions with optional filtering."""
        with self._lock:
            sessions = list(self._sessions.values())

            # Apply filters
            if session_type:
                sessions = [s for s in sessions if s.session_type == session_type]

            if status:
                sessions = [s for s in sessions if s.status == status]

            # Sort by creation time (newest first)
            sessions.sort(key=lambda x: x.created_at, reverse=True)

            # Apply limit
            if limit:
                sessions = sessions[:limit]

            return [asdict(session) for session in sessions]

    def cleanup_old_sessions(self, max_age_hours: int = 24) -> int:
        """Remove old sessions to prevent memory buildup."""
        with self._lock:
            current_time = time.time()
            max_age_seconds = max_age_hours * 3600

            old_session_ids = []
            for session_id, session in self._sessions.items():
                # Only cleanup completed, failed, or cancelled sessions
                if session.status in ["completed", "failed", "cancelled"]:
                    age = current_time - session.created_at
                    if age > max_age_seconds:
                        old_session_ids.append(session_id)

            # Remove old sessions
            for session_id in old_session_ids:
                del self._sessions[session_id]

            return len(old_session_ids)

    def get_session_count(self) -> Dict[str, int]:
        """Get count of sessions by status."""
        with self._lock:
            counts = {}
            for session in self._sessions.values():
                status = session.status
                counts[status] = counts.get(status, 0) + 1

            return counts

    def get_running_sessions(self) -> List[Dict[str, Any]]:
        """Get all currently running sessions."""
        return self.list_sessions(status="running")

    def _get_current_time(self) -> float:
        """Get current timestamp."""
        return time.time()


# Global session manager instance
_session_manager = None


def get_session_manager() -> ExecutionSessionManager:
    """Get the global session manager instance."""
    global _session_manager
    if _session_manager is None:
        _session_manager = ExecutionSessionManager()
    return _session_manager
