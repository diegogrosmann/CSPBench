"""
Utility functions for web services.
"""

import asyncio
import logging
import tempfile
from pathlib import Path
from typing import Any, Dict

from .session_manager import ExecutionSessionManager, get_session_manager

logger = logging.getLogger(__name__)


class SessionManager:
    """Manages temporary execution sessions."""

    def __init__(self):
        self.active_sessions = {}

    def create_session(self, session_id: str) -> Path:
        """Create a new session directory."""
        session_dir = Path(tempfile.mkdtemp(prefix=f"session_{session_id}_"))
        self.active_sessions[session_id] = {
            "directory": session_dir,
            "created_at": asyncio.get_event_loop().time(),
        }
        return session_dir

    def get_session_dir(self, session_id: str) -> Path:
        """Get session directory."""
        session_info = self.active_sessions.get(session_id)
        if session_info:
            return session_info["directory"]
        return None

    def cleanup_session(self, session_id: str):
        """Clean up session directory."""
        try:
            session_info = self.active_sessions.pop(session_id, None)
            if session_info:
                import shutil

                if session_info["directory"].exists():
                    shutil.rmtree(session_info["directory"])
                logger.info(f"Cleaned up session {session_id}")
        except Exception as e:
            logger.warning(f"Error cleaning up session {session_id}: {e}")

    def cleanup_old_sessions(self, max_age_hours: int = 24):
        """Clean up old sessions."""
        try:
            current_time = asyncio.get_event_loop().time()
            max_age_seconds = max_age_hours * 3600

            old_sessions = [
                session_id
                for session_id, info in self.active_sessions.items()
                if current_time - info["created_at"] > max_age_seconds
            ]

            for session_id in old_sessions:
                self.cleanup_session(session_id)

            if old_sessions:
                logger.info(f"Cleaned up {len(old_sessions)} old sessions")

        except Exception as e:
            logger.error(f"Error cleaning up old sessions: {e}")


class BackgroundTasks:
    """Manages background task execution."""

    def __init__(self):
        self.running_tasks = {}

    def start_task(self, task_id: str, coro):
        """Start a background task."""
        task = asyncio.create_task(coro)
        self.running_tasks[task_id] = task

        # Clean up when done
        def cleanup_task(task):
            self.running_tasks.pop(task_id, None)

        task.add_done_callback(cleanup_task)
        return task

    def get_task_status(self, task_id: str) -> str:
        """Get task status."""
        task = self.running_tasks.get(task_id)
        if not task:
            return "not_found"
        elif task.done():
            if task.exception():
                return "failed"
            return "completed"
        else:
            return "running"

    def cancel_task(self, task_id: str) -> bool:
        """Cancel a running task."""
        task = self.running_tasks.get(task_id)
        if task and not task.done():
            task.cancel()
            return True
        return False


# Global instances
session_manager = SessionManager()
background_tasks = BackgroundTasks()


async def cleanup_worker():
    """Background worker for cleanup tasks."""
    while True:
        try:
            await asyncio.sleep(3600)  # Run every hour
            session_manager.cleanup_old_sessions()
        except Exception as e:
            logger.error(f"Error in cleanup worker: {e}")


def start_background_workers():
    """Start background workers."""
    try:
        asyncio.create_task(cleanup_worker())
        logger.info("Background workers started")
    except Exception as e:
        logger.error(f"Error starting background workers: {e}")
