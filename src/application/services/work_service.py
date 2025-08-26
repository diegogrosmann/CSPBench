"""
WorkService - Unified Work Management Service
Provides application-wide persistent work management state.
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager
from typing import Any, Optional

from fastapi import FastAPI

from src.application.work.work_manager import WorkManager
from src.infrastructure.persistence.work_service_persistence import (
    WorkServicePersistence,
)

# Global state
_global_work_service: Optional[WorkManager] = None
_repository: Optional[WorkServicePersistence] = None

logger = logging.getLogger(__name__)


def get_work_service() -> WorkManager:
    """
    Get the global WorkService instance.

    Returns:
        Global WorkService instance

    Raises:
        RuntimeError: If WorkService not initialized
    """
    global _global_work_service
    if _global_work_service is None:
        # Auto-initialize if not already done
        return initialize_work_service()
    return _global_work_service


def initialize_work_service() -> WorkManager:
    """
    Initialize the global WorkService with persistent storage.
    Also pauses any running work items from previous sessions.

    Returns:
        Initialized WorkService instance
    """
    global _global_work_service, _repository

    if _global_work_service is not None:
        logger.info("WorkService already initialized")
        return _global_work_service

    try:
        # Initialize persistent repository
        _repository = WorkServicePersistence()

        # Create WorkManager with persistent repository
        _global_work_service = WorkManager(repository=_repository)

        logger.info("WorkService initialized successfully with persistent storage")
        return _global_work_service

    except Exception as e:
        logger.error(f"Failed to initialize WorkService: {e}")
        raise RuntimeError(f"WorkService initialization failed: {e}") from e


def pause_orphaned_running_work(work_service: WorkManager) -> None:
    """
    Pause any work items that are in 'running' state from previous sessions.
    This prevents orphaned work items that can't be controlled.
    """
    try:
        from src.domain.work import WorkStatus
        from src.domain.status import BaseStatus

        # Get all work items
        all_work = work_service.list()

        running_count = 0
        paused_count = 0

        for work_item in all_work:
            work_id = work_item.get("work_id") or work_item.get("id")
            status = work_item.get("status")

            if status == BaseStatus.RUNNING.value:
                running_count += 1
                try:
                    # Pause the orphaned running work
                    success = work_service.pause(work_id)
                    if success:
                        paused_count += 1
                        logger.info(f"Paused orphaned running work: {work_id}")
                    else:
                        logger.warning(f"Failed to pause orphaned work: {work_id}")
                except Exception as e:
                    logger.error(f"Error pausing orphaned work {work_id}: {e}")

        if running_count > 0:
            logger.info(
                f"Found {running_count} orphaned running work items, paused {paused_count}"
            )
        else:
            logger.debug("No orphaned running work items found")

    except Exception as e:
        logger.error(f"Error during orphaned work cleanup: {e}")
        # Don't fail initialization due to cleanup errors


def cleanup_work_service() -> None:
    """
    Cleanup the global WorkService and close connections.
    """
    global _global_work_service, _repository

    try:
        if _repository:
            _repository.close()
            logger.info("WorkService repository closed")

        _global_work_service = None
        _repository = None
        logger.info("WorkService cleanup completed")

    except Exception as e:
        logger.error(f"Error during WorkService cleanup: {e}")


@asynccontextmanager
async def work_service_lifespan(app: FastAPI):
    """
    FastAPI lifespan context manager for WorkService.

    Initializes WorkService on startup and cleans up on shutdown.
    """
    try:
        # Startup
        logger.info("Starting WorkService...")
        initialize_work_service()
        logger.info("WorkService started successfully")

        yield

    finally:
        # Shutdown
        logger.info("Shutting down WorkService...")
        cleanup_work_service()
        logger.info("WorkService shutdown completed")


def get_work_status(work_id: str):
    """Get work status from the global WorkService."""
    return get_work_service().get_status(work_id)


def get_work_details(work_id: str):
    """Get work details from the global WorkService."""
    return get_work_service().get(work_id)


def list_all_work():
    """List all work items from the global WorkService."""
    return get_work_service().list()


def control_work(work_id: str, action: str) -> bool:
    """Control work execution (pause, resume, cancel, restart)."""
    work_service = get_work_service()

    if action == "pause":
        return work_service.pause(work_id)
    elif action == "resume":
        return work_service.resume(work_id)
    elif action == "cancel":
        return work_service.cancel(work_id)
    elif action == "restart":
        return work_service.restart(work_id)
    else:
        return False
