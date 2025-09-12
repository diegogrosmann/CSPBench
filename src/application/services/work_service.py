"""
WorkService - Unified Work Management Service.

Provides application-wide persistent work management state with global
singleton pattern and FastAPI lifecycle integration.

This module implements a global work service that manages the lifecycle
of work items across the entire application, ensuring persistence and
proper cleanup of resources.

Features:
- Global singleton work service
- Persistent storage integration
- Orphaned work cleanup on startup
- FastAPI lifecycle management
- Thread-safe operations
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager
from typing import Optional

from fastapi import FastAPI

from src.application.work.work_manager import WorkManager
from src.infrastructure.persistence.work_state.core import WorkPersistence

# Global state
_global_work_service: Optional[WorkManager] = None

logger = logging.getLogger(__name__)

def get_work_service() -> WorkManager:
    """
    Get the global WorkService instance.

    Returns:
        WorkManager: Global WorkService instance

    Raises:
        RuntimeError: If WorkService not initialized
        
    Note:
        Auto-initializes if not already done for convenience.
    """
    global _global_work_service
    if _global_work_service is None:
        # Auto-initialize if not already done
        return initialize_work_service()
    return _global_work_service


def initialize_work_service() -> WorkManager:
    """
    Initialize the global WorkService with persistent storage.
    
    Also pauses any running work items from previous sessions to prevent
    orphaned work items that cannot be controlled.

    Returns:
        WorkManager: Initialized WorkService instance
        
    Raises:
        RuntimeError: If WorkService initialization fails
    """
    global _global_work_service

    if _global_work_service is not None:
        logger.info("WorkService already initialized")
        return _global_work_service

    try:
        _global_work_service = WorkManager()
        
        # Perform initialization cleanup and recovery
        init_summary = _global_work_service.initialize()
        logger.info(f"WorkManager initialization completed: {init_summary}")

        logger.info("WorkService initialized successfully with persistent storage")
        return _global_work_service

    except Exception as e:
        logger.error(f"Failed to initialize WorkService: {e}")
        raise RuntimeError(f"WorkService initialization failed: {e}") from e


def cleanup_work_service() -> None:
    """
    Cleanup the global WorkService and close connections.
    
    Properly closes repository connections and resets global state.
    Should be called during application shutdown.
    """
    global _global_work_service

    try:
        if _global_work_service:
            _global_work_service.close()
            logger.info("WorkService closed")

        _global_work_service = None
        logger.info("WorkService cleanup completed")

    except Exception as e:
        logger.error(f"Error during WorkService cleanup: {e}")


@asynccontextmanager
async def work_service_lifespan(app: FastAPI):
    """
    FastAPI lifespan context manager for WorkService.

    Initializes WorkService on startup and cleans up on shutdown.
    Provides proper resource management for FastAPI applications.
    
    Args:
        app: FastAPI application instance
        
    Yields:
        None: Application runs between startup and shutdown
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
    """
    Get work status from the global WorkService.
    
    Args:
        work_id: Unique work identifier
        
    Returns:
        str: Current work status, None if not found
    """
    return get_work_service().get_status(work_id)


def get_work_details(work_id: str):
    """
    Get work details from the global WorkService.
    
    Args:
        work_id: Unique work identifier
        
    Returns:
        WorkItem: Work item instance, None if not found
    """
    return get_work_service().get(work_id)


def list_all_work():
    """
    List all work items from the global WorkService.
    
    Returns:
        list: List of all work items as WorkItem instances
    """
    return get_work_service().list()


def control_work(work_id: str, action: str) -> bool:
    """
    Control work execution (pause, cancel, restart).
    
    Args:
        work_id: Unique work identifier
        action: Control action ('pause', 'cancel', 'restart')
        
    Returns:
        bool: True if action successful, False otherwise
    """
    work_service = get_work_service()

    if action == "pause":
        return work_service.pause(work_id)
    elif action == "cancel":
        return work_service.cancel(work_id)
    elif action == "restart":
        return work_service.restart(work_id)
    else:
        return False
