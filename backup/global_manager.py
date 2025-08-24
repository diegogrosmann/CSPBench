"""
Global WorkManager with FastAPI Lifespan Integration
Provides application-wide persistent work management state.
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager
from typing import Any, Optional

from fastapi import FastAPI

from application.work.work_manager import WorkManager
from infrastructure.persistence.persistent_repository import PersistentWorkRepository

# Global state
_global_work_manager: Optional[WorkManager] = None
_repository: Optional[PersistentWorkRepository] = None

logger = logging.getLogger(__name__)


def get_global_work_manager() -> WorkManager:
    """
    Get the global WorkManager instance.
    
    Returns:
        Global WorkManager instance
        
    Raises:
        RuntimeError: If WorkManager not initialized
    """
    global _global_work_manager
    if _global_work_manager is None:
        raise RuntimeError(
            "Global WorkManager not initialized. "
            "Make sure FastAPI lifespan is properly configured."
        )
    return _global_work_manager


def initialize_global_work_manager() -> WorkManager:
    """
    Initialize the global WorkManager with persistent storage.
    
    Returns:
        Initialized WorkManager instance
    """
    global _global_work_manager, _repository
    
    if _global_work_manager is not None:
        logger.warning("Global WorkManager already initialized")
        return _global_work_manager
    
    try:
        # Initialize persistent repository
        _repository = PersistentWorkRepository()
        logger.info(f"Persistent repository initialized: {_repository.db_path}")
        
        # Initialize WorkManager with persistent repository
        _global_work_manager = WorkManager(repository=_repository)
        logger.info("Global WorkManager initialized successfully")
        
        # Log current state
        item_count = _repository.count()
        logger.info(f"WorkManager loaded with {item_count} existing work items")
        
        return _global_work_manager
        
    except Exception as e:
        logger.error(f"Failed to initialize global WorkManager: {e}")
        raise


def cleanup_global_work_manager() -> None:
    """
    Cleanup the global WorkManager and close connections.
    """
    global _global_work_manager, _repository
    
    try:
        if _repository:
            # Optional: cleanup old items on shutdown
            removed_count = _repository.cleanup_old_items(older_than_days=7)
            if removed_count > 0:
                logger.info(f"Cleaned up {removed_count} old work items")
            
            # Close database connection
            _repository.close()
            logger.info("Repository connection closed")
            
        _global_work_manager = None
        _repository = None
        logger.info("Global WorkManager cleanup completed")
        
    except Exception as e:
        logger.error(f"Error during WorkManager cleanup: {e}")


@asynccontextmanager
async def work_manager_lifespan(app: FastAPI):
    """
    FastAPI lifespan context manager for WorkManager.
    
    Usage:
        app = FastAPI(lifespan=work_manager_lifespan)
    """
    # Startup
    logger.info("Starting WorkManager lifespan...")
    try:
        initialize_global_work_manager()
        logger.info("WorkManager lifespan startup completed")
        yield
    finally:
        # Shutdown
        logger.info("Shutting down WorkManager lifespan...")
        cleanup_global_work_manager()
        logger.info("WorkManager lifespan shutdown completed")


# Convenience functions for common operations
def submit_work(config: Any, extra: dict[str, Any] | None = None) -> tuple[str, dict[str, Any]]:
    """Submit work to global WorkManager."""
    return get_global_work_manager().submit(config=config, extra=extra)


def get_work_status(work_id: str) -> Optional[str]:
    """Get work status from global WorkManager."""
    return get_global_work_manager().get_status(work_id)


def get_work_details(work_id: str) -> dict[str, Any] | None:
    """Get work details from global WorkManager."""
    return get_global_work_manager().get(work_id)


def list_all_work() -> list[dict[str, Any]]:
    """List all work items from global WorkManager."""
    return get_global_work_manager().list()


def control_work(work_id: str, action: str) -> bool:
    """
    Control work execution (pause, resume, cancel, etc.).
    
    Args:
        work_id: Work item ID
        action: Control action ('pause', 'resume', 'cancel', 'restart')
        
    Returns:
        True if action was successful
    """
    manager = get_global_work_manager()
    
    if action == "pause":
        return manager.pause(work_id)
    elif action == "resume":
        return manager.resume(work_id)
    elif action == "cancel":
        return manager.cancel(work_id)
    elif action == "restart":
        return manager.restart(work_id)
    else:
        logger.warning(f"Unknown work control action: {action}")
        return False
