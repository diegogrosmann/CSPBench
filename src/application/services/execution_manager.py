"""
DEPRECATED: ExecutionManager.

This class has been moved to WorkManager as part of the architecture
simplification refactoring. The ExecutionManager was incorporated into
WorkManager to unify work management and pipeline execution.

Migration:
    - Replace ExecutionManager() with get_work_service()
    - The execute() method is available in WorkManager
    - The restart() method has been renamed to restart_execution()

Example:
    Before::
    
        from src.application.services.execution_manager import ExecutionManager
        manager = ExecutionManager()
        work_id = manager.execute(config)
        
    After::
    
        from src.application.services.work_service import get_work_service
        manager = get_work_service()
        work_id = manager.execute(config)

This file will be removed in future versions.
"""

import warnings
from typing import Any, Dict, Optional

from src.application.services.work_service import get_work_service
from src.domain.config import CSPBenchConfig


class ExecutionManager:
    """
    DEPRECATED: Use WorkManager via get_work_service() instead.
    
    This class is a compatibility wrapper that delegates to WorkManager.
    It will be removed in future versions.
    """
    
    def __init__(self, work_service=None):
        """
        Initialize ExecutionManager (deprecated).
        
        Args:
            work_service: Optional work service instance. If None, will use
                the global work service.
        """
        warnings.warn(
            "ExecutionManager is deprecated. Use get_work_service() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        if work_service is None:
            self._work_service = get_work_service()
        else:
            self._work_service = work_service
    
    def execute(self, config: CSPBenchConfig, extra: Optional[Dict[str, Any]] = None) -> str:
        """
        Delegate to WorkManager.execute().
        
        Args:
            config (CSPBenchConfig): Configuration to execute.
            extra (Optional[Dict[str, Any]]): Additional metadata.
            
        Returns:
            str: Work ID for tracking execution.
        """
        return self._work_service.execute(config, extra)
    
    def restart(self, work_id: str) -> bool:
        """
        Delegate to WorkManager.restart_execution().
        
        Args:
            work_id (str): Work ID to restart.
            
        Returns:
            bool: True if restart successful, False otherwise.
        """
        return self._work_service.restart_execution(work_id)