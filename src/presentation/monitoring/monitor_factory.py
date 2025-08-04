"""Factory for monitor creation."""

import os
from typing import Any, Dict, Optional, List

from .interfaces import MonitoringInterface, TaskType
from .simple_monitor import SimpleMonitor
from .web_monitor import WebMonitor

class MonitorFactory:
    """Factory for creating monitors based on configuration."""

    @staticmethod
    def create_monitor(config: Dict[str, Any], execution_context: Optional[str] = None) -> Optional[MonitoringInterface]:
        """
        Create a monitor based on configuration and execution context.
        
        For terminal execution, returns SimpleMonitor if monitoring enabled.
        For web execution, returns None (web interface has independent monitoring).

        Args:
            config: Batch configuration
            execution_context: Execution context ("web", "terminal", or None)

        Returns:
            MonitoringInterface or None if monitoring disabled
        """
        # Check if monitoring is enabled
        monitoring_config = config.get("monitoring", {})
        if not monitoring_config.get("enabled", True):  # Default to True for progress monitoring
            return None

        # Detect execution context if not provided
        if execution_context is None:
            # Check for web execution indicators
            is_web = _is_web_execution()
            execution_context = "web" if is_web else "terminal"
            
            # Debug logging
            import logging
            logger = logging.getLogger(__name__)
            logger.debug(f"MonitorFactory: Detected execution context: {execution_context}")

        # Use appropriate monitor based on execution context
        if execution_context == "web":
            import logging
            logger = logging.getLogger(__name__)
            logger.info("MonitorFactory: Web execution detected, using WebMonitor for progress tracking")
            return WebMonitor()

        # For terminal execution, check if hierarchical monitoring is requested
        monitor_type = monitoring_config.get("interface", "simple")
        
        # Check if this is a structured batch that could benefit from hierarchical monitoring
        has_structured_config = (
            config.get("task", {}).get("type") == "execution" and
            "execution" in config and 
            "executions" in config.get("execution", {})
        )
        
        # Use hierarchical monitor for structured batches or when explicitly requested
        if monitor_type == "hierarchical" or (has_structured_config and monitor_type != "simple"):
            import logging
            logger = logging.getLogger(__name__)
            logger.info("MonitorFactory: Using hierarchical monitor for structured batch")
            
            try:
                # Import the hierarchical monitor directly
                import sys
                import os
                current_dir = os.path.dirname(os.path.abspath(__file__))
                workspace_dir = os.path.join(current_dir, '..', '..', '..')
                sys.path.insert(0, workspace_dir)
                
                from hierarchical_simple_monitor import HierarchicalSimpleMonitor
                return HierarchicalSimpleMonitor()
            except ImportError as e:
                logger.warning(f"Failed to import hierarchical monitor: {e}, using legacy SimpleMonitor")
                # Fall through to legacy SimpleMonitor

        # Return legacy SimpleMonitor for compatibility
        try:
            return SimpleMonitor()
        except Exception as e:
            # If SimpleMonitor has issues, create a minimal monitor implementation
            import logging
            logger = logging.getLogger(__name__)
            logger.warning(f"SimpleMonitor instantiation failed: {e}, using basic monitor")
            
            # Create a basic monitor that implements all required methods
            class BasicMonitor(MonitoringInterface):
                def __init__(self):
                    self.task_type = None
                    self.task_name = ""
                    self.is_active = False
                
                def start_task(self, task_type: TaskType, task_name: str, config: Dict[str, Any]) -> None:
                    print(f"📋 Starting {task_type.value}: {task_name}")
                    self.task_type = task_type
                    self.task_name = task_name
                    self.is_active = True
                
                def update_hierarchy_level(self, level: Any, current: int, total: int, name: str = "") -> None:
                    pass
                
                def start_run(self, run: Any) -> None:
                    pass
                
                def finish_run(self, run_id: str, success: bool) -> None:
                    pass
                
                def update_run_progress(self, run_id: str, progress: float, status: str) -> None:
                    pass
                
                def algorithm_callback(self, run_id: str, algorithm_name: str, progress: float, message: str) -> None:
                    pass
                
                def show_error(self, error: str) -> None:
                    print(f"❌ Error: {error}")
                
                def finish_task(self, results: Dict[str, Any]) -> None:
                    print(f"✅ Task completed")
                    self.is_active = False
                
                def initialize_hierarchy(self, hierarchy: Any) -> None:
                    pass
                
                def get_full_status(self) -> Dict[str, Any]:
                    return {"status": "basic"}
                
                def get_active_runs(self) -> List[Any]:
                    return []
                
                def get_recent_callbacks(self) -> List[Any]:
                    return []
                
                def get_hierarchy_status(self) -> Dict[str, Any]:
                    return {}
                
                def get_algorithm_statistics(self) -> Dict[str, Any]:
                    return {}
                
                def get_real_time_data(self) -> Dict[str, Any]:
                    return {}
                
                def get_performance_metrics(self) -> Dict[str, Any]:
                    return {}
                
                def get_completed_runs_history(self) -> List[Any]:
                    return []
                
                def get_detailed_algorithm_data(self, algorithm_name: str) -> Dict[str, Any]:
                    return {}
                
                def update_hierarchy(self, level: Any, level_id: str, progress: float, message: str, data: Dict[str, Any]) -> None:
                    pass
                
                def start_item(self, item_type: str, item_name: str, **kwargs) -> None:
                    pass
                
                def update_item(self, item_id: str, progress: float, message: str = "") -> None:
                    pass
                
                def finish_item(self, item_id: str, success: bool, result: Any = None) -> None:
                    pass
                
                def item_callback(self, item_id: str, message: str, data: Optional[Dict[str, Any]] = None) -> None:
                    pass
                
                def update_execution_data(self, execution_data: Dict[str, Any]) -> None:
                    pass
                
                def get_current_status(self) -> Dict[str, Any]:
                    return {"status": "basic"}
                
                def get_summary(self) -> Dict[str, Any]:
                    return {"summary": "basic"}
                
                def stop(self) -> None:
                    self.is_active = False
                
                def get_session_id(self) -> Optional[str]:
                    return None
                
                def close(self) -> None:
                    self.is_active = False
            
            return BasicMonitor()

    @staticmethod
    def is_monitoring_enabled(config: Dict[str, Any]) -> bool:
        """
        Check if monitoring is enabled in configuration.
        
        For terminal execution, returns True if monitoring enabled for progress display.

        Args:
            config: Batch configuration

        Returns:
            True if monitoring enabled
        """
        monitoring_config = config.get("monitoring", {})
        return monitoring_config.get("enabled", True)


def _is_web_execution() -> bool:
    """
    Detect if execution is happening via web interface.
    
    Returns:
        True if web execution detected
    """
    import inspect
    import threading
    
    # Method 1: Check for FastAPI/uvicorn in the call stack
    frame = inspect.currentframe()
    try:
        while frame:
            filename = frame.f_code.co_filename
            if any(indicator in filename.lower() for indicator in ['fastapi', 'uvicorn', 'starlette', 'async']):
                return True
            if 'batch_execution.py' in filename:
                return True
            frame = frame.f_back
    finally:
        del frame
    
    # Method 2: Check thread name (uvicorn uses specific thread naming)
    thread_name = threading.current_thread().name.lower()
    if any(indicator in thread_name for indicator in ['asyncio', 'uvicorn', 'fastapi']):
        return True
    
    # Method 3: Check environment variables that might indicate web execution
    if os.environ.get('UVICORN_HOST') or os.environ.get('FASTAPI_ENV'):
        return True
    
    return False
