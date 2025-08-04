"""Web monitor for silent progress tracking."""

from typing import Any, Dict, Optional

from .interfaces import (
    ExecutionLevel,
    HierarchicalContext,
    MonitoringInterface,
    TaskType,
)


class WebMonitor(MonitoringInterface):
    """Silent monitor that tracks progress without terminal output."""

    def __init__(self):
        """Initialize web monitor."""
        self.is_active = False
        self.current_task = None
        self.current_progress = 0.0

    def start_task(
        self,
        task_type: TaskType,
        task_name: str,
        config: Dict[str, Any],
        total_items: int = 0,
    ) -> None:
        """Start monitoring a task silently."""
        self.is_active = True
        self.current_task = {
            "type": task_type,
            "name": task_name,
            "config": config,
            "total_items": total_items,
        }

    def finish_task(
        self,
        success: bool = True,
        final_results: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Finish task monitoring silently."""
        self.is_active = False
        self.current_progress = 100.0 if success else 0.0

    def update_hierarchy(
        self,
        level: ExecutionLevel,
        level_id: str,
        progress: float,
        message: str = "",
        data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Update hierarchical progress silently."""
        # Silent update - no terminal output
        self.current_progress = progress

    def update_item(
        self,
        item_id: str,
        progress: float,
        message: str = "",
        context: Optional[HierarchicalContext] = None,
    ) -> None:
        """Update individual item progress silently."""
        # Silent update - no terminal output
        pass

    def start_item(
        self,
        item_id: str,
        item_type: str = "repetition",
        context: Optional[HierarchicalContext] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Start item monitoring silently."""
        # Silent start - no terminal output
        pass

    def finish_item(
        self,
        item_id: str,
        success: bool = True,
        result: Optional[Dict[str, Any]] = None,
        error: Optional[str] = None,
    ) -> None:
        """Finish item monitoring silently."""
        # Silent finish - no terminal output
        pass

    def show_error(self, error: str) -> None:
        """Handle error silently."""
        # Silent error handling - no terminal output
        pass

    def close(self) -> None:
        """Close monitoring silently."""
        self.is_active = False

    def algorithm_callback(
        self,
        algorithm_name: str,
        progress: float,
        message: str,
        item_id: Optional[str] = None,
    ) -> None:
        """Algorithm callback silently."""
        # Silent callback - no terminal output
        pass

    def stop(self) -> None:
        """Stop monitoring."""
        self.is_active = False

    def get_summary(self) -> Dict[str, Any]:
        """Get monitoring summary."""
        return {
            "task": self.current_task,
            "progress": self.current_progress,
            "is_active": self.is_active,
        }
