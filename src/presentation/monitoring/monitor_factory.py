"""Factory for monitor creation."""

from typing import Any, Dict, Optional

from .interfaces import MonitoringInterface
from .simple_monitor import SimpleMonitor

# from .tui_monitor import TUIMonitor  # Temporariamente desabilitado


class MonitorFactory:
    """Factory for creating monitors based on configuration."""

    @staticmethod
    def create_monitor(config: Dict[str, Any]) -> Optional[MonitoringInterface]:
        """
        Create a monitor based on configuration.

        Args:
            config: Batch configuration

        Returns:
            MonitoringInterface or None if monitoring disabled
        """
        # Check if monitoring is enabled
        monitoring_config = config.get("monitoring", {})
        if not monitoring_config.get("enabled", True):
            return None

        # Determine interface type
        interface_type = monitoring_config.get("interface", "simple")

        if interface_type == "tui":
            # return TUIMonitor()  # Temporarily disabled
            return SimpleMonitor()  # Use SimpleMonitor as fallback
        elif interface_type == "simple":
            return SimpleMonitor()
        else:
            # Default to simple interface
            return SimpleMonitor()

    @staticmethod
    def is_monitoring_enabled(config: Dict[str, Any]) -> bool:
        """
        Check if monitoring is enabled in configuration.

        Args:
            config: Batch configuration

        Returns:
            True if monitoring enabled
        """
        monitoring_config = config.get("monitoring", {})
        return monitoring_config.get("enabled", True)
