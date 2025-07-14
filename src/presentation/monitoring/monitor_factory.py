"""Factory para criação de monitores."""

from typing import Any, Dict, Optional

from .interfaces import MonitoringInterface, TaskType
from .simple_monitor import SimpleMonitor

# from .tui_monitor import TUIMonitor  # Temporariamente desabilitado


class MonitorFactory:
    """Factory para criação de monitores baseado na configuração."""

    @staticmethod
    def create_monitor(config: Dict[str, Any]) -> Optional[MonitoringInterface]:
        """
        Cria um monitor baseado na configuração.

        Args:
            config: Configuração do batch

        Returns:
            MonitoringInterface ou None se monitoramento desabilitado
        """
        # Verifica se monitoramento está habilitado
        monitoring_config = config.get("monitoring", {})
        if not monitoring_config.get("enabled", True):
            return None

        # Determina tipo de interface
        interface_type = monitoring_config.get("interface", "simple")

        if interface_type == "tui":
            # return TUIMonitor()  # Temporariamente desabilitado
            return SimpleMonitor()  # Usar SimpleMonitor como fallback
        elif interface_type == "simple":
            return SimpleMonitor()
        else:
            # Default para interface simples
            return SimpleMonitor()

    @staticmethod
    def is_monitoring_enabled(config: Dict[str, Any]) -> bool:
        """
        Verifica se monitoramento está habilitado na configuração.

        Args:
            config: Configuração do batch

        Returns:
            True se monitoramento habilitado
        """
        monitoring_config = config.get("monitoring", {})
        return monitoring_config.get("enabled", True)
