"""Factory for creating monitoring components."""

from typing import Optional, Any
import logging

from .progress_broker import ProgressBroker
from .monitoring_service import MonitoringService


class MonitoringFactory:
    """
    Factory for creating monitoring system components.

    This factory creates the appropriate display adapters and monitoring service
    based on the execution context (terminal vs web).
    """

    @staticmethod
    def create_monitoring_system(
        session_id: Optional[str] = None,
        web_session_manager=None,
        enable_terminal_display: bool = True,
        terminal_verbose: bool = True,
    ) -> tuple[MonitoringService, ProgressBroker]:
        """
        Create a complete monitoring system.

        Args:
            session_id: Optional session ID for web sessions
            web_session_manager: Web session manager for web display
            enable_terminal_display: Whether to enable terminal display
            terminal_verbose: Whether terminal display should be verbose

        Returns:
            Tuple of (MonitoringService, ProgressBroker)
        """
        logger = logging.getLogger(__name__)

        # Create progress broker
        broker = ProgressBroker()

        # Create monitoring service
        monitoring_service = MonitoringService(broker, session_id)

        # Create and register display adapters
        displays_created = []

        # Terminal display (if enabled and not in web-only mode)
        if enable_terminal_display and not MonitoringFactory._is_web_only_mode():
            from src.presentation.display.live_terminal_display import (
                LiveTerminalDisplay,
            )

            terminal_display = LiveTerminalDisplay(verbose=terminal_verbose)
            broker.subscribe_all(terminal_display.handle_event)
            displays_created.append("Live Terminal")
            logger.info("Live terminal display registered")

        # Web display (if web session manager is provided)
        if web_session_manager and session_id:
            from src.presentation.display.web_display import WebDisplay

            web_display = WebDisplay(web_session_manager, session_id)
            broker.subscribe_all(web_display.handle_event)
            displays_created.append("Web")
            logger.info(f"Web display registered for session {session_id}")

        logger.info(
            f"Monitoring system created with displays: {', '.join(displays_created)}"
        )

        return monitoring_service, broker

    @staticmethod
    def _is_web_only_mode() -> bool:
        """
        Check if we're running in web-only mode.

        Returns:
            True if in web-only mode (no terminal display needed)
        """
        # This could check environment variables, execution context, etc.
        # For now, we'll use a simple heuristic
        import os

        return os.getenv("WEB_ONLY_MODE", "false").lower() == "true"

    @staticmethod
    def create_terminal_only_system(
        verbose: bool = True,
    ) -> tuple[MonitoringService, ProgressBroker]:
        """
        Create a terminal-only monitoring system.

        Args:
            verbose: Whether terminal display should be verbose

        Returns:
            Tuple of (MonitoringService, ProgressBroker)
        """
        return MonitoringFactory.create_monitoring_system(
            enable_terminal_display=True, terminal_verbose=verbose
        )

    @staticmethod
    def create_web_only_system(
        session_id: str, web_session_manager
    ) -> tuple[MonitoringService, ProgressBroker]:
        """
        Create a web-only monitoring system.

        Args:
            session_id: Session ID for web session
            web_session_manager: Web session manager

        Returns:
            Tuple of (MonitoringService, ProgressBroker)
        """
        return MonitoringFactory.create_monitoring_system(
            session_id=session_id,
            web_session_manager=web_session_manager,
            enable_terminal_display=False,
        )
