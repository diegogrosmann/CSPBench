"""
WebSocket configuration settings.
"""

from dataclasses import dataclass


@dataclass
class WebSocketConfig:
    """WebSocket configuration."""

    # Update intervals
    update_interval_ms: int = 500  # Base update interval
    heartbeat_interval_s: int = 15  # Heartbeat interval

    # Progress tracking
    progress_precision: int = 3  # Decimal places for progress values

    # Connection limits
    max_connections_per_work: int = 20
    max_message_queue_size: int = 100

    # Timeouts
    connection_timeout_s: int = 60
    grace_period_after_completion_s: int = 30

    # Logging
    debug_mode: bool = False


# Global configuration instance
websocket_config = WebSocketConfig()
