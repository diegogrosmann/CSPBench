"""
Logging Configuration for CSPBench

Centralized logging configurations without external console dependencies.
"""

import logging
import logging.handlers
from pathlib import Path
from typing import Optional


def setup_basic_logging(
    level: str = "INFO",
    log_dir: str = "outputs/logs",
    base_name: str = "cspbench",
    max_bytes: int = 10 * 1024 * 1024,  # 10MB
    backup_count: int = 5,
    log_file_path: Optional[str] = None,  # Full file path if specified
) -> None:
    """
    Configure basic system logging.

    Args:
        level: Log level (DEBUG, INFO, WARNING, ERROR)
        log_dir: Directory for log files
        base_name: Base name for log files
        max_bytes: Maximum file size before rotation
        backup_count: Number of backups to maintain
        log_file_path: Full file path (overrides log_dir + base_name)
    """
    # Determine file path
    if log_file_path:
        log_file = Path(log_file_path)
        log_file.parent.mkdir(parents=True, exist_ok=True)
    else:
        # Create directory if it doesn't exist
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        log_file = log_path / f"{base_name}.log"

    # Configure formatting
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Configure file handler with rotation
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
    )
    file_handler.setFormatter(formatter)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Add new handler
    root_logger.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    """
    Returns configured logger for specific module.

    Args:
        name: Module/logger name

    Returns:
        logging.Logger: Configured logger
    """
    return logging.getLogger(name)


class LoggerConfig:
    """Centralized logger configuration."""

    _initialized = False
    _log_level = "INFO"

    @classmethod
    def initialize(cls, level: str = "INFO", **kwargs) -> None:
        """Initialize logging configuration."""
        if not cls._initialized:
            setup_basic_logging(level=level, **kwargs)
            cls._log_level = level
            cls._initialized = True

    @classmethod
    def get_level(cls) -> str:
        """Return current logging level."""
        return cls._log_level

    @classmethod
    def set_level(cls, level: str) -> None:
        """Change logging level."""
        cls._log_level = level
        logging.getLogger().setLevel(getattr(logging, level.upper(), logging.INFO))

    @classmethod
    def is_initialized(cls) -> bool:
        """Check if logging was initialized."""
        return cls._initialized
