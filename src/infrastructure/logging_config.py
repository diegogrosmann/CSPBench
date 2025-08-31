"""
Logging Configuration for CSPBench

Centralized logging configurations without external console dependencies.
"""

import logging
import logging.handlers
import os
from pathlib import Path
from typing import Optional


def setup_basic_logging(
    level: Optional[str] = None,
    log_dir: Optional[str] = None,
    base_name: Optional[str] = None,
    max_bytes: Optional[int] = None,
    backup_count: Optional[int] = None,
    log_file_path: Optional[str] = None,  # Full file path if specified
    log_to_stdout: Optional[bool] = None,
) -> None:
    """
    Configure basic system logging using environment variables as defaults.

    Args:
        level: Log level (DEBUG, INFO, WARNING, ERROR) - defaults to LOG_LEVEL env var
        log_dir: Directory for log files - defaults to LOG_DIRECTORY env var
        base_name: Base name for log files - defaults to LOG_BASE_NAME env var
        max_bytes: Maximum file size before rotation - defaults to LOG_MAX_BYTES env var
        backup_count: Number of backups to maintain - defaults to LOG_BACKUP_COUNT env var
        log_file_path: Full file path (overrides log_dir + base_name)
    """
    # Use environment variables as defaults
    level = level or os.getenv("LOG_LEVEL", "INFO")
    log_dir = log_dir or os.getenv("LOG_DIRECTORY", "outputs/logs")
    base_name = base_name or os.getenv("LOG_BASE_NAME", "cspbench")
    max_bytes = max_bytes or int(os.getenv("LOG_MAX_BYTES", "10485760"))  # 10MB
    backup_count = backup_count or int(os.getenv("LOG_BACKUP_COUNT", "5"))
    if log_to_stdout is None:
        env_flag = os.getenv("LOG_TO_STDOUT", "false").lower()
        log_to_stdout = env_flag in {"1", "true", "yes", "on"}
    # Determine file path (only if not logging to stdout)
    if not log_to_stdout:
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

    handlers = []
    if log_to_stdout:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        handlers.append(stream_handler)
    else:
        # Configure file handler with rotation
        file_handler = logging.handlers.RotatingFileHandler(
            log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
        )
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Add new handlers
    for h in handlers:
        root_logger.addHandler(h)


def setup_logging_from_env() -> None:
    """
    Configure logging using only environment variables.

    This is a convenience function that sets up logging entirely from .env configuration.
    """
    setup_basic_logging()


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
    _log_level = None

    @classmethod
    def initialize(cls, level: Optional[str] = None, **kwargs) -> None:
        """Initialize logging configuration."""
        if not cls._initialized:
            # Use environment variable if level not provided
            level = level or os.getenv("LOG_LEVEL", "INFO")
            setup_basic_logging(level=level, **kwargs)
            cls._log_level = level
            cls._initialized = True

    @classmethod
    def get_level(cls) -> str:
        """Return current logging level."""
        return cls._log_level or os.getenv("LOG_LEVEL", "INFO")

    @classmethod
    def set_level(cls, level: str) -> None:
        """Change logging level."""
        cls._log_level = level
        logging.getLogger().setLevel(getattr(logging, level.upper(), logging.INFO))

    @classmethod
    def is_initialized(cls) -> bool:
        """Verifica se o logging foi inicializado."""
        return cls._initialized

    @classmethod
    def get_logger(cls, name: str) -> logging.Logger:
        """
        Obtém um logger configurado.

        Args:
            name: Nome do logger (tipicamente __name__)

        Returns:
            Logger configurado
        """
        # Se não foi inicializado, inicializa automaticamente
        if not cls._initialized:
            cls.initialize()

        return logging.getLogger(name)

        return logging.getLogger(name)
