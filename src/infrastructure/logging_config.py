"""
Logging Configuration for CSPBench

This module provides centralized logging configuration for the CSPBench framework.
It supports file-based and console logging with configurable levels, rotation,
and formatting. The configuration is environment-variable driven for flexibility
across different deployment scenarios.

Features:
    - Environment-variable driven configuration
    - Rotating file handlers with size limits
    - Configurable log levels and formatting
    - Console and file output options
    - Thread-safe logger initialization
    - Centralized logger factory

Environment Variables:
    LOG_LEVEL: Logging level (DEBUG, INFO, WARNING, ERROR)
    LOG_DIRECTORY: Directory for log files
    LOG_BASE_NAME: Base filename for log files
    LOG_MAX_BYTES: Maximum file size before rotation
    LOG_BACKUP_COUNT: Number of backup files to keep
    LOG_TO_STDOUT: Whether to log to console instead of file

Example:
    Basic setup using environment variables::

        import os
        os.environ['LOG_LEVEL'] = 'INFO'
        os.environ['LOG_DIRECTORY'] = 'logs'

        from src.infrastructure.logging_config import setup_logging_from_env
        setup_logging_from_env()

        logger = get_logger(__name__)
        logger.info("Application started")
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
    log_file_path: Optional[str] = None,
    log_to_stdout: Optional[bool] = None,
) -> None:
    """Configure basic system logging using environment variables as defaults.

    This function sets up the root logger with either file-based or console output,
    configurable formatting, and optional log rotation. All parameters have
    sensible defaults that can be overridden via environment variables.

    Args:
        level: Log level (DEBUG, INFO, WARNING, ERROR) - defaults to LOG_LEVEL env var
        log_dir: Directory for log files - defaults to LOG_DIRECTORY env var
        base_name: Base name for log files - defaults to LOG_BASE_NAME env var
        max_bytes: Maximum file size before rotation - defaults to LOG_MAX_BYTES env var
        backup_count: Number of backups to maintain - defaults to LOG_BACKUP_COUNT env var
        log_file_path: Full file path (overrides log_dir + base_name combination)
        log_to_stdout: Whether to log to console instead of file - defaults to LOG_TO_STDOUT env var

    Note:
        This function configures the root logger and will affect all loggers
        in the application. Call only once during application initialization.
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
    """Configure logging using only environment variables.

    This is a convenience function that sets up logging entirely from .env configuration.
    It's equivalent to calling setup_basic_logging() with no arguments, which will
    use all environment variable defaults.

    Environment Variables Used:
        - LOG_LEVEL: Logging level (default: INFO)
        - LOG_DIRECTORY: Log file directory (default: outputs/logs)
        - LOG_BASE_NAME: Log file base name (default: cspbench)
        - LOG_MAX_BYTES: Max file size for rotation (default: 10485760)
        - LOG_BACKUP_COUNT: Number of backup files (default: 5)
        - LOG_TO_STDOUT: Log to console instead of file (default: false)
    """
    # Configure logging using only environment variables.

    setup_basic_logging()


def get_logger(name: str) -> logging.Logger:
    """Get a configured logger for a specific module.

    Returns a logger instance with the specified name. The logger will inherit
    the configuration set up by setup_basic_logging() or setup_logging_from_env().

    Args:
        name: Module/logger name, typically __name__ from the calling module

    Returns:
        Configured logger instance ready for use

    Example:
        module_logger = get_logger(__name__)
        module_logger.info("Module initialized")
    """
    return logging.getLogger(name)


class LoggerConfig:
    """Centralized logger configuration manager.

    This class provides a singleton-like interface for managing logging
    configuration across the application. It ensures logging is initialized
    only once and provides methods for runtime configuration changes.

    Attributes:
        _initialized: Whether logging has been set up
        _log_level: Current logging level
    """

    _initialized = False
    _log_level = None

    @classmethod
    def initialize(cls, level: Optional[str] = None, **kwargs) -> None:
        """Initialize logging configuration.

        Sets up logging if not already initialized. This method is idempotent
        and safe to call multiple times.

        Args:
            level: Logging level to use (defaults to LOG_LEVEL env var)
            **kwargs: Additional arguments passed to setup_basic_logging()
        """
        if not cls._initialized:
            # Use environment variable if level not provided
            level = level or os.getenv("LOG_LEVEL", "INFO")
            setup_basic_logging(level=level, **kwargs)
            cls._log_level = level
            cls._initialized = True

    @classmethod
    def get_level(cls) -> str:
        """Get the current logging level.

        Returns:
            Current logging level as string (e.g., 'INFO', 'DEBUG')
        """
        return cls._log_level or os.getenv("LOG_LEVEL", "INFO")

    @classmethod
    def set_level(cls, level: str) -> None:
        """Change the logging level at runtime.

        Updates both the class state and the actual logger configuration.

        Args:
            level: New logging level (DEBUG, INFO, WARNING, ERROR)
        """
        cls._log_level = level
        logging.getLogger().setLevel(getattr(logging, level.upper(), logging.INFO))

    @classmethod
    def is_initialized(cls) -> bool:
        """Check if logging has been initialized.

        Returns:
            True if logging configuration has been set up, False otherwise
        """
        return cls._initialized

    @classmethod
    def get_logger(cls, name: str) -> logging.Logger:
        """Get a logger with automatic initialization.

        Returns a configured logger, automatically initializing the logging
        system if it hasn't been set up yet. This is the recommended way
        to get loggers in application code.

        Args:
            name: Logger name, typically __name__ from the calling module

        Returns:
            Configured logger instance

        Example:
            logger = LoggerConfig.get_logger(__name__)
            logger.info("Logger obtained with auto-initialization")
        """
        # Se não foi inicializado, inicializa automaticamente
        if not cls._initialized:
            cls.initialize()

        return logging.getLogger(name)
