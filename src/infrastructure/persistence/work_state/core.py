"""
SQLAlchemy-based work state persistence.

Orchestrates database access through specialized mixins for each table,
providing better code organization and maintainability.
"""

import logging
import os
import threading
from contextlib import contextmanager
from typing import Optional

from sqlalchemy import create_engine, event
from sqlalchemy.orm import Session, scoped_session, sessionmaker
from sqlalchemy.pool import StaticPool

from .mixins import (
    CombinationCRUDMixin,
    DatasetCRUDMixin,
    DatasetSequenceCRUDMixin,
    ErrorSummary,
    EventCRUDMixin,
    ExecutionCRUDMixin,
    ExecutionDetail,
    ExecutionProgressCRUDMixin,
    ProgressSummary,
    QueriesMixin,
    WorkCRUDMixin,
)
from .models import Base


class WorkPersistence(
    WorkCRUDMixin,
    DatasetCRUDMixin,
    DatasetSequenceCRUDMixin,
    CombinationCRUDMixin,
    ExecutionCRUDMixin,
    ExecutionProgressCRUDMixin,
    EventCRUDMixin,
    QueriesMixin,
):
    """
    Work state persistence using SQLAlchemy ORM with thread-safe operations.

    Orchestrates database access through specialized mixins for each table,
    providing better code organization and maintainability.

    Implements best practices for SQLAlchemy concurrency:
    - Scoped sessions for thread safety
    - Connection pooling for performance
    - Proper transaction management
    - Automatic session cleanup

    Attributes:
        _logger: Logger instance for this class
        _local: Thread-local storage
        _database_url: Database connection URL
        _engine: SQLAlchemy engine instance
        _session_factory: Scoped session factory

    Example:
        >>> persistence = WorkPersistence("sqlite:///work.db")
        >>> with persistence.session_scope() as session:
        ...     work = Work(...)
        ...     session.add(work)
        >>> persistence.close()
    """

    def __init__(self, database_url: Optional[str] = None):
        """
        Initialize SQLAlchemy persistence with database connection.

        Args:
            database_url: Database URL. If None, will try to get from environment
                         or use default SQLite path.

        Raises:
            ValueError: If database URL is invalid
            ConnectionError: If database connection fails
        """
        self._logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._local = threading.local()

        # Get database URL from environment or use default
        if database_url is None:
            database_url = self._get_database_url()

        # Store the database URL for db_path extraction
        self._database_url = database_url

        # Create engine with proper configuration for concurrency
        self._engine = self._create_engine(database_url)

        # Use scoped_session for thread safety
        session_factory = sessionmaker(bind=self._engine)
        self._session_factory = scoped_session(session_factory)

        # Create tables
        self._init_schema()

        self._logger.info(
            f"Initialized SQLAlchemy persistence with URL: {database_url}"
        )

    def _create_engine(self, database_url: str):
        """
        Create SQLAlchemy engine with appropriate configuration.

        Configures engine based on database type (SQLite, PostgreSQL, etc.)
        with optimal settings for concurrency and performance.

        Args:
            database_url: Database connection URL

        Returns:
            Engine: Configured SQLAlchemy engine

        Note:
            SQLite databases get WAL mode and busy timeout configuration.
            Other databases get connection pooling configuration.
        """
        engine_kwargs = {
            "echo": False,  # Set to True for SQL debugging
            "future": True,  # Use SQLAlchemy 2.0+ features
        }

        # Convert URL object to string if needed
        url_str = str(database_url)

        # Configure based on database type
        if url_str.startswith("sqlite:"):
            engine_kwargs.update(
                {
                    "poolclass": StaticPool,
                    "connect_args": {
                        "check_same_thread": False,  # Allow multi-threading
                        "timeout": 30,  # Connection timeout
                    },
                    "pool_pre_ping": True,  # Verify connections before use
                }
            )
        else:
            # For other databases (PostgreSQL, MySQL, etc.)
            engine_kwargs.update(
                {
                    "pool_size": 10,  # Connection pool size
                    "max_overflow": 20,  # Max additional connections
                    "pool_pre_ping": True,  # Verify connections before use
                    "pool_recycle": 3600,  # Recycle connections after 1 hour
                }
            )

        engine = create_engine(database_url, **engine_kwargs)

        # Add SQLite-specific event listeners
        if url_str.startswith("sqlite:"):
            self._configure_sqlite_pragmas(engine)

        return engine

    def _configure_sqlite_pragmas(self, engine) -> None:
        """
        Configure SQLite pragmas for better concurrency and performance.

        Sets up event listeners to configure SQLite-specific settings:
        - WAL mode for better concurrency
        - Normal synchronous mode for performance
        - Foreign key constraints enabled
        - Busy timeout for locked database handling

        Args:
            engine: SQLAlchemy engine for SQLite database
        """

        @event.listens_for(engine, "connect")
        def set_sqlite_pragma(dbapi_connection, connection_record):
            """Configure SQLite for better concurrency and performance."""
            cursor = dbapi_connection.cursor()
            # Enable WAL mode for better concurrency
            cursor.execute("PRAGMA journal_mode=WAL")
            # Set synchronous mode for better performance
            cursor.execute("PRAGMA synchronous=NORMAL")
            # Enable foreign key constraints
            cursor.execute("PRAGMA foreign_keys=ON")
            # Set timeout for busy database
            cursor.execute("PRAGMA busy_timeout=30000")
            cursor.close()

    def _get_database_url(self) -> str:
        """
        Get database URL from environment or use default.

        Checks environment variables in order:
        1. DB_URL - Full database URL
        2. DB_PATH - Path for SQLite database
        3. Default to 'work_state.db' in current directory

        Returns:
            str: Database connection URL
        """
        # Try environment variables first
        db_url = os.getenv("DB_URL")
        if db_url:
            return db_url

        # Default SQLite database
        db_path = os.getenv("DB_PATH", "work_state.db")
        return f"sqlite:///{db_path}"

    def _init_schema(self) -> None:
        """
        Initialize database schema.

        Creates all tables defined in the Base metadata if they don't exist.
        This is safe to call multiple times.
        """
        Base.metadata.create_all(self._engine)

    def get_session(self) -> Session:
        """
        Get a new database session.

        Returns a scoped session that is thread-safe and properly managed.
        The session is tied to the current thread and will be reused for
        subsequent calls within the same thread.

        Returns:
            Session: Thread-scoped SQLAlchemy session

        Note:
            Use session_scope() context manager for automatic transaction management.
        """
        return self._session_factory()

    def get_new_session(self) -> Session:
        """
        Get a new session instance (for testing purposes).

        This creates a new unscoped session instance that is not tied to
        the thread-local scoped session factory. Useful for testing scenarios
        where you need fresh sessions.

        Returns:
            Session: New SQLAlchemy session instance

        Warning:
            Remember to close this session manually when done.
        """
        return sessionmaker(bind=self._engine)()

    @contextmanager
    def session_scope(self):
        """
        Provide a transactional scope around a series of operations.

        This is the recommended way to handle database operations as it
        ensures proper transaction management and cleanup.

        Yields:
            Session: Database session within transaction scope

        Raises:
            Exception: Any exception raised within the context will cause rollback

        Example:
            >>> with persistence.session_scope() as session:
            ...     work = Work(...)
            ...     session.add(work)
            ...     # Automatic commit on success, rollback on exception
        """
        session = self.get_session()
        try:
            yield session
            session.commit()
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()

    def __enter__(self):
        """
        Context manager entry.

        Returns:
            WorkPersistence: Self instance for context manager protocol
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit.

        Ensures proper cleanup of resources when exiting context.

        Args:
            exc_type: Exception type if exception occurred
            exc_val: Exception value if exception occurred
            exc_tb: Exception traceback if exception occurred
        """
        self.close()

    def close(self):
        """
        Close the database engine and clean up resources.

        Removes scoped sessions and disposes of the engine connection pool.
        Call this when shutting down the application to ensure clean resource cleanup.
        """
        if hasattr(self, "_session_factory"):
            self._session_factory.remove()

        if hasattr(self, "_engine"):
            self._engine.dispose()

    # === COMPATIBILITY METHODS ===
    # Methods for compatibility with wrappers that expect specific methods

    def combination_error(
        self, work_id: str, combination_id: str, error: Exception
    ) -> None:
        """
        Register a combination error as an event.

        Compatibility method for CombinationScopedPersistence wrappers.
        Creates an error event associated with the combination for tracking purposes.

        Args:
            work_id: Work identifier
            combination_id: Combination identifier where error occurred
            error: Exception that occurred during combination processing

        Note:
            This method provides backward compatibility with existing wrapper interfaces.
        """
        self.event_create(
            work_id=work_id,
            event_type="error",
            event_category="combination",
            entity_data={
                "combination_id": combination_id,
                "error_type": error.__class__.__name__,
                "error_message": str(error),
            },
        )

    def clear_all_progress_for_non_finalized_executions(self) -> int:
        """
        Clear progress entries for all non-finalized executions system-wide.

        Removes all progress entries for executions that are not finalized
        (status: 'queued', 'running', 'paused') across all works in the system.

        Returns:
            int: Total number of progress entries removed

        Note:
            This is a system-wide cleanup operation that affects all works.
            Use with caution in multi-work environments.
        """
        return self.execution_progress_clear_for_non_finalized()


# Export for external imports
__all__ = [
    "WorkPersistence",
    "ExecutionDetail",
    "ProgressSummary",
    "ErrorSummary",
]
