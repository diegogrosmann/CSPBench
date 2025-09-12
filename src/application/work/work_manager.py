from __future__ import annotations

import logging
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Dict, Optional

from src.domain.config import CSPBenchConfig
from src.domain.status import BaseStatus
from src.domain.work import WorkItem
from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
from src.infrastructure.persistence.work_state import WorkScopedPersistence
from src.infrastructure.persistence.work_state.core import WorkPersistence
from src.infrastructure.utils.path_utils import get_output_base_directory

logger = logging.getLogger(__name__)


class WorkManager:
    """
    WorkManager - Unified Work Management and Execution Service.

    Provides centralized management of work items with persistence support
    and pipeline execution orchestration. Uses isolated database worker pool
    for enhanced thread-safety and process isolation.

    Features:
    - Process-isolated database operations via worker pool
    - Thread-safe operations with RLock
    - Automatic work directory creation
    - Status transition management
    - Persistence through repository pattern
    - Blocking wait functionality for terminal states
    - Pipeline execution orchestration
    - Configuration validation before execution
    - Work restart and resume capabilities
    - Background execution with proper error handling
    - Health monitoring and failure recovery

    Attributes:
        _repo: Repository for work item persistence
        _lock: Thread lock for concurrent access protection

    Example:
        >>> with WorkManager() as manager:
        ...     work_id = manager.submit(config=my_config)
        ...     final_status = manager.wait_until_terminal(work_id)
        ...     work = manager.get(work_id)
    """

    # ==============================================
    # INITIALIZATION AND LIFECYCLE MANAGEMENT
    # ==============================================

    def __init__(self, repository: Optional[WorkPersistence] = None):
        """
        Initialize WorkManager with persistence repository.

        Args:
            repository: WorkPersistence instance for data persistence (optional)
                       If None, a new instance will be created.
                       
        Note:
            The repository handles all database operations and should be
            thread-safe for concurrent access.
        """
        if repository is None:
            repository = WorkPersistence()
        self._repo = repository
        self._lock = threading.RLock()
        
        logger.info("WorkManager initialized with SQLAlchemy persistence")

    def initialize(self) -> Dict[str, Any]:
        """
        Initialize WorkManager and perform minimal cleanup operations.
        
        This method should be called once during application startup to ensure
        consistent state after potential abrupt shutdowns.
        
        Operations performed:
        - Pause works that were in RUNNING status (likely orphaned from abrupt shutdown)
        
        Individual restart operations will handle detailed cleanup when needed.
        
        Returns:
            Dict[str, Any]: Summary of initialization operations performed
                           containing 'works_paused' count and any 'errors'
        """
        with self._lock:
            logger.info("Starting WorkManager initialization...")
            
            summary = {
                "works_paused": 0,
                "errors": []
            }
            
            try:
                # Only pause works that were in RUNNING status (likely orphaned)
                summary["works_paused"] = self._pause_running_works()
                
                logger.info(f"WorkManager initialization completed successfully: {summary}")
                
            except Exception as e:
                error_msg = f"Error during WorkManager initialization: {e}"
                logger.error(error_msg, exc_info=True)
                summary["errors"].append(error_msg)
            
            return summary

    def __enter__(self):
        """
        Context manager entry.
        
        Returns:
            WorkManager: Self instance for context manager protocol
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit - cleanup resources.
        
        Args:
            exc_type: Exception type if exception occurred
            exc_val: Exception value if exception occurred
            exc_tb: Exception traceback if exception occurred
        """
        self.shutdown()

    def shutdown(self, wait: bool = True, timeout: Optional[float] = 30.0) -> None:
        """
        Shutdown WorkManager and cleanup resources.
        
        Gracefully shuts down the WorkManager, optionally waiting for
        pending operations to complete before closing resources.
        
        Args:
            wait: Whether to wait for pending operations
            timeout: Maximum time to wait for shutdown
        """
        with self._lock:            
            if hasattr(self._repo, 'close'):
                self._repo.close()
            
            logger.info("WorkManager shutdown completed")

    # ==============================================
    # CORE WORK ITEM CRUD OPERATIONS
    # ==============================================

    def submit(
        self, *, config: CSPBenchConfig, extra: dict[str, Any] | None = None
    ) -> str:
        """
        Submit new work item for execution.

        Creates a new work item with QUEUED status and persists it.
        Automatically generates work directory and unique ID.

        Args:
            config: CSPBench configuration for the work
            extra: Optional additional metadata to store with work item

        Returns:
            str: Unique work ID for tracking

        Thread Safety:
            This method is thread-safe.
            
        Raises:
            ValueError: If configuration is invalid
            RuntimeError: If work item creation fails
        """
        with self._lock:
            wid = self._new_id()
            work_dir = get_output_base_directory() / wid
            item = WorkItem(
                id=wid,
                config=config,
                status=BaseStatus.QUEUED,
                created_at=self._now(),
                updated_at=self._now(),
                output_path=str(work_dir),
                extra=extra or {},
            )
            self._repo.work_create(
                id=item.id,
                config=item.config.to_dict(),
                status=item.status,
                output_path=item.output_path,
                error=item.error,
                extra=item.extra,
                created_at=item.created_at,
                updated_at=item.updated_at,
            )

            return wid

    def get(self, work_id: str) -> WorkItem | None:
        """
        Retrieve work item by ID.

        Args:
            work_id: Unique work identifier

        Returns:
            WorkItem: Work item instance, None if not found

        Thread Safety:
            This method is thread-safe.
        """
        with self._lock:
            work_data = self._repo.work_get(work_id)
            if work_data is None:
                return None
            return WorkItem.from_dict(work_data)

    def list(
        self, 
        *, 
        filters: dict[str, Any] | None = None,
        order_by: str | None = None,
        order_desc: bool = False
    ) -> list[WorkItem]:
        """
        List all work items with optional filtering.

        This method handles pagination internally to return all work items
        that match the specified filters, not just a single page.

        Args:
            filters: Optional dictionary of field:value pairs for filtering
            order_by: Optional field to order by
            order_desc: Whether to order in descending order

        Returns:
            list[WorkItem]: List of WorkItem instances

        Thread Safety:
            This method is thread-safe.
            
        Note:
            Processes results in batches of 1000 to avoid memory issues
            with large result sets.
        """
        with self._lock:
            all_work_items = []
            offset = 0
            limit = 1000  # Process in batches to avoid memory issues
            
            while True:
                # Get a batch of work items
                work_data_list, total_count = self._repo.work_list(
                    filters=filters,
                    offset=offset,
                    limit=limit,
                    order_by=order_by,
                    order_desc=order_desc
                )
                
                # Convert dictionaries to WorkItem instances
                work_items = [WorkItem.from_dict(work_data) for work_data in work_data_list]
                all_work_items.extend(work_items)
                
                # If we got fewer items than the limit, we've reached the end
                if len(work_data_list) < limit:
                    break
                    
                offset += limit
                
                # Safety check to avoid infinite loops
                if offset > total_count:
                    break
            
            return all_work_items

    # ==============================================
    # STATUS MANAGEMENT OPERATIONS
    # ==============================================

    def get_status(self, work_id: str) -> Optional[str]:
        """
        Get current status of work item.

        Args:
            work_id: Unique work identifier

        Returns:
            str: Current status value, None if work not found

        Thread Safety:
            This method is thread-safe.
        """
        with self._lock:
            item_data = self._repo.work_get(work_id)
            return item_data.get('status') if item_data else None

    def mark_running(self, work_id: str) -> bool:
        """
        Mark work item as running.

        Transitions work item from QUEUED or PAUSED status to RUNNING.
        Updates the timestamp to track when execution started.

        Args:
            work_id: Unique work identifier

        Returns:
            bool: True if status change successful, False otherwise

        Thread Safety:
            This method is thread-safe.
        """
        with self._lock:
            item_data = self._repo.work_get(work_id)
            if not item_data:
                return False
            
            # Update status directly using repository
            try:
                return self._repo.work_update(
                    work_id,
                    status=BaseStatus.RUNNING.value,
                    updated_at=self._now()
                )
            except Exception as e:
                logger.error(f"Failed to update work {work_id} status: {e}")
                return False

    def mark_finished(self, work_id: str) -> bool:
        """
        Mark work item as finished/completed.

        Transitions work item to COMPLETED status indicating successful
        completion of all processing.

        Args:
            work_id: Unique work identifier

        Returns:
            bool: True if status change successful, False otherwise

        Thread Safety:
            This method is thread-safe.
        """
        with self._lock:
            item_data = self._repo.work_get(work_id)
            if not item_data:
                return False
            
            # Update status directly using repository
            try:
                return self._repo.work_update(
                    work_id,
                    status=BaseStatus.COMPLETED.value,
                    updated_at=self._now()
                )
            except Exception as e:
                logger.error(f"Failed to update work {work_id} status: {e}")
                return False

    def mark_completed(self, work_id: str) -> bool:
        """
        Mark work item as completed (alias for mark_finished).

        Args:
            work_id: Unique work identifier

        Returns:
            bool: True if status change successful, False otherwise

        Thread Safety:
            This method is thread-safe.
        """
        return self.mark_finished(work_id)

    def mark_error(self, work_id: str, error: str) -> bool:
        """
        Mark work item as failed with error message.

        Transitions work item to FAILED status and stores the error
        message for debugging and monitoring purposes.

        Args:
            work_id: Unique work identifier
            error: Error message to store

        Returns:
            bool: True if status change successful, False otherwise

        Thread Safety:
            This method is thread-safe.
        """
        with self._lock:
            item_data = self._repo.work_get(work_id)
            if not item_data:
                return False
            
            # Update status directly using repository
            try:
                return self._repo.work_update(
                    work_id,
                    status=BaseStatus.FAILED.value,
                    updated_at=self._now(),
                    error=error
                )
            except Exception as e:
                logger.error(f"Failed to update work {work_id} status: {e}")
                return False

    # ==============================================
    # EXECUTION CONTROL OPERATIONS
    # ==============================================

    def execute(
        self, config: CSPBenchConfig, extra: Optional[Dict[str, Any]] = None
    ) -> str:
        """
        Execute pipeline configuration with unified workflow.

        Validates configuration, submits work item, and starts execution
        in a separate thread for non-blocking operation.

        Args:
            config: Pipeline configuration to execute
            extra: Additional metadata to store with work item

        Returns:
            str: Work ID for tracking execution progress
            
        Raises:
            ValueError: If configuration validation fails
            RuntimeError: If work submission fails
        """
        if extra is None:
            extra = {}

        try:
            # Validate configuration before execution
            self._validate_config(config)

            # Submit work to persistent storage
            work_id = self.submit(config=config, extra=extra)

            logger.info(f"Work submitted: {work_id}")
            logger.info(f"Starting execution thread for work: {work_id}")

            thread = threading.Thread(
                target=self._execute_work, args=(work_id, config), daemon=True
            )
            thread.start()
            
            logger.info(f"Thread started for work: {work_id}")

            return work_id

        except Exception as e:
            logger.error(f"Failed to execute work: {e}")
            raise

    def pause(self, work_id: str) -> bool:
        """
        Pause work item execution.

        Transitions work item from RUNNING or QUEUED status to PAUSED.
        The execution can be resumed later.

        Args:
            work_id: Unique work identifier

        Returns:
            bool: True if status change successful, False otherwise

        Thread Safety:
            This method is thread-safe.
            
        Note:
            Only works in RUNNING or QUEUED status can be paused.
        """
        with self._lock:
            item_data = self._repo.work_get(work_id)
            if not item_data:
                return False
            
            # Check if current status allows pausing
            current_status = item_data.get('status')
            if current_status not in [BaseStatus.RUNNING.value, BaseStatus.QUEUED.value]:
                return False
            
            try:
                return self._repo.work_update(
                    work_id,
                    status=BaseStatus.PAUSED.value,
                    updated_at=self._now()
                )
            except Exception as e:
                logger.error(f"Failed to pause work {work_id}: {e}")
                return False

    def cancel(self, work_id: str) -> bool:
        """
        Cancel work item execution.

        Transitions work item to CANCELED status. This is typically
        irreversible and stops all processing.

        Args:
            work_id: Unique work identifier

        Returns:
            bool: True if status change successful, False otherwise

        Thread Safety:
            This method is thread-safe.
            
        Note:
            Already completed, failed, or canceled works cannot be canceled.
        """
        with self._lock:
            item_data = self._repo.work_get(work_id)
            if not item_data:
                return False
                
            # Check if current status allows canceling
            current_status = item_data.get('status')
            if current_status in [BaseStatus.COMPLETED.value, BaseStatus.FAILED.value, BaseStatus.CANCELED.value]:
                return False
            
            try:
                return self._repo.work_update(
                    work_id,
                    status=BaseStatus.CANCELED.value,
                    updated_at=self._now()
                )
            except Exception as e:
                logger.error(f"Failed to cancel work {work_id}: {e}")
                return False

    def restart(self, work_id: str) -> bool:
        """
        Restart an existing work item execution.

        Retrieves the work item configuration and starts a new execution
        thread after validating the configuration and resetting the status.
        Performs cleanup of non-finalized combinations and executions.

        Args:
            work_id: ID of the work item to restart

        Returns:
            bool: True if restart was successful, False otherwise
            
        Note:
            This operation resets all non-finalized combinations and executions
            to QUEUED status and clears associated progress data.
        """
        with self._lock:
            item_data = self._repo.work_get(work_id)
            if not item_data:
                logger.error(f"Work item {work_id} not found")
                return False
            
            work = WorkItem.from_dict(item_data)
            if not work:
                logger.error(f"Failed to parse work item {work_id}")
                return False
                
            # Check if current status allows restarting
            current_status = item_data.get('status')
            if current_status in [BaseStatus.COMPLETED.value, BaseStatus.FAILED.value]:
                logger.error(f"Cannot restart work {work_id} - status is {current_status}")
                return False
            
            try:
                logger.info(f"Starting restart cleanup for work {work_id}")
                
                # 1. Reset all non-finalized combinations for this work
                combinations_reset = self._reset_work_combinations(work_id)
                if combinations_reset > 0:
                    logger.info(f"Reset {combinations_reset} combinations for work {work_id}")
                
                # 2. Reset all non-finalized executions for this work
                executions_reset = self._reset_work_executions(work_id)
                if executions_reset > 0:
                    logger.info(f"Reset {executions_reset} executions for work {work_id}")
                
                # 3. Clear progress entries for non-finalized executions
                progress_cleared = self._repo.execution_progress_clear_for_non_finalized(work_id=work_id)
                if progress_cleared > 0:
                    logger.info(f"Cleared {progress_cleared} progress entries for work {work_id}")

                # Validate configuration before restart
                try:
                    self._validate_config(work.config)
                except ValueError as e:
                    logger.error(f"Configuration validation failed for work {work_id}: {e}")
                    return False

                logger.info(f"Work restarted: {work_id}")

                # Start execution thread
                thread = threading.Thread(
                    target=self._execute_work, args=(work_id, work.config), daemon=True
                )
                thread.start()

                return True

            except Exception as e:
                logger.error(f"Failed to restart work {work_id}: {e}")
                return False

    def wait_until_terminal(
        self, work_id: str, timeout: float | None = None, poll_interval: float = 0.5
    ) -> Optional[str]:
        """
        Block until work reaches terminal state or timeout.

        Polls work status until it reaches a terminal state (COMPLETED, FAILED, CANCELED)
        or until the specified timeout is reached.

        Args:
            work_id: Unique work identifier
            timeout: Maximum time to wait in seconds, None for no timeout
            poll_interval: Time between status checks in seconds

        Returns:
            str: Final status if terminal state reached
            None: If work item not found
            str: Last observed status if timeout reached (non-terminal)

        Thread Safety:
            This method is thread-safe.
            
        Example:
            >>> final_status = manager.wait_until_terminal("work_123", timeout=3600)
            >>> if final_status == "completed":
            ...     print("Work completed successfully")
        """
        start = self._now()
        last_status: Optional[str] = None
        terminal_statuses = [
            BaseStatus.COMPLETED.value,
            BaseStatus.FAILED.value,
            BaseStatus.CANCELED.value,
        ]

        while True:
            status = self.get_status(work_id)
            if status is None:
                return None
            last_status = status
            if status in terminal_statuses:
                return status

            if timeout is not None and (self._now() - start) >= timeout:
                return last_status

            time.sleep(poll_interval)

    # ==============================================
    # CONFIGURATION AND VALIDATION UTILITIES
    # ==============================================

    def _validate_config(self, config: CSPBenchConfig) -> None:
        """
        Validate configuration before execution.

        Ensures the configuration has all required sections and parameters
        for successful pipeline execution.

        Args:
            config: Configuration to validate

        Raises:
            ValueError: If configuration is invalid or missing required sections
        """
        if not config:
            raise ValueError("Configuration is required")

        # Validate required configuration sections
        if not hasattr(config, 'algorithms') or not config.algorithms:
            raise ValueError("At least one algorithm must be configured")

        if not hasattr(config, 'datasets') or not config.datasets:
            raise ValueError("At least one dataset must be configured")

        # Additional validation can be added here
        logger.debug("Configuration validation passed")

    def get_config_summary(self, config: CSPBenchConfig) -> Dict[str, Any]:
        """
        Generate a summary of the configuration for monitoring and logging.

        Creates a compact summary of the configuration showing algorithms,
        datasets, and total combination count for monitoring purposes.

        Args:
            config: Configuration to summarize

        Returns:
            dict: Configuration summary containing algorithms, datasets,
                 and total_combinations count
        """
        try:
            summary = {
                "algorithms": [],
                "datasets": [],
                "total_combinations": 0
            }

            if hasattr(config, 'algorithms') and config.algorithms:
                summary["algorithms"] = [
                    {
                        "name": algo.name if hasattr(algo, 'name') else str(algo),
                        "type": type(algo).__name__
                    }
                    for algo in config.algorithms
                ]

            if hasattr(config, 'datasets') and config.datasets:
                summary["datasets"] = [
                    {
                        "name": dataset.name if hasattr(dataset, 'name') else str(dataset),
                        "path": getattr(dataset, 'path', 'unknown')
                    }
                    for dataset in config.datasets
                ]

            # Calculate total combinations
            num_algorithms = len(summary["algorithms"])
            num_datasets = len(summary["datasets"])
            summary["total_combinations"] = num_algorithms * num_datasets

            return summary

        except Exception as e:
            logger.warning(f"Failed to generate config summary: {e}")
            return {"error": f"Failed to generate summary: {e}"}

    # ==============================================
    # INTERNAL EXECUTION ENGINE
    # ==============================================

    def _execute_work(self, work_id: str, config: CSPBenchConfig) -> Dict[str, Any]:
        """
        Execute a work item with proper error handling and status management.

        Internal method that runs in a separate thread to execute pipeline
        configuration. Handles status transitions and error reporting.

        Args:
            work_id: Work item identifier
            config: Pipeline configuration

        Returns:
            dict: Execution result summary containing work_id and status
        """
        try:
            logger.info(f"=== Starting work execution thread: {work_id} ===")
            logger.info(f"Thread ID: {threading.current_thread().ident}")
            
            # Get work item
            work = self.get(work_id)
            if not work:
                raise RuntimeError(f"Work item {work_id} not found")

            # Mark as running
            logger.info(f"Marking work {work_id} as running")
            if not self.mark_running(work_id):
                raise RuntimeError(f"Failed to mark work {work_id} as running")

            # Execute the pipeline
            logger.info(f"Starting pipeline execution for work {work_id}")
            self._run_pipeline(work, config)

            return {"work_id": work_id, "status": "completed"}

        except Exception as e:
            error_msg = f"Work execution failed: {e}"
            logger.error(f"Work {work_id} failed: {e}")
            logger.error(f"Full traceback:", exc_info=True)
            
            # Mark as failed
            if not self.mark_error(work_id, error_msg):
                logger.error(f"Failed to mark work {work_id} as failed")

            return {
                "work_id": work_id,
                "status": "failed",
                "error": error_msg
            }

    def _run_pipeline(self, work: WorkItem, config: CSPBenchConfig) -> None:
        """
        Run the pipeline for a work item.

        Creates work-scoped persistence and pipeline runner to execute
        the configured pipeline with proper resource management.

        Args:
            work: Work item to execute
            config: Pipeline configuration
            
        Raises:
            Exception: Any pipeline execution errors are re-raised
        """
        try:
            logger.info(f"Running pipeline for work {work.id}")
            
            # Create work-scoped persistence for the pipeline
            scoped_persistence = WorkScopedPersistence(
                work_id=work.id,
                store=self._repo
            )
            
            # Create pipeline runner and execute
            runner = PipelineRunner(work_store=scoped_persistence)
            
            runner.run(config)
            logger.info(f"Pipeline completed for work {work.id}")
            
        except Exception as e:
            logger.error(f"Pipeline failed for work {work.id}: {e}")
            raise

    # ==============================================
    # INTERNAL CLEANUP AND RESET OPERATIONS
    # ==============================================

    def _pause_running_works(self) -> int:
        """
        Pause all work items that are not in finalized status.
        
        This is used during initialization to handle orphaned works from
        abrupt application shutdowns. Pauses all works that are not in
        finalized states (COMPLETED, FAILED, ERROR).
        
        Returns:
            int: Number of work items paused
        """
        try:
            from src.domain.status import BaseStatus, INCOMPLETE_STATUSES
            
            # Get all non-finalized works
            non_finalized_status_values = [status.value for status in INCOMPLETE_STATUSES]
            works_data, _ = self._repo.work_list(
                filters={'status': non_finalized_status_values}
            )
            
            paused_count = 0
            
            # Update each work using the proper mixin method
            for work_data in works_data:
                work_id = work_data['id']
                current_status = work_data['status']
                
                # Only update if not already paused
                if current_status != BaseStatus.PAUSED.value:
                    success = self._repo.work_update(
                        work_id,
                        status=BaseStatus.PAUSED.value,
                        updated_at=self._now()
                    )
                    
                    if success:
                        paused_count += 1
                        logger.debug(f"Paused work {work_id} (was {current_status})")
                    else:
                        logger.warning(f"Failed to pause work {work_id}")
            
            if paused_count > 0:
                logger.info(f"Paused {paused_count} orphaned work items")
            else:
                logger.debug("No orphaned work items found to pause")
            
            return paused_count
                
        except Exception as e:
            logger.error(f"Error pausing non-finalized works: {e}")
            return 0

    def _clear_non_finalized_progress(self) -> int:
        """
        Clear all progress entries for non-finalized executions.
        
        Removes progress data for executions that are not in final states,
        useful for cleanup during restart operations.
        
        Returns:
            int: Number of progress entries cleared
        """
        try:
            # Use the existing repository method for clearing progress
            cleared_count = self._repo.execution_progress_clear_for_non_finalized()
            logger.info(f"Cleared {cleared_count} progress entries for non-finalized executions")
            return cleared_count
            
        except Exception as e:
            logger.error(f"Error clearing non-finalized progress: {e}")
            return 0

    def _reset_work_combinations(self, work_id: str) -> int:
        """
        Reset all non-finalized combinations for a specific work to queued status.
        
        Resets combinations in RUNNING or PAUSED status back to QUEUED
        and clears timing information for clean restart.
        
        Args:
            work_id: Work identifier to reset combinations for
            
        Returns:
            int: Number of combinations reset
        """
        try:
            from src.domain.status import BaseStatus
            
            non_finalized_statuses = [
                BaseStatus.RUNNING.value,
                BaseStatus.PAUSED.value,
            ]
            
            # Use bulk update for better performance
            reset_count = self._repo.combination_bulk_update(
                filters={
                    'work_id': work_id, 
                    'status': non_finalized_statuses
                },
                update_fields={
                    'status': BaseStatus.QUEUED.value,
                    'started_at': None,
                    'finished_at': None
                }
            )
            
            if reset_count > 0:
                logger.info(f"Reset {reset_count} combinations to queued status for work {work_id}")
            else:
                logger.debug(f"No combinations to reset for work {work_id}")
            
            return reset_count
                
        except Exception as e:
            logger.error(f"Error resetting combinations for work {work_id}: {e}")
            return 0

    def _reset_work_executions(self, work_id: str) -> int:
        """
        Reset all non-finalized executions for a specific work to queued status 
        and clear runtime data.
        
        Resets executions in RUNNING or PAUSED status back to QUEUED,
        clears timing information and results for clean restart.
        
        Args:
            work_id: Work identifier to reset executions for
            
        Returns:
            int: Number of executions reset
        """
        try:
            from src.domain.status import BaseStatus
            
            non_finalized_statuses = [
                BaseStatus.RUNNING.value,
                BaseStatus.PAUSED.value,
            ]
            
            # First, get combination IDs for this work using the proper mixin method
            combinations_data, _ = self._repo.combination_list(
                filters={'work_id': work_id}
            )
            
            if not combinations_data:
                return 0
                
            combination_ids = [combo['id'] for combo in combinations_data]
            
            # Use bulk update for better performance
            reset_count = self._repo.execution_bulk_update(
                filters={
                    'combination_id': combination_ids,
                    'status': non_finalized_statuses
                },
                update_fields={
                    'status': BaseStatus.QUEUED.value,
                    'started_at': None,
                    'finished_at': None,
                    'result': {},  # This will be handled as result_json in the mixin
                    'objective': None
                }
            )
            
            if reset_count > 0:
                logger.info(f"Reset {reset_count} executions to queued status for work {work_id}")
            else:
                logger.debug(f"No executions to reset for work {work_id}")
            
            return reset_count
                
        except Exception as e:
            logger.error(f"Error resetting executions for work {work_id}: {e}")
            return 0

    # ==============================================
    # UTILITY HELPER METHODS
    # ==============================================

    def _now(self) -> float:
        """
        Get current timestamp.
        
        Returns:
            float: Current time as Unix timestamp
        """
        return time.time()

    def _new_id(self) -> str:
        """
        Generate new unique work ID.
        
        Creates a random 12-character hexadecimal identifier
        suitable for use as a work ID.
        
        Returns:
            str: Unique work identifier
        """
        return uuid.uuid4().hex[:12]
