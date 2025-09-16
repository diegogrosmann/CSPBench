"""
Execution Controller

Manages resources, timeouts and status control for task execution.
Centralizes all execution control logic in a single point.
"""

import signal
import threading
import time
from contextlib import contextmanager
from typing import Any, Dict, List, Optional

import psutil

from src.domain.config import ResourcesConfig
from src.domain.status import BaseStatus
from src.infrastructure.logging_config import get_logger
from src.infrastructure.persistence.work_state.core import WorkPersistence


class ExecutionLimitError(Exception):
    """
    Raised when execution limits are exceeded.

    This exception is raised when resource limits (memory, CPU, time)
    are exceeded during execution, allowing for graceful handling
    of resource constraint violations.
    """

    pass


class ExecutionController:
    """
    Controls resources, timeouts and status for task execution.

    Responsibilities:
    - Apply CPU and memory limits
    - Manage timeouts (per-item and total-batch)
    - Control execution status (pause/cancel)
    - Monitor resources during execution

    Attributes:
        _resources: Resource configuration
        _work_id: Work identifier for status checking
        _logger: Logger instance
        _active_processes: List of monitored processes
        _batch_start_time: Batch execution start timestamp
        work_service: Work service for status management
        _max_workers: Maximum number of worker processes
        _exclusive_cores: Whether to use exclusive CPU cores
        _internal_jobs: Number of internal algorithm jobs
        _max_memory_gb: Maximum memory limit in GB
        _timeout_per_item: Timeout per item in seconds
        _timeout_total_batch: Total batch timeout in seconds
        _monitor_thread: Resource monitoring thread
        _stop_monitoring: Event to stop monitoring
        _current_workers: Current number of active workers

    Example:
        >>> controller = ExecutionController("work_123", resources_config)
        >>> with controller.item_timeout(300):
        ...     # Execute item with 300s timeout
        ...     process_item()
        >>> controller.cleanup()
    """

    def __init__(self, work_id: str, resources: Optional[ResourcesConfig] = None):
        """
        Initialize execution controller.

        Args:
            work_id: Work ID to check status from WorkService (required)
            resources: Resource configuration (can be None for defaults)

        Raises:
            ValueError: If work_id is not found
            ImportError: If required system modules are not available
        """
        self._resources = resources
        self._work_id = work_id
        self._logger = get_logger(__name__)
        self._active_processes: List[psutil.Process] = []
        self._batch_start_time = time.time()

        from src.application.services.work_service import get_work_service

        self.work_service = get_work_service()

        if resources is None:
            # Get work details to access output_path
            work_details = self.work_service.get(work_id)
            if work_details is None:
                raise ValueError(f"Work {work_id} not found")

            # Try to get resources from state.db if it exists
            try:
                base_store = WorkPersistence()
                work_data = base_store.work_get(work_id)
                if work_data and work_data.get("config_json"):
                    from src.domain.config import CSPBenchConfig

                    config = CSPBenchConfig.from_dict(work_data["config_json"])
                    if config and config.resources:
                        self._resources = config.resources
                    else:
                        # Use defaults if no config found
                        self._resources = None
                else:
                    # Use defaults if no config found
                    self._resources = None
            except Exception:
                # If state.db doesn't exist or has issues, use defaults
                self._resources = None

        # Extract configurations with defaults
        self._max_workers = None
        self._exclusive_cores = False
        self._internal_jobs = 1
        self._max_memory_gb = None
        self._timeout_per_item = 3600  # 1 hour default
        self._timeout_total_batch = 86400  # 24 hours default

        if resources:
            if resources.cpu:
                self._max_workers = resources.cpu.max_workers
                self._exclusive_cores = resources.cpu.exclusive_cores
                self._internal_jobs = resources.cpu.internal_jobs

            if resources.memory:
                self._max_memory_gb = resources.memory.max_memory_gb

            if resources.timeouts:
                self._timeout_per_item = resources.timeouts.timeout_per_item
                self._timeout_total_batch = resources.timeouts.timeout_total_batch

        # Apply defaults
        if self._max_workers is None:
            total_cores = psutil.cpu_count() or 1
            if self._exclusive_cores and total_cores > 1:
                # Reserve core 0 for main application, use remaining cores for algorithms
                self._max_workers = total_cores - 1
            else:
                self._max_workers = total_cores

        # Apply CPU affinity to main application if exclusive_cores is enabled
        if self._exclusive_cores:
            self._apply_main_process_affinity()

        # Monitor thread
        self._monitor_thread = None
        self._stop_monitoring = threading.Event()
        self._current_workers = (
            self._max_workers
        )  # Track current active workers for memory control

        self._logger.info(
            f"[EXECUTION] Controller initialized - "
            f"max_workers={self._max_workers}, memory={self._max_memory_gb}GB, "
            f"exclusive_cores={self._exclusive_cores}, internal_jobs={self._internal_jobs}, "
            f"timeout_item={self._timeout_per_item}s, timeout_batch={self._timeout_total_batch}s"
        )

    @property
    def max_workers(self) -> int:
        """
        Get max workers for executor parallelization (considering memory limits).

        Returns the current maximum number of workers, which may be less than
        the configured maximum if memory constraints have forced a reduction.

        Returns:
            int: Current maximum number of workers
        """
        return self._current_workers

    @property
    def configured_max_workers(self) -> int:
        """
        Get originally configured max workers (ignoring memory adjustments).

        Returns the originally configured maximum workers before any
        dynamic adjustments due to memory constraints.

        Returns:
            int: Originally configured maximum workers
        """
        return self._max_workers

    @property
    def internal_jobs(self) -> int:
        """
        Get internal jobs for algorithm parallelization.

        Returns internal_jobs limited by current memory constraints.
        When memory forces reduced workers, algorithms should use proportional
        internal parallelism.

        Returns:
            int: Current number of internal jobs (limited by current workers)
        """
        # Limit internal_jobs to current_workers to prevent thread creation beyond capacity
        return min(self._internal_jobs, self._current_workers)

    @classmethod
    def from_config(cls, work_id: str, config: dict[str, Any]) -> "ExecutionController":
        """
        Create ExecutionController from configuration dictionary.

        Factory method to construct ExecutionController from a dictionary
        configuration, reconstructing ResourcesConfig objects from nested dicts.

        Args:
            work_id: Work identifier for status checking
            config: Configuration dictionary with cpu, memory, timeouts sections

        Returns:
            ExecutionController: Configured controller instance

        Example:
            >>> config = {
            ...     "cpu": {"max_workers": 4, "exclusive_cores": True},
            ...     "memory": {"max_memory_gb": 8},
            ...     "timeouts": {"timeout_per_algorithm": 3600}
            ... }
            >>> controller = ExecutionController.from_config("work_123", config)
        """
        from src.domain.config import (
            CPUConfig,
            MemoryConfig,
            ResourcesConfig,
            TimeoutsConfig,
        )

        # Reconstruct ResourcesConfig from dictionary
        cpu_cfg = None
        if "cpu" in config:
            cpu_cfg = CPUConfig(
                max_workers=config["cpu"].get("max_workers"),
                exclusive_cores=config["cpu"].get("exclusive_cores", False),
                internal_jobs=config["cpu"].get("internal_jobs", 1),
            )

        memory_cfg = None
        if "memory" in config:
            memory_cfg = MemoryConfig(
                max_memory_gb=config["memory"].get("max_memory_gb")
            )

        timeouts_cfg = None
        if "timeouts" in config:
            timeouts_cfg = TimeoutsConfig(
                timeout_per_item=config["timeouts"].get("timeout_per_algorithm", 3600),
                timeout_total_batch=config["timeouts"].get(
                    "timeout_total_batch", 86400
                ),
            )

        resources = ResourcesConfig(
            cpu=cpu_cfg,
            memory=memory_cfg,
            timeouts=timeouts_cfg,
        )

        return cls(work_id=work_id, resources=resources)

    @property
    def timeout_per_item(self) -> int:
        """
        Get timeout per item in seconds.

        Returns:
            int: Timeout for individual item processing
        """
        return self._timeout_per_item

    @property
    def timeout_total_batch(self) -> int:
        """
        Get total batch timeout in seconds.

        Returns:
            int: Total timeout for entire batch processing
        """
        return self._timeout_total_batch

    def get_batch_elapsed_time(self) -> float:
        """
        Get elapsed time since batch start.

        Returns:
            float: Elapsed time in seconds since batch execution started
        """
        return time.time() - self._batch_start_time

    def check_status(self) -> BaseStatus:
        """
        Check current execution status.

        Queries the WorkService for the current status of the work item.
        Falls back to RUNNING status if status check fails.

        Returns:
            BaseStatus: Current execution status from WorkService

        Note:
            Returns RUNNING as default if status cannot be determined.
        """
        try:
            status = self.work_service.get_status(self._work_id)
            if status is not None:
                # Convert string status to BaseStatus if needed
                if isinstance(status, str):
                    try:
                        return BaseStatus(status)
                    except ValueError:
                        self._logger.warning(
                            f"Invalid status from WorkService: {status}"
                        )
                        return BaseStatus.RUNNING
                return status
        except Exception as e:
            self._logger.warning(f"Error checking WorkService status: {e}")

        # Default to RUNNING if status check fails
        return BaseStatus.RUNNING

    def wait_for_continue(self) -> BaseStatus:
        """
        Wait while paused, return when status changes.

        Blocks execution while work is in PAUSED status, polling
        at regular intervals until status changes to something else.

        Returns:
            BaseStatus: New execution status (RUNNING or CANCELED)
        """
        while self.check_status() == BaseStatus.PAUSED:
            time.sleep(0.3)

        return self.check_status()

    def check_batch_timeout(self) -> None:
        """
        Check if total batch timeout has been exceeded.

        Raises:
            ExecutionLimitError: If batch timeout exceeded

        Note:
            Also triggers cleanup of all active processes when timeout is reached.
        """
        elapsed = self.get_batch_elapsed_time()
        if elapsed > self._timeout_total_batch:
            # Force cleanup and termination of all active processes
            self._force_cleanup_on_timeout()
            raise ExecutionLimitError(
                f"Batch execution exceeded {self._timeout_total_batch} seconds"
            )

    def _force_cleanup_on_timeout(self) -> None:
        """
        Force cleanup and termination of all processes when batch timeout is reached.

        This method is called when batch timeout is exceeded to ensure
        all active processes are properly terminated and cleaned up.
        """
        self._logger.warning(
            "[EXECUTION] Batch timeout reached - forcing cleanup of all processes"
        )

        # Stop monitoring
        self._stop_monitoring.set()

        # Terminate all registered processes
        terminated_count = 0
        for process in self._active_processes:
            try:
                if process.is_running():
                    process.terminate()
                    terminated_count += 1
                    # Give process 1 second to terminate gracefully
                    try:
                        process.wait(timeout=1)
                    except psutil.TimeoutExpired:
                        # Force kill if not terminated
                        process.kill()
                        self._logger.warning(
                            f"[EXECUTION] Force killed process {process.pid}"
                        )
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass

        if terminated_count > 0:
            self._logger.warning(
                f"[EXECUTION] Terminated {terminated_count} active processes due to batch timeout"
            )

        # Clear the list
        self._active_processes.clear()

    def get_cpu_affinity(self) -> Optional[List[int]]:
        """
        Calculate CPU affinity based on exclusive_cores setting.

        When exclusive_cores is enabled, reserves core 0 for the main application
        and returns the remaining cores for algorithm processes.

        Returns:
            List[int]: List of CPU cores to use for algorithm processes,
                      or None if no CPU restriction should be applied

        Note:
            Returns cores 1 to max_workers when exclusive_cores is True,
            or None when exclusive_cores is False.
        """
        if not self._exclusive_cores:
            return None

        try:
            total_cores = psutil.cpu_count()
            if total_cores <= 1:
                return [0]  # Only one core available

            # Exclude core 0 for main process when exclusive_cores is enabled
            # Algorithm processes use cores 1 to max_workers
            return list(range(1, min(total_cores, self._max_workers + 1)))
        except Exception as e:
            self._logger.warning(f"[EXECUTION] Cannot calculate CPU affinity: {e}")
            return None

    def _apply_main_process_affinity(self) -> None:
        """
        Apply CPU affinity to main application process when exclusive_cores is enabled.

        Forces main application to run only on core 0 when exclusive_cores
        configuration is active, leaving other cores for algorithm execution.
        """
        try:
            total_cores = psutil.cpu_count()
            if total_cores <= 1:
                return  # Nothing to do with single core

            main_process = psutil.Process()
            main_process.cpu_affinity([0])  # Force main application to core 0
            self._logger.info(
                "[EXECUTION] Main application restricted to core 0 (exclusive_cores=True)"
            )
        except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError) as e:
            self._logger.warning(
                f"[EXECUTION] Cannot apply main process CPU affinity: {e}"
            )

    def apply_cpu_limits(self, process: Optional[psutil.Process] = None) -> None:
        """
        Apply CPU limits to current process or specified process.

        Sets CPU affinity and priority based on the controller configuration.
        Applies exclusive core assignments and priority adjustments.

        Args:
            process: Process to apply limits to (current process if None)

        Note:
            Logs warnings if CPU limits cannot be applied due to system restrictions.
        """
        if process is None:
            process = psutil.Process()

        try:
            # Apply CPU affinity if exclusive_cores is enabled
            cpu_affinity = self.get_cpu_affinity()

            if cpu_affinity:
                available_cpus = list(range(psutil.cpu_count()))
                valid_affinity = [cpu for cpu in cpu_affinity if cpu in available_cpus]
                if valid_affinity:
                    process.cpu_affinity(valid_affinity)
                    self._logger.info(
                        f"[EXECUTION] CPU affinity set to: {valid_affinity}"
                    )
                else:
                    self._logger.warning(
                        f"[EXECUTION] Invalid CPU affinity: {cpu_affinity}"
                    )

            # Apply CPU priority adjustment if max_workers is limited
            if self._max_workers < psutil.cpu_count():
                # Lower priority for processes when core limit is set
                process.nice(10)  # Higher nice value = lower priority
                self._logger.info(
                    f"[EXECUTION] CPU priority adjusted for max_workers={self._max_workers}"
                )

        except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError) as e:
            self._logger.warning(f"[EXECUTION] Cannot apply CPU limits: {e}")

    def apply_memory_limits(self, process: Optional[psutil.Process] = None) -> None:
        """
        Apply memory limits using system resource limits and start preventive monitoring.

        Sets system-level memory limits and starts a background monitoring thread
        that dynamically adjusts worker count based on memory usage.

        Args:
            process: Process to monitor (current process if None)

        Note:
            Uses system resource limits where available and falls back to
            preventive monitoring that reduces workers instead of killing processes.
        """
        if not self._max_memory_gb:
            return

        try:
            import resource

            # Convert GB to bytes
            max_memory_bytes = int(self._max_memory_gb * 1024 * 1024 * 1024)

            # Set memory limit using resource module
            resource.setrlimit(resource.RLIMIT_AS, (max_memory_bytes, max_memory_bytes))

            self._logger.info(
                f"[EXECUTION] Memory limit set to {self._max_memory_gb}GB ({max_memory_bytes} bytes)"
            )

        except (ImportError, OSError, ValueError) as e:
            self._logger.warning(f"[EXECUTION] Cannot apply memory limit: {e}")

        # Always start preventive memory monitoring
        self._start_preventive_memory_monitor(process)

    def _start_preventive_memory_monitor(
        self, process: Optional[psutil.Process] = None
    ) -> None:
        """
        Start preventive memory monitoring thread that reduces workers instead of killing process.

        Monitors memory usage and dynamically adjusts the number of workers
        to prevent memory limit violations. Provides graceful degradation
        instead of process termination.

        Args:
            process: Process to monitor (current process if None)

        Note:
            The monitor runs in a background thread and adjusts _current_workers
            based on memory pressure thresholds (85% warning, 95% critical).
        """
        if not self._max_memory_gb:
            return

        if process is None:
            process = psutil.Process()

        def monitor_memory():
            """Monitor memory usage and reduce workers when approaching limit."""
            max_memory_bytes = self._max_memory_gb * 1024 * 1024 * 1024
            warning_threshold = 0.85  # 85% of limit
            critical_threshold = 0.95  # 95% of limit

            while not self._stop_monitoring.is_set():
                try:
                    memory_info = process.memory_info()
                    current_memory = memory_info.rss  # Resident Set Size
                    memory_usage_ratio = current_memory / max_memory_bytes

                    if memory_usage_ratio >= critical_threshold:
                        # Critical: Reduce workers immediately to minimum (1)
                        if self._current_workers > 1:
                            old_workers = self._current_workers
                            self._current_workers = 1
                            self._logger.warning(
                                f"[EXECUTION] CRITICAL memory usage {memory_usage_ratio:.1%} "
                                f"({current_memory / (1024**3):.2f}GB / {self._max_memory_gb}GB). "
                                f"Reducing workers: {old_workers} → {self._current_workers}"
                            )
                    elif memory_usage_ratio >= warning_threshold:
                        # Warning: Reduce workers by half (minimum 1)
                        if self._current_workers > 1:
                            old_workers = self._current_workers
                            self._current_workers = max(1, self._current_workers // 2)
                            self._logger.warning(
                                f"[EXECUTION] High memory usage {memory_usage_ratio:.1%} "
                                f"({current_memory / (1024**3):.2f}GB / {self._max_memory_gb}GB). "
                                f"Reducing workers: {old_workers} → {self._current_workers}"
                            )
                    elif (
                        memory_usage_ratio < 0.5
                        and self._current_workers < self._max_workers
                    ):
                        # Memory usage low: Can safely increase workers back
                        old_workers = self._current_workers
                        self._current_workers = min(
                            self._max_workers, self._current_workers * 2
                        )
                        self._logger.info(
                            f"[EXECUTION] Memory usage normalized {memory_usage_ratio:.1%}. "
                            f"Increasing workers: {old_workers} → {self._current_workers}"
                        )

                    time.sleep(2)  # Check every 2 seconds

                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    break
                except Exception as e:
                    self._logger.warning(f"[EXECUTION] Memory monitor error: {e}")
                    break

        self._monitor_thread = threading.Thread(target=monitor_memory, daemon=True)
        self._monitor_thread.start()
        self._logger.info("[EXECUTION] Preventive memory monitoring started")

    @contextmanager
    def item_timeout(self, timeout_seconds: Optional[int] = None):
        """
        Context manager for individual item timeout enforcement.

        Uses SIGALRM to enforce timeout on individual item processing.
        Automatically restores previous signal handler when done.

        Args:
            timeout_seconds: Timeout in seconds (uses config default if None)

        Yields:
            None: Context for timeout-controlled execution

        Raises:
            TimeoutError: If item execution exceeds timeout

        Example:
            >>> with controller.item_timeout(300):
            ...     process_long_running_item()
        """
        if timeout_seconds is None:
            timeout_seconds = self._timeout_per_item

        def timeout_handler(signum, frame):
            raise TimeoutError(f"Item execution exceeded {timeout_seconds} seconds")

        # Set alarm for timeout
        old_handler = signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout_seconds)

        try:
            self._logger.debug(f"[EXECUTION] Item timeout set: {timeout_seconds}s")
            yield
        finally:
            # Cancel alarm and restore handler
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)

    def create_worker_config(self) -> Dict[str, Any]:
        """
        Create configuration for worker processes.

        Builds a configuration dictionary containing current resource limits
        and settings for use by ProcessPool workers.

        Returns:
            Dict[str, Any]: Configuration dictionary for ProcessPool workers
                           containing cpu, memory, and timeout settings
        """
        return {
            "cpu": {
                "max_workers": self._current_workers,  # Use current workers (considering memory limits)
                "exclusive_cores": self._exclusive_cores,
                "affinity": self.get_cpu_affinity(),
            },
            "memory": {"max_memory_gb": self._max_memory_gb},
            "timeouts": {
                "timeout_per_algorithm": self._timeout_per_item,
                "timeout_total_batch": self._timeout_total_batch,
            },
        }

    def get_worker_config(self) -> Dict[str, Any]:
        """
        Alias for create_worker_config for backward compatibility.

        Returns:
            Dict[str, Any]: Configuration dictionary for ProcessPool workers
        """
        return self.create_worker_config()

    def register_process(self, process: psutil.Process) -> None:
        """
        Register a process for resource monitoring.

        Adds process to the list of monitored processes and applies
        CPU limits to it. Registered processes will be cleaned up
        on timeout or controller shutdown.

        Args:
            process: Process to register and monitor
        """
        self._active_processes.append(process)
        self.apply_cpu_limits(process)

    def cleanup(self) -> None:
        """
        Cleanup execution controller and stop monitoring.

        Stops all monitoring threads and terminates any remaining
        active processes. Should be called when execution is complete
        or when shutting down.
        """
        self._stop_monitoring.set()

        if self._monitor_thread and self._monitor_thread.is_alive():
            self._monitor_thread.join(timeout=2.0)

        # Terminate any remaining processes
        for process in self._active_processes:
            try:
                if process.is_running():
                    process.terminate()
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass

        self._active_processes.clear()
        self._logger.info("[EXECUTION] Execution controller cleanup completed")

    def get_resource_info(self) -> Dict[str, Any]:
        """
        Get current resource usage information.

        Provides detailed information about current resource usage,
        limits, and execution state for monitoring and debugging.

        Returns:
            Dict[str, Any]: Resource information including CPU usage,
                           memory usage, thread count, elapsed time,
                           and all configured limits
        """
        try:
            process = psutil.Process()
            memory_info = process.memory_info()

            return {
                "cpu_percent": process.cpu_percent(),
                "memory_mb": memory_info.rss / (1024 * 1024),
                "memory_percent": process.memory_percent(),
                "num_threads": process.num_threads(),
                "batch_elapsed_time": self.get_batch_elapsed_time(),
                "status": self.check_status().value,
                "limits": {
                    "configured_max_workers": self._max_workers,
                    "current_max_workers": self._current_workers,
                    "exclusive_cores": self._exclusive_cores,
                    "configured_internal_jobs": self._internal_jobs,
                    "current_internal_jobs": self.internal_jobs,  # Show both configured and current values
                    "max_memory_gb": self._max_memory_gb,
                    "timeout_per_item": self._timeout_per_item,
                    "timeout_total_batch": self._timeout_total_batch,
                },
            }
        except Exception as e:
            self._logger.warning(f"[EXECUTION] Cannot get resource info: {e}")
            return {"error": str(e)}
