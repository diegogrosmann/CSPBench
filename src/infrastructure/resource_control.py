"""
Resource Control Manager

Implements actual resource limits enforcement for CSPBench.
This module provides mechanisms to apply CPU, memory, and timeout constraints.
"""

import signal
import threading
import time
from contextlib import contextmanager
from typing import Any, Dict, List, Optional

import psutil

from src.infrastructure.logging_config import get_logger


class ResourceLimitError(Exception):
    """Raised when resource limits are exceeded."""

    pass


class ResourceController:
    """Controls and enforces resource limits for algorithm execution."""

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize resource controller.

        Args:
            config: Resource configuration dictionary containing:
                - cpu: {max_cores, affinity, max_workers, exclusive_cores}
                - memory: {max_memory_gb}
                - parallel: {enabled, max_workers, internal_jobs}
                - timeouts: {timeout_per_algorithm, timeout_total_batch}
        """
        self._config = config
        self._logger = get_logger(__name__)
        self._active_processes: List[psutil.Process] = []
        self._start_time = time.time()

        # Extract configurations
        self._cpu_config = config.get("cpu", {})
        self._memory_config = config.get("memory", {})
        self._parallel_config = config.get("parallel", {})
        self._timeout_config = config.get("timeouts", {})

        # Resource limits with defaults
        self._max_cores = self._cpu_config.get("max_cores") or self._cpu_config.get("max_workers")
        if self._max_cores is None:
            # Default to number of CPU cores when max_workers is None
            self._max_cores = psutil.cpu_count()
        
        self._cpu_affinity = self._cpu_config.get("affinity")
        self._exclusive_cores = self._cpu_config.get("exclusive_cores", False)
        self._max_memory_gb = self._memory_config.get("max_memory_gb")
        self._timeout_per_algorithm = self._timeout_config.get(
            "timeout_per_algorithm", 3600
        )
        self._timeout_total_batch = self._timeout_config.get(
            "timeout_total_batch", 86400
        )

        # Monitor thread
        self._monitor_thread = None
        self._stop_monitoring = threading.Event()

        self._logger.info(
            f"[RESOURCE] Controller initialized with limits: "
            f"cores={self._max_cores}, memory={self._max_memory_gb}GB, "
            f"exclusive_cores={self._exclusive_cores}, "
            f"timeout_algo={self._timeout_per_algorithm}s, timeout_batch={self._timeout_total_batch}s"
        )

    def get_default_cpu_affinity(self) -> Optional[List[int]]:
        """
        Calculate default CPU affinity based on exclusive_cores setting.
        
        Returns:
            List of CPU cores to use, or None if no restriction
        """
        if not self._exclusive_cores:
            return None
            
        try:
            total_cores = psutil.cpu_count()
            if total_cores <= 1:
                return [0]  # Only one core available
            
            # Exclude core 0 for main process when exclusive_cores is enabled
            return list(range(1, min(total_cores, self._max_cores + 1)))
        except Exception as e:
            self._logger.warning(f"[RESOURCE] Cannot calculate CPU affinity: {e}")
            return None

    def apply_cpu_limits(self, process: Optional[psutil.Process] = None) -> None:
        """
        Apply CPU limits to current process or specified process.

        Args:
            process: Process to apply limits to (current process if None)
        """
        if process is None:
            process = psutil.Process()

        try:
            # Apply CPU affinity - either specified or calculated from exclusive_cores
            affinity_to_use = self._cpu_affinity or self.get_default_cpu_affinity()
            
            if affinity_to_use:
                available_cpus = list(range(psutil.cpu_count()))
                valid_affinity = [
                    cpu for cpu in affinity_to_use if cpu in available_cpus
                ]
                if valid_affinity:
                    process.cpu_affinity(valid_affinity)
                    self._logger.info(
                        f"[RESOURCE] CPU affinity set to: {valid_affinity}"
                    )
                else:
                    self._logger.warning(
                        f"[RESOURCE] Invalid CPU affinity: {affinity_to_use}"
                    )

            # Apply CPU core limit (nice value - not perfect but helps)
            if self._max_cores and self._max_cores < psutil.cpu_count():
                # Lower priority for processes when core limit is set
                process.nice(10)  # Higher nice value = lower priority
                self._logger.info(
                    f"[RESOURCE] CPU limit applied: max_cores={self._max_cores}"
                )

        except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError) as e:
            self._logger.warning(f"[RESOURCE] Cannot apply CPU limits: {e}")

    def apply_memory_limits(self, process: Optional[psutil.Process] = None) -> None:
        """
        Apply memory limits using system resource limits.

        Args:
            process: Process to monitor (current process if None)
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
                f"[RESOURCE] Memory limit set to {self._max_memory_gb}GB ({max_memory_bytes} bytes)"
            )

        except (ImportError, OSError, ValueError) as e:
            self._logger.warning(f"[RESOURCE] Cannot apply memory limit: {e}")
            # Fallback to monitoring approach
            self._start_memory_monitor(process)

    def _start_memory_monitor(self, process: Optional[psutil.Process] = None) -> None:
        """Start memory monitoring thread for enforcement."""
        if not self._max_memory_gb:
            return

        if process is None:
            process = psutil.Process()

        def monitor_memory():
            """Monitor memory usage and terminate if exceeded."""
            max_memory_bytes = self._max_memory_gb * 1024 * 1024 * 1024

            while not self._stop_monitoring.is_set():
                try:
                    memory_info = process.memory_info()
                    current_memory = memory_info.rss  # Resident Set Size

                    if current_memory > max_memory_bytes:
                        self._logger.error(
                            f"[RESOURCE] Memory limit exceeded: "
                            f"{current_memory / (1024**3):.2f}GB > {self._max_memory_gb}GB"
                        )

                        # Terminate the process
                        process.terminate()
                        time.sleep(1)
                        if process.is_running():
                            process.kill()

                        raise ResourceLimitError(
                            f"Memory limit exceeded: {self._max_memory_gb}GB"
                        )

                    time.sleep(1)  # Check every second

                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    break
                except Exception as e:
                    self._logger.warning(f"[RESOURCE] Memory monitor error: {e}")
                    break

        self._monitor_thread = threading.Thread(target=monitor_memory, daemon=True)
        self._monitor_thread.start()

    @contextmanager
    def algorithm_timeout(self, timeout_seconds: Optional[int] = None):
        """
        Context manager for algorithm timeout enforcement.

        Args:
            timeout_seconds: Timeout in seconds (uses config default if None)
        """
        if timeout_seconds is None:
            timeout_seconds = self._timeout_per_algorithm

        def timeout_handler(signum, frame):
            raise TimeoutError(
                f"Algorithm execution exceeded {timeout_seconds} seconds"
            )

        # Set alarm for timeout
        old_handler = signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout_seconds)

        try:
            self._logger.debug(f"[RESOURCE] Algorithm timeout set: {timeout_seconds}s")
            yield
        finally:
            # Cancel alarm and restore handler
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)

    def check_batch_timeout(self) -> None:
        """Check if total batch timeout has been exceeded."""
        elapsed = time.time() - self._start_time
        if elapsed > self._timeout_total_batch:
            raise TimeoutError(
                f"Batch execution exceeded {self._timeout_total_batch} seconds"
            )

    def register_process(self, process: psutil.Process) -> None:
        """Register a process for resource monitoring."""
        self._active_processes.append(process)
        self.apply_cpu_limits(process)

    def cleanup(self) -> None:
        """Cleanup resource controller and stop monitoring."""
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
        self._logger.info("[RESOURCE] Resource controller cleanup completed")

    def get_resource_info(self) -> Dict[str, Any]:
        """Get current resource usage information."""
        try:
            process = psutil.Process()
            memory_info = process.memory_info()

            return {
                "cpu_percent": process.cpu_percent(),
                "memory_mb": memory_info.rss / (1024 * 1024),
                "memory_percent": process.memory_percent(),
                "num_threads": process.num_threads(),
                "elapsed_time": time.time() - self._start_time,
                "limits": {
                    "max_cores": self._max_cores,
                    "cpu_affinity": self._cpu_affinity,
                    "max_memory_gb": self._max_memory_gb,
                    "timeout_per_algorithm": self._timeout_per_algorithm,
                    "timeout_total_batch": self._timeout_total_batch,
                },
            }
        except Exception as e:
            self._logger.warning(f"[RESOURCE] Cannot get resource info: {e}")
            return {"error": str(e)}


def create_resource_controller(config: Dict[str, Any]) -> ResourceController:
    """
    Factory function to create resource controller.

    Args:
        config: Resource configuration

    Returns:
        ResourceController instance
    """
    return ResourceController(config)
