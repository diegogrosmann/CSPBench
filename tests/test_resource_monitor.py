"""Tests for ResourceMonitor module."""

from unittest.mock import MagicMock

from src.utils.resource_monitor import (
    ResourceLimits,
    ResourceMonitor,
    check_algorithm_feasibility,
    estimate_algorithm_memory,
    force_garbage_collection,
    get_available_memory_mb,
    get_safe_memory_limit,
)


class TestResourceLimits:
    """Test ResourceLimits functionality."""

    def test_from_config(self):
        """Test creating ResourceLimits from config."""
        limits = ResourceLimits.from_config()
        assert limits.max_memory_mb is not None
        assert limits.max_iterations is not None
        assert limits.check_interval is not None
        assert limits.gc_threshold is not None

    def test_from_config_fields(self):
        """Test that ResourceLimits has expected fields."""
        limits = ResourceLimits.from_config()
        assert hasattr(limits, "max_memory_mb")
        assert hasattr(limits, "max_iterations")
        assert hasattr(limits, "check_interval")
        assert hasattr(limits, "gc_threshold")
        assert hasattr(limits, "memory_usage_ratio")
        assert hasattr(limits, "max_safe_memory_mb")
        assert hasattr(limits, "min_free_memory_mb")
        assert hasattr(limits, "memory_check_aggressive")
        assert hasattr(limits, "gc_auto_collect")
        assert hasattr(limits, "gc_frequency")
        assert hasattr(limits, "gc_force_on_limit")
        assert hasattr(limits, "gc_threshold_ratio")


class TestResourceMonitor:
    """Test ResourceMonitor functionality."""

    def setup_method(self):
        """Setup test method."""
        self.limits = ResourceLimits.from_config()
        self.monitor = ResourceMonitor(self.limits)

    def test_init(self):
        """Test ResourceMonitor initialization."""
        monitor = ResourceMonitor()
        assert monitor.limits is not None
        assert monitor.monitoring is False
        assert monitor.monitor_thread is None

    def test_init_with_limits(self):
        """Test ResourceMonitor initialization with custom limits."""
        limits = ResourceLimits.from_config()
        monitor = ResourceMonitor(limits)
        assert monitor.limits == limits
        assert monitor.monitoring is False

    def test_set_violation_callback(self):
        """Test setting violation callback."""
        callback = MagicMock()
        self.monitor.set_violation_callback(callback)
        assert self.monitor.violation_callback == callback

    def test_start_monitoring(self):
        """Test starting monitoring."""
        self.monitor.start_monitoring()
        assert self.monitor.monitoring is True

    def test_stop_monitoring(self):
        """Test stopping monitoring."""
        self.monitor.start_monitoring()
        self.monitor.stop_monitoring()
        assert self.monitor.monitoring is False


class TestResourceFunctions:
    """Test standalone resource functions."""

    def test_get_available_memory_mb(self):
        """Test getting available memory."""
        mem_mb = get_available_memory_mb()
        assert isinstance(mem_mb, int | float)
        assert mem_mb > 0

    def test_get_safe_memory_limit(self):
        """Test getting safe memory limit."""
        safe_limit = get_safe_memory_limit()
        assert isinstance(safe_limit, int | float)
        assert safe_limit > 0

    def test_estimate_algorithm_memory(self):
        """Test estimating algorithm memory usage."""
        n, L = 10, 100
        estimated_mb = estimate_algorithm_memory(n, L, "baseline")
        assert isinstance(estimated_mb, int | float)
        assert estimated_mb > 0

    def test_check_algorithm_feasibility(self):
        """Test checking algorithm feasibility."""
        n, L = 10, 100
        is_feasible, reason = check_algorithm_feasibility(n, L, "baseline")

        assert isinstance(is_feasible, bool)
        assert isinstance(reason, str)

    def test_force_garbage_collection(self):
        """Test forcing garbage collection."""
        result = force_garbage_collection()
        # Function might return None or memory info
        # Just check that it doesn't crash
        assert True
