"""
Basic tests for ExecutionManager to validate unification approach.
"""

import pytest
from unittest.mock import Mock, patch
from pathlib import Path

from src.application.services.execution_manager import ExecutionManager, ProgressReporter
from src.domain.config import CSPBenchConfig, MetadataConfig, TasksGroupConfig
from src.infrastructure.monitoring.monitor_interface import NoOpMonitor


class TestExecutionManager:
    """Test cases for ExecutionManager unified execution."""
    
    def test_execution_manager_creation(self):
        """Test that ExecutionManager can be created."""
        manager = ExecutionManager()
        assert manager is not None
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_sync_execution_mode(self, mock_get_work_service):
        """Test synchronous execution mode."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work123", {"output_path": "/tmp/test"})
        mock_wm.mark_running.return_value = True
        mock_wm.mark_finished.return_value = True
        mock_wm.get.return_value = {
            "status": "completed", 
            "output_path": "/tmp/test",
            "error": None,
            "summary": {}
        }
        mock_get_work_service.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=MetadataConfig(
                name="test_config", 
                description="test",
                author="test",
                version="1.0",
                creation_date="2025-01-01"
            ),
            tasks=TasksGroupConfig()
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Mock the pipeline execution
        with patch.object(manager, '_execute_pipeline_direct') as mock_execute:
            mock_execute.return_value = {"status": "completed", "results": {}}
            
            # Execute in sync mode
            work_id, results = manager.execute(
                config=config,
                monitor=monitor,
                mode="sync",
                extra={"test": "data"}
            )
            
            # Verify results
            assert work_id == "work123"
            assert results["status"] == "completed"
            assert results["output_path"] == "/tmp/test"
            
            # Verify WorkManager calls
            mock_wm.submit.assert_called_once()
            mock_wm.mark_running.assert_called_once_with("work123")
            mock_wm.mark_finished.assert_called_once_with("work123")
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_async_execution_mode(self, mock_get_work_service):
        """Test asynchronous execution mode."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work456", {"output_path": "/tmp/async"})
        mock_get_work_service.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=MetadataConfig(
                name="async_config", 
                description="test",
                author="test", 
                version="1.0",
                creation_date="2025-01-01"
            ),
            tasks=TasksGroupConfig()
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Execute in async mode
        work_id = manager.execute(
            config=config,
            monitor=monitor,
            mode="async",
            extra={"test": "async_data"}
        )
        
        # Verify async returns only work_id
        assert work_id == "work456"
        
        # Verify WorkManager submit was called
        mock_wm.submit.assert_called_once()
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_execution_error_handling(self, mock_get_work_service):
        """Test error handling during execution."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work789", {"output_path": "/tmp/error"})
        mock_wm.mark_running.return_value = True
        mock_wm.mark_error.return_value = True
        mock_get_work_service.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=MetadataConfig(
                name="error_config",
                description="test", 
                author="test",
                version="1.0",
                creation_date="2025-01-01"
            ),
            tasks=TasksGroupConfig()
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Mock pipeline execution to raise error
        with patch.object(manager, '_execute_pipeline_direct') as mock_execute:
            mock_execute.side_effect = Exception("Pipeline error")
            
            # Execute and expect exception
            with pytest.raises(Exception) as exc_info:
                manager.execute(
                    config=config,
                    monitor=monitor,
                    mode="sync"
                )
            
            assert "Pipeline error" in str(exc_info.value)
            
            # Verify error marking
            mock_wm.mark_error.assert_called_once_with("work789", "Pipeline error")


class TestProgressReporter:
    """Test cases for ProgressReporter monitor."""
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_progress_reporter_creation(self, mock_get_work_service):
        """Test ProgressReporter creation and basic functionality."""
        mock_wm = Mock()
        mock_get_work_service.return_value = mock_wm
        
        reporter = ProgressReporter("test_work_id")
        assert reporter.work_id == "test_work_id"
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_progress_reporter_on_progress(self, mock_get_work_service):
        """Test progress reporting functionality."""
        mock_wm = Mock()
        mock_get_work_service.return_value = mock_wm
        
        reporter = ProgressReporter("progress_test")
        
        # Test on_progress call
        reporter.on_progress(0.5, {"step": "processing"})
        
        # Verify WorkManager was called with progress
        # Note: This test is basic since we haven't implemented mark_progress yet
        assert reporter.work_id == "progress_test"

import pytest
from unittest.mock import Mock, patch
from pathlib import Path

from src.application.services.execution_manager import ExecutionManager, ProgressReporter
from src.domain.config import CSPBenchConfig, MetadataConfig, TasksGroupConfig
from src.infrastructure.monitoring.monitor_interface import NoOpMonitor


class TestExecutionManager:
    """Test cases for ExecutionManager unified execution."""
    
    def test_execution_manager_creation(self):
        """Test that ExecutionManager can be created."""
        manager = ExecutionManager()
        assert manager is not None
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_sync_execution_mode(self, mock_get_work_service):
        """Test synchronous execution mode."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work123", {"output_path": "/tmp/test"})
        mock_wm.mark_running.return_value = True
        mock_wm.mark_finished.return_value = True
        mock_wm.get.return_value = {
            "status": "completed", 
            "output_path": "/tmp/test",
            "error": None,
            "summary": {}
        }
        mock_get_work_service.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=MetadataConfig(
                name="test_config", 
                description="test",
                author="test",
                version="1.0",
                creation_date="2025-01-01"
            ),
            tasks=TasksGroupConfig()
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Mock the pipeline execution
        with patch.object(manager, '_execute_pipeline_direct') as mock_execute:
            mock_execute.return_value = {"status": "completed", "results": {}}
            
            # Execute in sync mode
            work_id, results = manager.execute(
                config=config,
                monitor=monitor,
                mode="sync",
                extra={"test": "data"}
            )
            
            # Verify results
            assert work_id == "work123"
            assert results["status"] == "completed"
            assert results["output_path"] == "/tmp/test"
            
            # Verify WorkManager calls
            mock_wm.submit.assert_called_once()
            mock_wm.mark_running.assert_called_once_with("work123")
            mock_wm.mark_finished.assert_called_once_with("work123")
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_async_execution_mode(self, mock_get_work_service):
        """Test asynchronous execution mode."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work456", {"output_path": "/tmp/async"})
        mock_get_work_service.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=MetadataConfig(
                name="async_config", 
                description="test",
                author="test", 
                version="1.0",
                creation_date="2025-01-01"
            ),
            tasks=TasksGroupConfig()
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Execute in async mode
        work_id = manager.execute(
            config=config,
            monitor=monitor,
            mode="async",
            extra={"test": "async_data"}
        )
        
        # Verify async returns only work_id
        assert work_id == "work456"
        
        # Verify WorkManager submit was called
        mock_wm.submit.assert_called_once()
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_execution_error_handling(self, mock_get_work_service):
        """Test error handling during execution."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work789", {"output_path": "/tmp/error"})
        mock_wm.mark_running.return_value = True
        mock_wm.mark_error.return_value = True
        mock_get_work_service.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=MetadataConfig(
                name="error_config",
                description="test", 
                author="test",
                version="1.0",
                creation_date="2025-01-01"
            ),
            tasks=TasksGroupConfig()
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Mock pipeline execution to raise error
        with patch.object(manager, '_execute_pipeline_direct') as mock_execute:
            mock_execute.side_effect = Exception("Pipeline error")
            
            # Execute and expect exception
            with pytest.raises(Exception) as exc_info:
                manager.execute(
                    config=config,
                    monitor=monitor,
                    mode="sync"
                )
            
            assert "Pipeline error" in str(exc_info.value)
            
            # Verify error marking
            mock_wm.mark_error.assert_called_once_with("work789", "Pipeline error")


class TestProgressReporter:
    """Test cases for ProgressReporter monitor."""
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_progress_reporter_creation(self, mock_get_work_service):
        """Test ProgressReporter creation and basic functionality."""
        mock_wm = Mock()
        mock_get_work_service.return_value = mock_wm
        
        reporter = ProgressReporter("test_work_id")
        assert reporter.work_id == "test_work_id"
    
    @patch('src.application.services.execution_manager.get_work_service')
    def test_progress_reporter_on_progress(self, mock_get_work_service):
        """Test progress reporting functionality."""
        mock_wm = Mock()
        mock_get_work_service.return_value = mock_wm
        
        reporter = ProgressReporter("progress_test")
        
        # Test on_progress call
        reporter.on_progress(0.5, {"step": "processing"})
        
        # Verify WorkManager was called with progress
        # Note: This test is basic since we haven't implemented mark_progress yet
        assert reporter.work_id == "progress_test"

import pytest
from unittest.mock import Mock, patch
from pathlib import Path

from src.application.services.execution_manager import ExecutionManager, ProgressReporter
from src.domain.config import CSPBenchConfig, MetadataConfig, TasksGroupConfig
from src.infrastructure.monitoring.monitor_interface import NoOpMonitor


class TestExecutionManager:
    """Test cases for ExecutionManager unified execution."""
    
    def test_execution_manager_creation(self):
        """Test that ExecutionManager can be created."""
        manager = ExecutionManager()
        assert manager is not None
    
    @patch('src.application.services.execution_manager.get_global_work_manager')
    def test_sync_execution_mode(self, mock_get_work_manager):
        """Test synchronous execution mode."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work123", {"output_path": "/tmp/test"})
        mock_wm.mark_running.return_value = True
        mock_wm.mark_finished.return_value = True
        mock_wm.get.return_value = {
            "status": "completed", 
            "output_path": "/tmp/test",
            "error": None,
            "summary": {}
        }
        mock_get_work_manager.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=CSPBenchMetadata(name="test_config", version="1.0"),
            tasks=TasksGroup(experiments=[], optimizations=[], sensitivities=[])
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Mock the pipeline execution
        with patch.object(manager, '_execute_pipeline_direct') as mock_execute:
            mock_execute.return_value = {"status": "completed", "results": {}}
            
            # Execute in sync mode
            work_id, results = manager.execute(
                config=config,
                monitor=monitor,
                mode="sync",
                extra={"test": "data"}
            )
            
            # Verify results
            assert work_id == "work123"
            assert results["status"] == "completed"
            assert results["output_path"] == "/tmp/test"
            
            # Verify WorkManager calls
            mock_wm.submit.assert_called_once()
            mock_wm.mark_running.assert_called_once_with("work123")
            mock_wm.mark_finished.assert_called_once_with("work123")
    
    @patch('src.application.services.execution_manager.get_global_work_manager')
    def test_async_execution_mode(self, mock_get_work_manager):
        """Test asynchronous execution mode."""
        # Setup mocks
        mock_wm = Mock()
        mock_wm.submit.return_value = ("work456", {"output_path": "/tmp/async"})
        mock_get_work_manager.return_value = mock_wm
        
        # Create test config
        config = CSPBenchConfig(
            metadata=CSPBenchMetadata(name="async_config", version="1.0"),
            tasks=TasksGroup(experiments=[], optimizations=[], sensitivities=[])
        )
        
        monitor = NoOpMonitor()
        manager = ExecutionManager()
        
        # Execute in async mode
        work_id = manager.execute(
            config=config,
            monitor=monitor,
            mode="async",
            extra={"async": "test"}
        )
        
        # Verify immediate return
        assert work_id == "work456"
        
        # Verify submission
        mock_wm.submit.assert_called_once()
        submit_args = mock_wm.submit.call_args
        assert submit_args[1]["config"] == config
        assert submit_args[1]["extra"]["async"] == "test"
        assert submit_args[1]["extra"]["mode"] == "async"


class TestProgressReporter:
    """Test cases for ProgressReporter monitor wrapper."""
    
    def test_progress_reporter_creation(self):
        """Test ProgressReporter creation."""
        mock_monitor = Mock()
        mock_work_service = Mock()
        
        reporter = ProgressReporter("work123", mock_monitor, mock_work_service)
        
        assert reporter.work_id == "work123"
        assert reporter.delegate == mock_monitor
        assert reporter.work_service == mock_work_service
    
    def test_progress_delegation(self):
        """Test that progress calls are delegated to both WorkService and monitor."""
        mock_monitor = Mock()
        mock_work_service = Mock()
        mock_work_service.mark_progress = Mock()
        
        reporter = ProgressReporter("work123", mock_monitor, mock_work_service)
        
        # Call on_progress
        reporter.on_progress(0.5, {"step": "processing"})
        
        # Verify WorkService call
        mock_work_service.mark_progress.assert_called_once_with(
            "work123", 0.5, {"step": "processing"}
        )
        
        # Verify delegate call
        mock_monitor.on_progress.assert_called_once_with(0.5, {"step": "processing"})
    
    def test_should_stop_delegation(self):
        """Test cancellation checking."""
        mock_monitor = Mock()
        mock_monitor.should_stop.return_value = False
        
        mock_work_service = Mock()
        mock_work_service.get_control_flags.return_value = {"cancelled": True}
        
        reporter = ProgressReporter("work123", mock_monitor, mock_work_service)
        
        # Should return True due to WorkService cancellation
        assert reporter.should_stop() is True
        
        # Test fallback to monitor
        mock_work_service.get_control_flags.side_effect = Exception("No control flags")
        mock_monitor.should_stop.return_value = True
        
        assert reporter.should_stop() is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
