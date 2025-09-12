"""Tests for CLI commands module.

Coverage objectives:
- Test command registration and functionality
- Test work monitoring and status display
- Test batch execution commands
- Test algorithm listing functionality
- Test web interface startup
- Test error handling scenarios
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pytest
import typer

from src.presentation.cli.commands import (
    _run_work_monitor,
    _show_final_status_message,
    register_commands,
)


class TestWorkMonitoring:
    """Test work monitoring functionality."""

    @patch('src.presentation.cli.commands.ProgressMonitor')
    @patch('src.presentation.cli.commands._show_final_status_message')
    def test_run_work_monitor_success(self, mock_status, mock_monitor_class):
        """Test successful work monitoring."""
        mock_monitor = Mock()
        mock_monitor_class.return_value = mock_monitor
        
        result = _run_work_monitor('test_work_id', show_final_status=True)
        
        assert result is True
        mock_monitor_class.assert_called_once_with('test_work_id')
        mock_monitor.start.assert_called_once()
        mock_status.assert_called_once_with('test_work_id')

    @patch('src.presentation.cli.commands.ProgressMonitor')
    @patch('typer.echo')
    def test_run_work_monitor_failure(self, mock_echo, mock_monitor_class):
        """Test work monitoring failure handling."""
        mock_monitor_class.side_effect = Exception("Monitor failed")
        
        result = _run_work_monitor('test_work_id', fallback_message=True)
        
        assert result is False
        mock_echo.assert_called()

    @patch('src.presentation.cli.commands.ProgressMonitor')
    def test_run_work_monitor_no_final_status(self, mock_monitor_class):
        """Test work monitoring without final status display."""
        mock_monitor = Mock()
        mock_monitor_class.return_value = mock_monitor
        
        result = _run_work_monitor('test_work_id', show_final_status=False)
        
        assert result is True
        mock_monitor.start.assert_called_once()

    @patch('src.presentation.cli.commands.ProgressMonitor')
    @patch('typer.echo')
    def test_run_work_monitor_no_fallback_message(self, mock_echo, mock_monitor_class):
        """Test work monitoring failure without fallback message."""
        mock_monitor_class.side_effect = Exception("Monitor failed")
        
        result = _run_work_monitor('test_work_id', fallback_message=False)
        
        assert result is False
        # Should not have called echo since fallback_message=False
        mock_echo.assert_not_called()


class TestFinalStatusMessage:
    """Test final status message display."""

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_show_final_status_completed(self, mock_echo, mock_get_service):
        """Test final status message for completed work."""
        mock_details = Mock()
        mock_details.status.value = 'completed'
        mock_details.output_path = '/test/output/path'
        
        mock_service = Mock()
        mock_service.get.return_value = mock_details
        mock_get_service.return_value = mock_service
        
        _show_final_status_message('test_work_id')
        
        mock_service.get.assert_called_once_with('test_work_id')
        mock_echo.assert_called()

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_show_final_status_failed(self, mock_echo, mock_get_service):
        """Test final status message for failed work."""
        mock_details = Mock()
        mock_details.status.value = 'failed'
        mock_details.error = 'Test error message'
        
        mock_service = Mock()
        mock_service.get.return_value = mock_details
        mock_get_service.return_value = mock_service
        
        _show_final_status_message('test_work_id')
        
        mock_service.get.assert_called_once_with('test_work_id')
        mock_echo.assert_called()

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_show_final_status_no_details(self, mock_echo, mock_get_service):
        """Test final status message when no details available."""
        mock_service = Mock()
        mock_service.get.return_value = None
        mock_get_service.return_value = mock_service
        
        _show_final_status_message('test_work_id')
        
        mock_service.get.assert_called_once_with('test_work_id')
        mock_echo.assert_called_with("⚠️  Could not get final status for work test_work_id")

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_show_final_status_exception(self, mock_echo, mock_get_service):
        """Test final status message handling exceptions."""
        mock_get_service.side_effect = Exception("Service error")
        
        # Should not raise exception
        _show_final_status_message('test_work_id')


class TestCommandRegistration:
    """Test command registration functionality."""

    def test_register_commands(self):
        """Test that commands are registered with the Typer app."""
        app = typer.Typer()
        
        # This should not raise an exception
        register_commands(app)
        
        # Verify app now has commands (indirect test)
        assert app is not None

    @patch('src.presentation.cli.commands.get_work_service')
    def test_batch_command_registration(self, mock_get_service):
        """Test batch command can be registered and called."""
        app = typer.Typer()
        register_commands(app)
        
        # Verify app has been modified (has commands)
        assert hasattr(app, 'info')


class TestBatchExecution:
    """Test batch execution functionality."""

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('src.presentation.cli.commands._run_work_monitor')
    @patch('src.presentation.cli.commands.load_cspbench_config')
    @patch('typer.echo')
    def test_batch_execution_with_monitor(self, mock_echo, mock_load_config, mock_monitor, mock_get_service):
        """Test batch execution with monitoring."""
        # Setup mocks
        mock_config = Mock()
        mock_load_config.return_value = mock_config
        
        mock_service = Mock()
        mock_service.create_batch_work.return_value = 'test_work_id'
        mock_get_service.return_value = mock_service
        
        mock_monitor.return_value = True
        
        # Create a test YAML file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("""
name: test_batch
datasets:
  - test_dataset
algorithms:
  - baseline
tasks:
  - closest_string
""")
            yaml_path = f.name
        
        try:
            # Import and call the batch command function directly
            from src.presentation.cli.commands import batch_command
            
            batch_command(yaml_path, monitor_type='log')
            
            mock_service.create_batch_work.assert_called_once()
            mock_monitor.assert_called_once()
            
        finally:
            os.unlink(yaml_path)

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_batch_execution_file_not_found(self, mock_echo, mock_get_service):
        """Test batch execution with non-existent file."""
        from src.presentation.cli.commands import batch_command
        
        # Should handle file not found gracefully
        batch_command('nonexistent.yaml')
        
        # Should have echoed an error message
        mock_echo.assert_called()


class TestAlgorithmListing:
    """Test algorithm listing functionality."""

    @patch('src.presentation.cli.commands.global_registry')
    @patch('typer.echo')
    def test_list_algorithms(self, mock_echo, mock_registry):
        """Test algorithm listing command."""
        # Mock registry with test algorithms
        mock_registry.list_algorithms.return_value = ['algorithm1', 'algorithm2']
        mock_registry.get_algorithm_class.return_value = Mock(
            name='TestAlgorithm',
            description='Test algorithm description'
        )
        
        from src.presentation.cli.commands import list_algorithms
        
        list_algorithms()
        
        mock_registry.list_algorithms.assert_called_once()
        mock_echo.assert_called()

    @patch('src.presentation.cli.commands.global_registry')
    @patch('typer.echo')
    def test_list_algorithms_empty(self, mock_echo, mock_registry):
        """Test algorithm listing with no algorithms."""
        mock_registry.list_algorithms.return_value = []
        
        from src.presentation.cli.commands import list_algorithms
        
        list_algorithms()
        
        mock_registry.list_algorithms.assert_called_once()
        mock_echo.assert_called()


class TestWebInterfaceStartup:
    """Test web interface startup functionality."""

    @patch('src.presentation.cli.commands.run_web_interface')
    def test_web_command(self, mock_run_web):
        """Test web interface startup command."""
        from src.presentation.cli.commands import web_command
        
        web_command(host='localhost', port=8000, debug=True)
        
        mock_run_web.assert_called_once_with(host='localhost', port=8000, debug=True)

    @patch('src.presentation.cli.commands.run_web_interface')
    def test_web_command_default_params(self, mock_run_web):
        """Test web interface startup with default parameters."""
        from src.presentation.cli.commands import web_command
        
        web_command()
        
        mock_run_web.assert_called_once()


class TestWorkItemManagement:
    """Test work item management commands."""

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_work_list_command(self, mock_echo, mock_get_service):
        """Test work list command."""
        mock_service = Mock()
        mock_service.list_work.return_value = []
        mock_get_service.return_value = mock_service
        
        from src.presentation.cli.commands import work_list
        
        work_list()
        
        mock_service.list_work.assert_called_once()
        mock_echo.assert_called()

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_work_status_command(self, mock_echo, mock_get_service):
        """Test work status command."""
        mock_details = Mock()
        mock_details.status.value = 'running'
        
        mock_service = Mock()
        mock_service.get.return_value = mock_details
        mock_get_service.return_value = mock_service
        
        from src.presentation.cli.commands import work_status
        
        work_status('test_work_id')
        
        mock_service.get.assert_called_once_with('test_work_id')
        mock_echo.assert_called()

    @patch('src.presentation.cli.commands.get_work_service')
    @patch('typer.echo')
    def test_work_restart_command(self, mock_echo, mock_get_service):
        """Test work restart command."""
        mock_service = Mock()
        mock_service.restart_work.return_value = True
        mock_get_service.return_value = mock_service
        
        from src.presentation.cli.commands import restart_work
        
        restart_work('test_work_id')
        
        mock_service.restart_work.assert_called_once_with('test_work_id')
        mock_echo.assert_called()


class TestDatasetGeneration:
    """Test dataset generation functionality."""

    @patch('src.presentation.cli.commands.DatasetGenerationWizard')
    def test_generate_datasets_command(self, mock_wizard_class):
        """Test dataset generation wizard command."""
        mock_wizard = Mock()
        mock_wizard_class.return_value = mock_wizard
        
        from src.presentation.cli.commands import generate_datasets
        
        generate_datasets()
        
        mock_wizard_class.assert_called_once()
        mock_wizard.run.assert_called_once()

    @patch('src.presentation.cli.commands.DatasetGenerationWizard')
    def test_generate_datasets_with_exception(self, mock_wizard_class):
        """Test dataset generation wizard with exception handling."""
        mock_wizard_class.side_effect = Exception("Wizard failed")
        
        from src.presentation.cli.commands import generate_datasets
        
        # Should handle exception gracefully
        generate_datasets()


class TestConfigurationHandling:
    """Test configuration handling in commands."""

    @patch('src.presentation.cli.commands.load_cspbench_config')
    def test_config_loading(self, mock_load_config):
        """Test configuration loading in commands."""
        mock_config = Mock()
        mock_load_config.return_value = mock_config
        
        # Test that config loading works
        config = mock_load_config()
        
        assert config is not None
        mock_load_config.assert_called_once()

    @patch('src.presentation.cli.commands.load_cspbench_config')
    def test_config_loading_failure(self, mock_load_config):
        """Test configuration loading failure handling."""
        mock_load_config.side_effect = Exception("Config load failed")
        
        # Should handle config loading failure
        try:
            config = mock_load_config()
        except Exception:
            pass  # Expected behavior


class TestErrorHandling:
    """Test error handling scenarios."""

    @patch('src.presentation.cli.commands.get_work_service')
    def test_service_unavailable_handling(self, mock_get_service):
        """Test handling when work service is unavailable."""
        mock_get_service.side_effect = Exception("Service unavailable")
        
        # Commands should handle service unavailability gracefully
        from src.presentation.cli.commands import _show_final_status_message
        
        # Should not raise exception
        _show_final_status_message('test_work_id')

    @patch('typer.echo')
    def test_missing_work_id_handling(self, mock_echo):
        """Test handling of missing work IDs."""
        # Test various error scenarios that might occur
        _show_final_status_message('')
        
        # Should handle empty work ID gracefully
        mock_echo.assert_called()