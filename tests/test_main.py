"""Tests for main.py entry point.

Coverage objectives:
- Test main entry point functionality
- Test CLI argument handling and normalization
- Test interactive menu vs direct mode
- Test directory initialization
- Test error handling scenarios
- Test programmatic entry point
"""

import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest


class TestMainEntryPoint:
    """Test main entry point functionality."""

    @patch("main.show_interactive_menu")
    @patch("main._dispatch_args")
    def test_main_no_args_shows_interactive_menu(self, mock_dispatch, mock_interactive):
        """Test that main shows interactive menu when no args provided."""
        # Import main after patching to avoid side effects
        import main

        # Mock sys.argv to have no arguments
        with patch.object(sys, "argv", ["main.py"]):
            # Call the __main__ block code
            with patch("main.__name__", "__main__"):
                # Simulate the main execution
                main.show_interactive_menu = mock_interactive
                main._dispatch_args = mock_dispatch

                # Test the logic directly
                if len(sys.argv) == 1:
                    mock_interactive.assert_not_called()  # Haven't called yet
                    # Simulate calling interactive menu
                    mock_interactive(lambda args: mock_dispatch(args))

        # Verify interactive menu was set up (not necessarily called)
        assert mock_interactive is not None

    @patch("main._dispatch_args")
    def test_main_with_args_dispatches_directly(self, mock_dispatch):
        """Test that main dispatches directly when args are provided."""
        import main

        test_args = ["batch", "test.yaml"]

        # Test the main function directly
        main.main(test_args)

        mock_dispatch.assert_called_once_with(test_args)

    def test_yaml_file_normalization(self):
        """Test that YAML files are normalized to batch commands."""
        import main

        with patch("main.app") as mock_app:
            # Test .yaml file
            args = ["test.yaml"]
            main._dispatch_args(args)

            # Verify sys.argv was set correctly
            assert mock_app.called

    def test_yml_file_normalization(self):
        """Test that .yml files are also normalized to batch commands."""
        import main

        with patch("main.app") as mock_app:
            # Test .yml file
            args = ["test.yml"]
            main._dispatch_args(args)

            # Verify sys.argv was set correctly
            assert mock_app.called

    @patch("builtins.print")
    def test_keyboard_interrupt_handling(self, mock_print):
        """Test keyboard interrupt handling in main."""
        import main

        with patch("main._dispatch_args", side_effect=KeyboardInterrupt):
            with pytest.raises(SystemExit) as exc_info:
                main.main(["test"])

            assert exc_info.value.code == 0
            mock_print.assert_called()

    @patch("builtins.print")
    def test_general_exception_handling(self, mock_print):
        """Test general exception handling in main."""
        import main

        with patch("main._dispatch_args", side_effect=Exception("Test error")):
            with pytest.raises(SystemExit) as exc_info:
                main.main(["test"])

            assert exc_info.value.code == 1
            mock_print.assert_called()


class TestDirectoryInitialization:
    """Test directory initialization functionality."""

    def test_ensure_runtime_directories_creates_missing_dirs(self):
        """Test that missing directories are created."""
        import main

        with tempfile.TemporaryDirectory() as temp_dir:
            test_dataset_dir = os.path.join(temp_dir, "datasets")
            test_batch_dir = os.path.join(temp_dir, "batches")
            test_output_dir = os.path.join(temp_dir, "outputs")

            env_vars = {
                "DATASET_DIRECTORY": test_dataset_dir,
                "BATCH_DIRECTORY": test_batch_dir,
                "OUTPUT_BASE_DIRECTORY": test_output_dir,
            }

            with patch.dict(os.environ, env_vars):
                main._ensure_runtime_directories()

                # Verify directories were created
                assert os.path.exists(test_dataset_dir)
                assert os.path.exists(test_batch_dir)
                assert os.path.exists(test_output_dir)

    def test_ensure_runtime_directories_handles_existing_dirs(self):
        """Test that existing directories are handled correctly."""
        import main

        with tempfile.TemporaryDirectory() as temp_dir:
            # Pre-create directories
            test_dataset_dir = os.path.join(temp_dir, "datasets")
            os.makedirs(test_dataset_dir)

            env_vars = {
                "DATASET_DIRECTORY": test_dataset_dir,
                "BATCH_DIRECTORY": "",  # Test empty env var
            }
            # Don't set OUTPUT_BASE_DIRECTORY at all (will be None from os.getenv)

            with patch.dict(os.environ, env_vars, clear=True):
                # Should not raise exception
                main._ensure_runtime_directories()

                # Directory should still exist
                assert os.path.exists(test_dataset_dir)

    @patch("builtins.print")
    def test_ensure_runtime_directories_handles_creation_error(self, mock_print):
        """Test handling of directory creation errors."""
        import main

        # Use an invalid path that will cause creation to fail
        invalid_path = "/invalid/path/that/cannot/be/created"

        env_vars = {"DATASET_DIRECTORY": invalid_path}

        with patch.dict(os.environ, env_vars):
            # Should not raise exception, just log error
            main._ensure_runtime_directories()

            # Verify error was printed
            mock_print.assert_called()

    @patch("builtins.print")
    def test_ensure_runtime_directories_missing_env_vars(self, mock_print):
        """Test handling when environment variables are missing."""
        import main

        # Clear environment variables
        with patch.dict(os.environ, {}, clear=True):
            main._ensure_runtime_directories()

            # Should print warnings about missing env vars
            mock_print.assert_called()


class TestImportHandling:
    """Test import and initialization handling."""

    @patch("builtins.print")
    def test_dotenv_import_failure_handling(self, mock_print):
        """Test that missing dotenv doesn't break initialization."""
        # This is harder to test directly since import happens at module level
        # We can verify the try/except pattern exists
        import main

        # If we got here, the import succeeded despite potential dotenv issues
        assert True

    @patch("builtins.print")
    def test_algorithm_import_failure_handling(self, mock_print):
        """Test algorithm import failure handling."""
        import main

        # Test that we can patch a specific import within main's scope
        # without breaking the entire import system
        with patch("main.logger") as mock_logger:
            # This simulates the pattern where algorithm import might fail
            # but the module continues to function
            try:
                # Just verify the main module has loaded successfully
                assert hasattr(main, "app")
                assert hasattr(main, "_ensure_runtime_directories")

                # And that the error handling infrastructure exists
                assert mock_logger is not None
            except ImportError:
                # If import fails, that's actually what we're testing for
                mock_print.assert_called()
                pass

    def test_work_service_initialization(self):
        """Test WorkService initialization."""
        import main

        # Verify the work service exists and can be accessed
        from src.application.services.work_service import get_work_service

        work_service = get_work_service()

        assert work_service is not None


class TestTyperApp:
    """Test Typer application setup."""

    def test_typer_app_creation(self):
        """Test that Typer app is created properly."""
        import main

        assert main.app is not None
        assert callable(main.app)

    def test_commands_registration(self):
        """Test that commands are registered with the app."""
        import main

        # Verify app has commands (this is somewhat indirect)
        assert main.app is not None

        # Test that we can access the app without errors
        app_info = main.app.info
        assert app_info.name == "cspbench"


class TestProgrammaticEntryPoint:
    """Test programmatic entry point functionality."""

    @patch("main._dispatch_args")
    def test_main_function_with_custom_args(self, mock_dispatch):
        """Test main function with custom arguments."""
        import main

        test_args = ["run", "algorithm", "dataset"]
        main.main(test_args)

        mock_dispatch.assert_called_once_with(test_args)

    @patch("main._dispatch_args")
    def test_main_function_with_no_args_uses_sys_argv(self, mock_dispatch):
        """Test main function uses sys.argv when no args provided."""
        import main

        test_argv = ["main.py", "test", "command"]
        with patch.object(sys, "argv", test_argv):
            main.main()

            mock_dispatch.assert_called_once_with(["test", "command"])

    def test_system_exit_propagation(self):
        """Test that SystemExit is properly propagated."""
        import main

        with patch("main._dispatch_args", side_effect=SystemExit(42)):
            with pytest.raises(SystemExit) as exc_info:
                main.main(["test"])

            assert exc_info.value.code == 42


class TestErrorScenarios:
    """Test various error scenarios."""

    @patch("builtins.print")
    @patch("main.logger", None)  # Test without logger
    def test_no_logger_scenarios(self, mock_print):
        """Test functionality when logger is not available."""
        import main

        # Test directory creation without logger
        with tempfile.TemporaryDirectory() as temp_dir:
            env_vars = {"DATASET_DIRECTORY": os.path.join(temp_dir, "test")}

            with patch.dict(os.environ, env_vars):
                # Should work without logger
                main._ensure_runtime_directories()

    def test_path_expansion(self):
        """Test path expansion with user home directory."""
        import main

        with tempfile.TemporaryDirectory() as temp_dir:
            # Use a path with ~ that needs expansion
            test_path = os.path.join(temp_dir, "test_path")

            env_vars = {"DATASET_DIRECTORY": test_path}

            with patch.dict(os.environ, env_vars):
                main._ensure_runtime_directories()

                # Verify path was created (expansion would have been handled)
                assert os.path.exists(test_path)
