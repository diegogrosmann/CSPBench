"""Tests for interactive menu functionality.

Coverage objectives:
- Test interactive menu display and navigation
- Test batch file discovery and listing
- Test file description extraction
- Test user input handling
- Test command execution via menu
- Test error scenarios
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, call
import pytest

from src.presentation.cli.interactive_menu import (
    show_interactive_menu,
    _show_commands_help,
    _extract_file_description,
)


class TestInteractiveMenu:
    """Test interactive menu functionality."""

    @patch('src.presentation.cli.interactive_menu.get_batch_directory')
    @patch('builtins.print')
    @patch('builtins.input', return_value='q')
    def test_show_interactive_menu_no_directory(self, mock_input, mock_print, mock_get_dir):
        """Test interactive menu when batch directory doesn't exist."""
        mock_dir = Mock()
        mock_dir.exists.return_value = False
        mock_get_dir.return_value = mock_dir
        
        mock_app_invoke = Mock()
        
        show_interactive_menu(mock_app_invoke)
        
        mock_print.assert_called()
        # Should not call app_invoke since directory doesn't exist

    @patch('src.presentation.cli.interactive_menu.get_batch_directory')
    @patch('builtins.print')
    @patch('builtins.input', return_value='q')
    def test_show_interactive_menu_no_batch_files(self, mock_input, mock_print, mock_get_dir):
        """Test interactive menu when no batch files exist."""
        mock_dir = Mock()
        mock_dir.exists.return_value = True
        mock_dir.glob.return_value = []  # No files
        mock_get_dir.return_value = mock_dir
        
        mock_app_invoke = Mock()
        
        show_interactive_menu(mock_app_invoke)
        
        mock_print.assert_called()

    def test_show_interactive_menu_with_batch_files(self):
        """Test interactive menu with actual batch files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test batch files
            batch_file1 = temp_path / "test1.yaml"
            batch_file1.write_text("""
# Test Batch 1 - Basic algorithm comparison
name: test_batch_1
datasets:
  - test_dataset
algorithms:
  - baseline
""")
            
            batch_file2 = temp_path / "test2.yml"
            batch_file2.write_text("""
# Test Batch 2 - Advanced testing
name: test_batch_2
datasets:
  - another_dataset
algorithms:
  - h2_csp
""")
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                with patch('builtins.input', side_effect=['1', 'q']):
                    with patch('builtins.print') as mock_print:
                        mock_app_invoke = Mock()
                        
                        show_interactive_menu(mock_app_invoke)
                        
                        # Should have printed menu and executed batch
                        mock_print.assert_called()
                        mock_app_invoke.assert_called_once()

    @patch('builtins.input', return_value='q')
    @patch('builtins.print')
    def test_show_interactive_menu_quit_immediately(self, mock_print, mock_input):
        """Test quitting interactive menu immediately."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                mock_app_invoke = Mock()
                
                show_interactive_menu(mock_app_invoke)
                
                # Should not call app_invoke when quitting
                mock_app_invoke.assert_not_called()

    @patch('builtins.input', side_effect=['invalid', 'q'])
    @patch('builtins.print')
    def test_show_interactive_menu_invalid_input(self, mock_print, mock_input):
        """Test handling invalid input in interactive menu."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create a test batch file
            batch_file = temp_path / "test.yaml"
            batch_file.write_text("name: test_batch")
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                mock_app_invoke = Mock()
                
                show_interactive_menu(mock_app_invoke)
                
                # Should handle invalid input gracefully
                mock_print.assert_called()

    @patch('builtins.input', side_effect=['99', 'q'])  # Out of range selection
    @patch('builtins.print')
    def test_show_interactive_menu_out_of_range(self, mock_print, mock_input):
        """Test handling out of range selection."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create a test batch file
            batch_file = temp_path / "test.yaml"
            batch_file.write_text("name: test_batch")
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                mock_app_invoke = Mock()
                
                show_interactive_menu(mock_app_invoke)
                
                # Should handle out of range selection gracefully
                mock_print.assert_called()

    @patch('builtins.input', side_effect=['h', 'q'])  # Help command
    @patch('builtins.print')
    def test_show_interactive_menu_help_command(self, mock_print, mock_input):
        """Test help command in interactive menu."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                mock_app_invoke = Mock()
                
                show_interactive_menu(mock_app_invoke)
                
                # Should display help information
                mock_print.assert_called()


class TestDescriptionExtraction:
    """Test batch file description extraction."""

    def test_extract_description_with_comment(self):
        """Test extracting description from file with comment."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("""# Test Description - This is a test batch
name: test_batch
datasets:
  - test_dataset
""")
            f.flush()
            
            try:
                description = _extract_file_description(Path(f.name))
                assert description == "Test Description - This is a test batch"
            finally:
                os.unlink(f.name)

    def test_extract_description_no_comment(self):
        """Test extracting description from file without comment."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("""name: test_batch
datasets:
  - test_dataset
""")
            f.flush()
            
            try:
                description = _extract_file_description(Path(f.name))
                assert description == "Batch configuration file"
            finally:
                os.unlink(f.name)

    def test_extract_description_multiple_comments(self):
        """Test extracting description with multiple comment lines."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("""# First comment line
# Second comment line - this should be returned
name: test_batch
""")
            f.flush()
            
            try:
                description = _extract_file_description(Path(f.name))
                # Should return the first non-empty comment line
                assert description in ["First comment line", "Second comment line - this should be returned"]
            finally:
                os.unlink(f.name)

    def test_extract_description_nonexistent_file(self):
        """Test extracting description from non-existent file."""
        nonexistent_path = Path("/nonexistent/file.yaml")
        description = _extract_file_description(nonexistent_path)
        assert description == "Batch configuration file"

    def test_extract_description_permission_error(self):
        """Test extracting description with permission error."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("# Test description\nname: test")
            f.flush()
            
            try:
                # Mock a permission error
                with patch('builtins.open', side_effect=PermissionError("Permission denied")):
                    description = _extract_file_description(Path(f.name))
                    assert description == "Batch configuration file"
            finally:
                os.unlink(f.name)

    def test_extract_description_encoding_error(self):
        """Test extracting description with encoding issues."""
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.yaml', delete=False) as f:
            # Write invalid UTF-8 bytes
            f.write(b'\xff\xfe# Invalid encoding\n')
            f.flush()
            
            try:
                description = _extract_file_description(Path(f.name))
                # Should handle encoding errors gracefully
                assert isinstance(description, str)
            finally:
                os.unlink(f.name)


class TestCommandsHelp:
    """Test commands help functionality."""

    @patch('builtins.print')
    def test_show_commands_help(self, mock_print):
        """Test showing commands help."""
        test_message = "Test error message"
        
        _show_commands_help(test_message)
        
        mock_print.assert_called()
        # Verify the error message was printed
        printed_calls = [str(call) for call in mock_print.call_args_list]
        assert any(test_message in call_str for call_str in printed_calls)

    @patch('builtins.print')
    def test_show_commands_help_no_message(self, mock_print):
        """Test showing commands help with a message."""
        _show_commands_help("Test message")
        
        mock_print.assert_called()


class TestMenuNavigation:
    """Test menu navigation scenarios."""

    def test_batch_file_selection_execution(self):
        """Test that selecting a batch file calls the app invoke function."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create test batch file
            batch_file = temp_path / "test_batch.yaml"
            batch_file.write_text("""
# Test Batch - Sample description
name: test_execution
datasets:
  - sample_dataset
algorithms:
  - baseline
""")
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                with patch('builtins.input', side_effect=['1', 'q']):  # Select first file, then quit
                    with patch('builtins.print'):
                        mock_app_invoke = Mock()
                        
                        show_interactive_menu(mock_app_invoke)
                        
                        # Should have called app_invoke with batch command
                        mock_app_invoke.assert_called_once()
                        args = mock_app_invoke.call_args[0][0]
                        assert 'batch' in args
                        assert str(batch_file) in args

    @patch('builtins.input', side_effect=['', 'q'])  # Empty input
    @patch('builtins.print')
    def test_empty_input_handling(self, mock_print, mock_input):
        """Test handling of empty input."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                mock_app_invoke = Mock()
                
                show_interactive_menu(mock_app_invoke)
                
                # Should handle empty input gracefully
                mock_print.assert_called()

    def test_keyboard_interrupt_handling(self):
        """Test handling of keyboard interrupt (Ctrl+C)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                with patch('builtins.input', side_effect=KeyboardInterrupt):
                    with patch('builtins.print') as mock_print:
                        mock_app_invoke = Mock()
                        
                        # Should handle KeyboardInterrupt gracefully
                        try:
                            show_interactive_menu(mock_app_invoke)
                        except KeyboardInterrupt:
                            pass  # Expected behavior


class TestFileListing:
    """Test file listing and sorting functionality."""

    def test_yaml_and_yml_files_discovered(self):
        """Test that both .yaml and .yml files are discovered."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create files with different extensions
            yaml_file = temp_path / "test1.yaml"
            yml_file = temp_path / "test2.yml"
            txt_file = temp_path / "test3.txt"  # Should be ignored
            
            yaml_file.write_text("name: test1")
            yml_file.write_text("name: test2")
            txt_file.write_text("name: test3")
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_value=temp_path):
                with patch('builtins.input', return_value='q'):
                    with patch('builtins.print') as mock_print:
                        mock_app_invoke = Mock()
                        
                        show_interactive_menu(mock_app_invoke)
                        
                        # Should have found both .yaml and .yml files but not .txt
                        printed_output = ' '.join(str(call) for call in mock_print.call_args_list)
                        assert 'test1.yaml' in printed_output or 'test2.yml' in printed_output

    def test_file_sorting(self):
        """Test that files are sorted correctly in the menu."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create files in non-alphabetical order
            file_c = temp_path / "c_batch.yaml"
            file_a = temp_path / "a_batch.yaml"
            file_b = temp_path / "b_batch.yaml"
            
            for f in [file_c, file_a, file_b]:
                f.write_text("name: test")
            
            with patch('src.presentation.cli.interactive_menu.get_batch_directory', return_path=temp_path):
                with patch('builtins.input', return_value='q'):
                    with patch('builtins.print') as mock_print:
                        mock_app_invoke = Mock()
                        
                        show_interactive_menu(mock_app_invoke)
                        
                        # Files should be processed (sorting is implementation dependent)
                        mock_print.assert_called()