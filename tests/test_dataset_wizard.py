"""Tests for dataset wizard functionality.

Coverage objectives:
- Test dataset generation wizard interface
- Test user input collection and validation
- Test synthetic dataset parameter collection
- Test real dataset parameter collection
- Test filename generation
- Test error handling and validation
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
import pytest

from src.presentation.cli.dataset_wizard import (
    DatasetWizard,
    _ask,
    _ask_int,
    _ask_float,
    _ask_optional_int,
)


class TestInputHelpers:
    """Test input helper functions."""

    @patch('builtins.input', return_value='test_value')
    def test_ask_with_input(self, mock_input):
        """Test _ask function with user input."""
        result = _ask("Enter value", "default")
        assert result == "test_value"

    @patch('builtins.input', return_value='')
    def test_ask_with_default(self, mock_input):
        """Test _ask function with default value."""
        result = _ask("Enter value", "default")
        assert result == "default"

    @patch('builtins.input', return_value='')
    def test_ask_no_default(self, mock_input):
        """Test _ask function without default value."""
        result = _ask("Enter value")
        assert result == ""

    @patch('builtins.input', return_value='  whitespace  ')
    def test_ask_strips_whitespace(self, mock_input):
        """Test that _ask strips whitespace from input."""
        result = _ask("Enter value")
        assert result == "whitespace"

    @patch('builtins.input', return_value='42')
    def test_ask_int_valid_input(self, mock_input):
        """Test _ask_int with valid integer."""
        result = _ask_int("Enter number", 10)
        assert result == 42

    @patch('builtins.input', return_value='')
    def test_ask_int_default(self, mock_input):
        """Test _ask_int with default value."""
        result = _ask_int("Enter number", 10)
        assert result == 10

    @patch('builtins.input', side_effect=['invalid', '42'])
    @patch('builtins.print')
    def test_ask_int_invalid_then_valid(self, mock_print, mock_input):
        """Test _ask_int with invalid then valid input."""
        result = _ask_int("Enter number", 10)
        assert result == 42
        assert mock_print.called

    @patch('builtins.input', return_value='3.14')
    def test_ask_float_valid_input(self, mock_input):
        """Test _ask_float with valid float."""
        result = _ask_float("Enter float", 1.0)
        assert result == 3.14

    @patch('builtins.input', return_value='')
    def test_ask_float_default(self, mock_input):
        """Test _ask_float with default value."""
        result = _ask_float("Enter float", 1.0)
        assert result == 1.0

    @patch('builtins.input', side_effect=['invalid', '3.14'])
    @patch('builtins.print')
    def test_ask_float_invalid_then_valid(self, mock_print, mock_input):
        """Test _ask_float with invalid then valid input."""
        result = _ask_float("Enter float", 1.0)
        assert result == 3.14
        assert mock_print.called

    @patch('builtins.input', return_value='42')
    def test_ask_optional_int_valid_input(self, mock_input):
        """Test _ask_optional_int with valid integer."""
        result = _ask_optional_int("Enter number", 10)
        assert result == 42

    @patch('builtins.input', return_value='')
    def test_ask_optional_int_default(self, mock_input):
        """Test _ask_optional_int with default value."""
        result = _ask_optional_int("Enter number", 10)
        assert result == 10

    @patch('builtins.input', return_value='')
    def test_ask_optional_int_none_default(self, mock_input):
        """Test _ask_optional_int with None default."""
        result = _ask_optional_int("Enter number", None)
        assert result is None


class TestDatasetWizard:
    """Test DatasetWizard class functionality."""

    def test_wizard_initialization(self):
        """Test DatasetWizard can be initialized."""
        wizard = DatasetWizard()
        assert wizard is not None

    @patch('builtins.input', return_value='1')
    @patch('builtins.print')
    def test_show_main_menu_synthetic(self, mock_print, mock_input):
        """Test main menu selection for synthetic datasets."""
        wizard = DatasetWizard()
        result = wizard.show_main_menu()
        assert result == "synthetic"

    @patch('builtins.input', return_value='2')
    @patch('builtins.print')
    def test_show_main_menu_real(self, mock_print, mock_input):
        """Test main menu selection for real datasets."""
        wizard = DatasetWizard()
        result = wizard.show_main_menu()
        assert result == "real"

    @patch('builtins.input', return_value='0')
    @patch('builtins.print')
    def test_show_main_menu_exit(self, mock_print, mock_input):
        """Test main menu selection for exit."""
        wizard = DatasetWizard()
        result = wizard.show_main_menu()
        assert result == "exit"

    @patch('builtins.input', return_value='invalid')
    @patch('builtins.print')
    def test_show_main_menu_invalid_choice(self, mock_print, mock_input):
        """Test main menu with invalid choice defaults to exit."""
        wizard = DatasetWizard()
        result = wizard.show_main_menu()
        assert result == "exit"

    @patch('src.presentation.cli.dataset_wizard._ask', return_value='random')
    @patch('src.presentation.cli.dataset_wizard._ask_int', return_value=100)
    @patch('src.presentation.cli.dataset_wizard._ask_float', return_value=0.5)
    @patch('src.presentation.cli.dataset_wizard._ask_optional_int', return_value=None)
    @patch('builtins.print')
    def test_collect_synthetic_params(self, mock_print, mock_opt_int, mock_float, mock_int, mock_ask):
        """Test collecting synthetic dataset parameters."""
        wizard = DatasetWizard()
        # Mock the input sequence for synthetic parameters
        mock_ask.side_effect = ['random', 'A', 'C', 'G', 'T']  # method, alphabet
        mock_int.side_effect = [100, 500, 50, 200]  # n, length, min_length, max_length
        mock_float.side_effect = [0.25, 0.25, 0.25, 0.25, 0.1]  # frequencies + noise
        mock_opt_int.side_effect = [None, None]  # Optional parameters
        
        result = wizard.collect_synthetic_params()
        
        assert isinstance(result, dict)
        assert 'method' in result
        assert 'n' in result

    @patch('src.presentation.cli.dataset_wizard._ask', return_value='human')
    @patch('src.presentation.cli.dataset_wizard._ask_int', return_value=100)
    @patch('builtins.print')
    def test_collect_real_params(self, mock_print, mock_int, mock_ask):
        """Test collecting real dataset parameters."""
        wizard = DatasetWizard()
        # Mock the input sequence for real parameters
        mock_ask.side_effect = ['human', 'protein', 'mitochondrion']  # query, database, organism
        mock_int.side_effect = [100, 50, 200]  # max_sequences, min_length, max_length
        
        result = wizard.collect_real_params()
        
        assert isinstance(result, dict)
        assert 'query' in result
        assert 'database' in result
        assert 'organism' in result
        assert 'max_sequences' in result

    def test_generate_default_filename_synthetic(self):
        """Test generating default filename for synthetic datasets."""
        wizard = DatasetWizard()
        params = {
            'method': 'random',
            'n': 100,
            'length': 500
        }
        result = wizard.generate_default_filename('synthetic', params)
        assert result == 'synthetic_random_n100_L500.fasta'

    def test_generate_default_filename_real(self):
        """Test generating default filename for real datasets."""
        wizard = DatasetWizard()
        params = {
            'query': 'human mitochondria'
        }
        result = wizard.generate_default_filename('real', params)
        assert result == 'real_human_mitochondria.fasta'

    def test_generate_default_filename_real_safe_chars(self):
        """Test generating default filename with special characters."""
        wizard = DatasetWizard()
        params = {
            'query': 'human/mitochondria[homo sapiens]'
        }
        result = wizard.generate_default_filename('real', params)
        assert 'real_' in result
        assert '.fasta' in result
        # Special characters should be replaced with underscores
        assert '[' not in result
        assert ']' not in result
        assert '/' not in result

    def test_generate_default_filename_unknown_kind(self):
        """Test generating default filename for unknown dataset kind."""
        wizard = DatasetWizard()
        params = {}
        result = wizard.generate_default_filename('unknown', params)
        assert result == 'dataset.fasta'

    @patch('src.presentation.cli.dataset_wizard._ask', return_value='custom_name')
    @patch('builtins.print')
    def test_get_output_filename_with_extension(self, mock_print, mock_ask):
        """Test getting output filename that already has .fasta extension."""
        wizard = DatasetWizard()
        mock_ask.return_value = 'custom_name.fasta'
        
        result = wizard.get_output_filename('default.fasta')
        assert result == 'custom_name.fasta'

    @patch('src.presentation.cli.dataset_wizard._ask', return_value='custom_name')
    @patch('builtins.print')
    def test_get_output_filename_without_extension(self, mock_print, mock_ask):
        """Test getting output filename that needs .fasta extension added."""
        wizard = DatasetWizard()
        mock_ask.return_value = 'custom_name'
        
        result = wizard.get_output_filename('default.fasta')
        assert result == 'custom_name.fasta'

    @patch('src.presentation.cli.dataset_wizard._ask', return_value='')
    @patch('builtins.print')
    def test_get_output_filename_use_default(self, mock_print, mock_ask):
        """Test getting output filename using default."""
        wizard = DatasetWizard()
        mock_ask.return_value = 'default.fasta'
        
        result = wizard.get_output_filename('default.fasta')
        assert result == 'default.fasta'


class TestIntegration:
    """Test integration scenarios for the dataset wizard."""

    @patch('builtins.input')
    @patch('builtins.print')
    def test_full_synthetic_workflow(self, mock_print, mock_input):
        """Test complete synthetic dataset generation workflow."""
        # Simulate user choosing synthetic, then providing all parameters
        mock_input.side_effect = [
            '1',  # Choose synthetic
            'random',  # Method
            '100',  # Number of sequences
            '500',  # Sequence length
            'ACGT',  # Alphabet
            '0.25', '0.25', '0.25', '0.25',  # Base frequencies
            '0.1',  # Noise level
            '',  # Seed (default)
            '',  # Min length (default)
            'my_dataset.fasta'  # Output filename
        ]
        
        wizard = DatasetWizard()
        
        # Test menu selection
        menu_choice = wizard.show_main_menu()
        assert menu_choice == 'synthetic'
        
        # Test parameter collection
        params = wizard.collect_synthetic_params()
        assert params['method'] == 'random'
        assert params['n'] == 100
        
        # Test filename generation
        default_name = wizard.generate_default_filename('synthetic', params)
        assert 'synthetic_random' in default_name
        
        # Test filename selection
        output_name = wizard.get_output_filename(default_name)
        assert output_name.endswith('.fasta')

    @patch('builtins.input')
    @patch('builtins.print')
    def test_full_real_workflow(self, mock_print, mock_input):
        """Test complete real dataset retrieval workflow."""
        # Simulate user choosing real, then providing all parameters
        mock_input.side_effect = [
            '2',  # Choose real
            'human mitochondria',  # Query
            'nucleotide',  # Database
            'homo sapiens',  # Organism
            '50',  # Max sequences
            '100',  # Min length
            '2000',  # Max length
            'real_dataset.fasta'  # Output filename
        ]
        
        wizard = DatasetWizard()
        
        # Test menu selection
        menu_choice = wizard.show_main_menu()
        assert menu_choice == 'real'
        
        # Test parameter collection
        params = wizard.collect_real_params()
        assert params['query'] == 'human mitochondria'
        assert params['database'] == 'nucleotide'
        
        # Test filename generation
        default_name = wizard.generate_default_filename('real', params)
        assert 'real_' in default_name
        
        # Test filename selection
