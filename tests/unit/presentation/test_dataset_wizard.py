"""
Unit tests for DatasetWizard (CLI).

Tests the command-line interface for dataset generation
including user input collection and parameter validation.
"""

from unittest.mock import patch, MagicMock
from io import StringIO
import sys

import pytest

from src.presentation.cli.dataset_wizard import (
    DatasetWizard,
    _ask,
    _ask_int,
    _ask_float,
)


class TestDatasetWizardHelperFunctions:
    """Tests for helper functions used by DatasetWizard."""

    def test_ask_with_default(self):
        """Test _ask function with default value."""
        with patch("builtins.input", return_value=""):
            result = _ask("Enter value", default="default_value")
            assert result == "default_value"

    def test_ask_with_user_input(self):
        """Test _ask function with user input."""
        with patch("builtins.input", return_value="user_input"):
            result = _ask("Enter value", default="default_value")
            assert result == "user_input"

    def test_ask_no_default(self):
        """Test _ask function without default."""
        with patch("builtins.input", return_value="user_input"):
            result = _ask("Enter value")
            assert result == "user_input"

    def test_ask_int_valid(self):
        """Test _ask_int with valid integer."""
        with patch("builtins.input", return_value="42"):
            result = _ask_int("Enter number", default=10)
            assert result == 42

    def test_ask_int_default(self):
        """Test _ask_int using default value."""
        with patch("builtins.input", return_value=""):
            result = _ask_int("Enter number", default=10)
            assert result == 10

    def test_ask_int_invalid_then_valid(self):
        """Test _ask_int with invalid input followed by valid."""
        with patch("builtins.input", side_effect=["invalid", "42"]):
            with patch("builtins.print"):  # Suppress warning output
                result = _ask_int("Enter number", default=10)
                assert result == 42

    def test_ask_int_below_minimum(self):
        """Test _ask_int with value below minimum."""
        with patch("builtins.input", side_effect=["0", "5"]):
            with patch("builtins.print"):  # Suppress warning output
                result = _ask_int("Enter number", default=10, min_value=1)
                assert result == 5

    def test_ask_float_valid(self):
        """Test _ask_float with valid float."""
        with patch("builtins.input", return_value="0.5"):
            result = _ask_float("Enter rate", default=0.2)
            assert result == 0.5

    def test_ask_float_default(self):
        """Test _ask_float using default value."""
        with patch("builtins.input", return_value=""):
            result = _ask_float("Enter rate", default=0.2)
            assert result == 0.2

    def test_ask_float_out_of_bounds(self):
        """Test _ask_float with out of bounds values."""
        with patch("builtins.input", side_effect=["-0.1", "1.5", "0.3"]):
            with patch("builtins.print"):  # Suppress warning output
                result = _ask_float(
                    "Enter rate", default=0.2, min_value=0.0, max_value=1.0
                )
                assert result == 0.3

    def test_ask_float_invalid_then_valid(self):
        """Test _ask_float with invalid input followed by valid."""
        with patch("builtins.input", side_effect=["not_a_number", "0.7"]):
            with patch("builtins.print"):  # Suppress warning output
                result = _ask_float("Enter rate", default=0.2)
                assert result == 0.7


class TestDatasetWizardMainMenu:
    """Tests for main menu functionality."""

    def test_show_main_menu_synthetic(self):
        """Test selecting synthetic dataset generation."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="1"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "synthetic"

    def test_show_main_menu_real(self):
        """Test selecting real dataset download."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="2"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "real"

    def test_show_main_menu_exit(self):
        """Test selecting exit option."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="0"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "exit"

    def test_show_main_menu_default(self):
        """Test default selection (synthetic)."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value=""):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "synthetic"

    def test_show_main_menu_invalid_input(self):
        """Test invalid input defaults to exit."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="invalid"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "exit"


class TestDatasetWizardSyntheticParams:
    """Tests for synthetic parameter collection."""

    def test_collect_synthetic_params_random(self):
        """Test collecting parameters for random synthetic dataset."""
        wizard = DatasetWizard()

        inputs = [
            "1",  # method (random)
            "10",  # n
            "8",  # length
            "ACGT",  # alphabet
            "42",  # seed
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

        expected = {
            "method": "random",
            "n": 10,
            "length": 8,
            "alphabet": "ACGT",
            "seed": 42,
        }
        assert params == expected

    def test_collect_synthetic_params_noise(self):
        """Test collecting parameters for noise-based synthetic dataset."""
        wizard = DatasetWizard()

        inputs = [
            "2",  # method (noise)
            "5",  # n
            "6",  # length
            "ACGT",  # alphabet
            "",  # seed (empty, so None)
            "0.3",  # noise_rate
            "ACGTAC",  # center_sequence
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

        expected = {
            "method": "noise",
            "n": 5,
            "length": 6,
            "alphabet": "ACGT",
            "noise_rate": 0.3,
            "center_sequence": "ACGTAC",
        }
        assert params == expected

    def test_collect_synthetic_params_clustered(self):
        """Test collecting parameters for clustered synthetic dataset."""
        wizard = DatasetWizard()

        inputs = [
            "3",  # method (clustered)
            "12",  # n
            "10",  # length
            "ACGT",  # alphabet
            "123",  # seed
            "3",  # num_clusters
            "0.1",  # cluster_noise
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

        expected = {
            "method": "clustered",
            "n": 12,
            "length": 10,
            "alphabet": "ACGT",
            "seed": 123,
            "num_clusters": 3,
            "cluster_noise": 0.1,
        }
        assert params == expected

    def test_collect_synthetic_params_mutations(self):
        """Test collecting parameters for mutations-based synthetic dataset."""
        wizard = DatasetWizard()

        inputs = [
            "4",  # method (mutations)
            "8",  # n
            "12",  # length
            "ACGT",  # alphabet
            "",  # seed (empty)
            "0.15",  # mutation_rate
            "ACGTACGTACGT",  # base_sequence
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

        expected = {
            "method": "mutations",
            "n": 8,
            "length": 12,
            "alphabet": "ACGT",
            "mutation_rate": 0.15,
            "base_sequence": "ACGTACGTACGT",
        }
        assert params == expected

    def test_collect_synthetic_params_defaults(self):
        """Test collecting parameters using all defaults."""
        wizard = DatasetWizard()

        # Use empty strings to trigger defaults
        inputs = ["", "", "", "", ""]  # All defaults

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                with patch(
                    "src.presentation.cli.dataset_wizard.SYNTHETIC_DEFAULTS",
                    {"n": 20, "length": 50, "alphabet": "ACGT"},
                ):
                    params = wizard.collect_synthetic_params()

        assert params["method"] == "random"  # Default method
        assert params["n"] == 20
        assert params["length"] == 50
        assert params["alphabet"] == "ACGT"
        assert "seed" not in params  # Empty seed should not be included

    def test_collect_synthetic_params_invalid_method(self):
        """Test handling invalid method selection."""
        wizard = DatasetWizard()

        inputs = [
            "invalid",  # invalid method choice
            "5",  # n
            "8",  # length
            "ACGT",  # alphabet
            "",  # seed
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

        # Should default to 'random' for invalid method
        assert params["method"] == "random"


class TestDatasetWizardRealParams:
    """Tests for real dataset parameter collection."""

    def test_collect_real_params_ncbi(self):
        """Test collecting parameters for NCBI download."""
        wizard = DatasetWizard()

        inputs = [
            "1",  # source (ncbi)
            "COX1 AND mitochondrion",  # query
            "25",  # max_sequences
            "100",  # min_length
            "200",  # max_length
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_real_params()

        expected = {
            "source": "ncbi",
            "query": "COX1 AND mitochondrion",
            "max_sequences": 25,
            "min_length": 100,
            "max_length": 200,
        }
        assert params == expected

    def test_collect_real_params_file(self):
        """Test collecting parameters for file import."""
        wizard = DatasetWizard()

        inputs = ["2", "data/sequences.fasta"]  # source (file)  # file_path

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_real_params()

        expected = {"source": "file", "file_path": "data/sequences.fasta"}
        assert params == expected

    def test_collect_real_params_defaults(self):
        """Test collecting real parameters with defaults."""
        wizard = DatasetWizard()

        inputs = ["", "", "", "", ""]  # All defaults

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                with patch(
                    "src.presentation.cli.dataset_wizard.ENTREZ_DEFAULTS",
                    {"query": "test query", "max_sequences": 30},
                ):
                    params = wizard.collect_real_params()

        assert params["source"] == "ncbi"  # Default source
        assert "query" in params
        assert "max_sequences" in params

    def test_collect_real_params_invalid_source(self):
        """Test handling invalid source selection."""
        wizard = DatasetWizard()

        inputs = [
            "invalid",  # invalid source
            "test query",  # query (will be used since defaults to ncbi)
            "10",  # max_sequences
            "50",  # min_length
            "100",  # max_length
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_real_params()

        # Should default to 'ncbi' for invalid source
        assert params["source"] == "ncbi"


class TestDatasetWizardUtilityMethods:
    """Tests for utility methods."""

    def test_generate_default_filename_synthetic_random(self):
        """Test generating filename for synthetic random dataset."""
        wizard = DatasetWizard()

        params = {"method": "random", "n": 10, "length": 8}

        filename = wizard.generate_default_filename("synthetic", params)

        assert filename == "synthetic_random_n10_L8.fasta"

    def test_generate_default_filename_synthetic_noise(self):
        """Test generating filename for synthetic noise dataset."""
        wizard = DatasetWizard()

        params = {"method": "noise", "n": 5, "length": 12}

        filename = wizard.generate_default_filename("synthetic", params)

        assert filename == "synthetic_noise_n5_L12.fasta"

    def test_generate_default_filename_real_ncbi(self):
        """Test generating filename for real NCBI dataset."""
        wizard = DatasetWizard()

        params = {"source": "ncbi", "query": "COX1 AND mitochondrion"}

        filename = wizard.generate_default_filename("real", params)

        assert filename == "real_COX1_AND_mitochondrion.fasta"

    def test_generate_default_filename_real_file(self):
        """Test generating filename for real file import."""
        wizard = DatasetWizard()

        params = {"source": "file"}

        filename = wizard.generate_default_filename("real", params)

        assert filename == "real_import.fasta"

    def test_generate_default_filename_special_characters(self):
        """Test filename generation with special characters in query."""
        wizard = DatasetWizard()

        params = {"source": "ncbi", "query": 'gene[ALL] AND species:"Homo sapiens"'}

        filename = wizard.generate_default_filename("real", params)

        # Special characters should be replaced with underscores
        assert "gene_ALL__AND_species" in filename
        assert filename.endswith(".fasta")

    def test_generate_default_filename_long_query(self):
        """Test filename generation with very long query."""
        wizard = DatasetWizard()

        params = {"source": "ncbi", "query": "a" * 50}  # Very long query

        filename = wizard.generate_default_filename("real", params)

        # Should be truncated to 30 characters plus extension
        assert len(filename) <= 36  # 30 chars + ".fasta" = 36

    def test_generate_default_filename_fallback(self):
        """Test filename generation fallback."""
        wizard = DatasetWizard()

        # Invalid kind should return generic filename
        filename = wizard.generate_default_filename("invalid", {})

        assert filename == "dataset.fasta"

    def test_get_output_filename_custom(self):
        """Test getting custom output filename."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="my_custom_dataset.fasta"):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default.fasta")

        assert filename == "my_custom_dataset.fasta"

    def test_get_output_filename_default(self):
        """Test using default output filename."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value=""):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default_name.fasta")

        assert filename == "default_name.fasta"

    def test_get_output_filename_auto_extension(self):
        """Test automatic addition of .fasta extension."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="my_dataset"):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default.fasta")

        assert filename == "my_dataset.fasta"

    def test_get_output_filename_preserve_extension(self):
        """Test preserving existing .fasta extension."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="my_dataset.FASTA"):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default.fasta")

        # Should preserve the extension (case insensitive check)
        assert filename == "my_dataset.FASTA"


class TestDatasetWizardIntegration:
    """Integration tests for DatasetWizard workflow."""

    def test_complete_synthetic_workflow(self):
        """Test complete synthetic dataset generation workflow."""
        wizard = DatasetWizard()

        # Simulate complete user interaction
        with patch(
            "builtins.input",
            side_effect=[
                "1",  # synthetic choice in main menu
                "1",  # random method
                "10",  # n
                "8",  # length
                "ACGT",  # alphabet
                "42",  # seed
                "",  # use default filename
            ],
        ):
            with patch("builtins.print"):
                # Main menu
                choice = wizard.show_main_menu()
                assert choice == "synthetic"

                # Collect parameters
                params = wizard.collect_synthetic_params()
                assert params["method"] == "random"
                assert params["n"] == 10

                # Generate filename
                default_name = wizard.generate_default_filename("synthetic", params)
                final_name = wizard.get_output_filename(default_name)
                assert final_name == "synthetic_random_n10_L8.fasta"

    def test_complete_real_workflow(self):
        """Test complete real dataset download workflow."""
        wizard = DatasetWizard()

        with patch(
            "builtins.input",
            side_effect=[
                "2",  # real choice
                "1",  # ncbi source
                "COX1",  # query
                "20",  # max_sequences
                "50",  # min_length
                "100",  # max_length
                "my_cox1_data",  # custom filename
            ],
        ):
            with patch("builtins.print"):
                # Main menu
                choice = wizard.show_main_menu()
                assert choice == "real"

                # Collect parameters
                params = wizard.collect_real_params()
                assert params["source"] == "ncbi"
                assert params["query"] == "COX1"

                # Generate and get filename
                default_name = wizard.generate_default_filename("real", params)
                final_name = wizard.get_output_filename(default_name)
                assert final_name == "my_cox1_data.fasta"

    def test_user_exit_workflow(self):
        """Test user exit workflow."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="0"):
            with patch("builtins.print"):
                choice = wizard.show_main_menu()
                assert choice == "exit"


class TestDatasetWizardEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_empty_inputs_handling(self):
        """Test handling of empty inputs throughout workflow."""
        wizard = DatasetWizard()

        # Test with mostly empty inputs (should use defaults)
        inputs = ["", "", "", "", "", ""]  # All empty

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                # Should not crash and use reasonable defaults
                choice = wizard.show_main_menu()
                assert choice == "synthetic"  # Default choice

    def test_case_sensitivity(self):
        """Test case sensitivity in inputs."""
        wizard = DatasetWizard()

        inputs = [
            "1",  # method
            "5",  # n
            "6",  # length
            "acgt",  # lowercase alphabet
            "",  # seed
        ]

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

        # Alphabet should be converted to uppercase
        assert params["alphabet"] == "ACGT"

    def test_whitespace_handling(self):
        """Test handling of whitespace in inputs."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="  spaced query  "):
            result = _ask("Enter query")
            assert result == "spaced query"  # Should be stripped
