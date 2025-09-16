"""
Comprehensive tests for DatasetWizard CLI module.

This module provides comprehensive test coverage for the dataset wizard CLI,
testing all user interaction functions, validation logic, parameter collection
methods, and edge cases.
"""

from unittest.mock import MagicMock, patch

import pytest

from src.presentation.cli.dataset_wizard import (
    DatasetWizard,
    _ask,
    _ask_float,
    _ask_int,
    _ask_optional_int,
)


class TestHelperFunctions:
    """Test helper functions used by DatasetWizard."""

    def test_ask_with_default(self):
        """Test _ask function with default value."""
        with patch("builtins.input", return_value=""):
            result = _ask("Enter value", "default_value")
            assert result == "default_value"

    def test_ask_with_user_input(self):
        """Test _ask function with user input."""
        with patch("builtins.input", return_value="user_input"):
            result = _ask("Enter value", "default_value")
            assert result == "user_input"

    def test_ask_no_default_empty_input(self):
        """Test _ask function with no default and empty input."""
        with patch("builtins.input", return_value=""):
            result = _ask("Enter value")
            assert result == ""

    def test_ask_strips_whitespace(self):
        """Test _ask function strips whitespace."""
        with patch("builtins.input", return_value="  value  "):
            result = _ask("Enter value")
            assert result == "value"

    def test_ask_int_valid_input(self):
        """Test _ask_int with valid integer input."""
        with patch("builtins.input", return_value="42"):
            result = _ask_int("Enter number", 10)
            assert result == 42

    def test_ask_int_default_value(self):
        """Test _ask_int with empty input returns default."""
        with patch("builtins.input", return_value=""):
            result = _ask_int("Enter number", 15)
            assert result == 15

    def test_ask_int_invalid_then_valid(self):
        """Test _ask_int with invalid input followed by valid input."""
        with patch("builtins.input", side_effect=["invalid", "25"]):
            with patch("builtins.print") as mock_print:
                result = _ask_int("Enter number", 10)
                assert result == 25
                mock_print.assert_called_with("⚠️  Please enter a valid integer number.")

    def test_ask_int_below_minimum(self):
        """Test _ask_int with value below minimum."""
        with patch("builtins.input", side_effect=["0", "5"]):
            with patch("builtins.print") as mock_print:
                result = _ask_int("Enter number", 10, min_value=1)
                assert result == 5
                mock_print.assert_called_with("⚠️  Value must be >= 1.")

    def test_ask_int_no_minimum_constraint(self):
        """Test _ask_int with no minimum constraint."""
        with patch("builtins.input", return_value="-5"):
            result = _ask_int("Enter number", 10, min_value=None)
            assert result == -5

    def test_ask_float_valid_input(self):
        """Test _ask_float with valid float input."""
        with patch("builtins.input", return_value="0.5"):
            result = _ask_float("Enter float", 0.1)
            assert result == 0.5

    def test_ask_float_default_value(self):
        """Test _ask_float with empty input returns default."""
        with patch("builtins.input", return_value=""):
            result = _ask_float("Enter float", 0.8)
            assert result == 0.8

    def test_ask_float_invalid_then_valid(self):
        """Test _ask_float with invalid input followed by valid input."""
        with patch("builtins.input", side_effect=["not_a_float", "0.3"]):
            with patch("builtins.print") as mock_print:
                result = _ask_float("Enter float", 0.1)
                assert result == 0.3
                mock_print.assert_called_with("⚠️  Please enter a valid numeric value.")

    def test_ask_float_below_minimum(self):
        """Test _ask_float with value below minimum."""
        with patch("builtins.input", side_effect=["-1.0", "0.2"]):
            with patch("builtins.print") as mock_print:
                result = _ask_float("Enter float", 0.1, min_value=0.0)
                assert result == 0.2
                mock_print.assert_called_with("⚠️  Value must be >= 0.0.")

    def test_ask_float_above_maximum(self):
        """Test _ask_float with value above maximum."""
        with patch("builtins.input", side_effect=["2.0", "0.8"]):
            with patch("builtins.print") as mock_print:
                result = _ask_float("Enter float", 0.5, max_value=1.0)
                assert result == 0.8
                mock_print.assert_called_with("⚠️  Value must be <= 1.0.")

    def test_ask_float_no_constraints(self):
        """Test _ask_float with no min/max constraints."""
        with patch("builtins.input", return_value="5.5"):
            result = _ask_float("Enter float", 1.0, min_value=None, max_value=None)
            assert result == 5.5

    def test_ask_optional_int_valid_input(self):
        """Test _ask_optional_int with valid integer input."""
        with patch("builtins.input", return_value="100"):
            result = _ask_optional_int("Enter optional number", 50)
            assert result == 100

    def test_ask_optional_int_empty_returns_default(self):
        """Test _ask_optional_int with empty input returns default."""
        with patch("builtins.input", return_value=""):
            result = _ask_optional_int("Enter optional number", 25)
            assert result == 25

    def test_ask_optional_int_empty_returns_none(self):
        """Test _ask_optional_int with empty input and None default."""
        with patch("builtins.input", return_value=""):
            result = _ask_optional_int("Enter optional number", None)
            assert result is None

    def test_ask_optional_int_invalid_then_valid(self):
        """Test _ask_optional_int with invalid input followed by valid input."""
        with patch("builtins.input", side_effect=["invalid", "75"]):
            with patch("builtins.print") as mock_print:
                result = _ask_optional_int("Enter optional number", 10)
                assert result == 75
                mock_print.assert_called_with(
                    "⚠️  Please enter a valid integer number or leave empty for None."
                )


class TestDatasetWizard:
    """Test DatasetWizard class methods."""

    def test_show_main_menu_synthetic(self):
        """Test main menu selection for synthetic dataset."""
        wizard = DatasetWizard()
        with patch("builtins.input", return_value="1"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "synthetic"

    def test_show_main_menu_real(self):
        """Test main menu selection for real dataset."""
        wizard = DatasetWizard()
        with patch("builtins.input", return_value="2"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "real"

    def test_show_main_menu_exit(self):
        """Test main menu selection for exit."""
        wizard = DatasetWizard()
        with patch("builtins.input", return_value="0"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "exit"

    def test_show_main_menu_invalid_choice(self):
        """Test main menu with invalid choice defaults to exit."""
        wizard = DatasetWizard()
        with patch("builtins.input", return_value="invalid"):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "exit"

    def test_show_main_menu_default_choice(self):
        """Test main menu with empty input uses default (synthetic)."""
        wizard = DatasetWizard()
        with patch("builtins.input", return_value=""):
            with patch("builtins.print"):
                result = wizard.show_main_menu()
                assert result == "synthetic"

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_collect_synthetic_params_random_method(self):
        """Test collecting synthetic parameters for random method."""
        wizard = DatasetWizard()
        inputs = ["1", "200", "75", "ATCG", "42"]  # method, n, length, alphabet, seed

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

                assert params["method"] == "random"
                assert params["n"] == 200
                assert params["length"] == 75
                assert params["alphabet"] == "ATCG"
                assert params["seed"] == 42

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_collect_synthetic_params_noise_method(self):
        """Test collecting synthetic parameters for noise method."""
        wizard = DatasetWizard()
        inputs = [
            "2",
            "150",
            "60",
            "ATCG",
            "",
            "0.2",
            "ATCGATCG",
        ]  # method=noise, n, length, alphabet, seed, noise_rate, center

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

                assert params["method"] == "noise"
                assert params["n"] == 150
                assert params["length"] == 60
                assert params["alphabet"] == "ATCG"
                assert "seed" not in params  # empty seed should not be included
                assert params["noise_rate"] == 0.2
                assert params["center_sequence"] == "ATCGATCG"

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_collect_synthetic_params_noise_method_no_center(self):
        """Test collecting synthetic parameters for noise method without center sequence."""
        wizard = DatasetWizard()
        inputs = [
            "2",
            "150",
            "60",
            "ATCG",
            "",
            "0.2",
            "",
        ]  # method=noise, n, length, alphabet, seed, noise_rate, empty center

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

                assert params["method"] == "noise"
                assert params["noise_rate"] == 0.2
                assert (
                    "center_sequence" not in params
                )  # empty center should not be included

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_collect_synthetic_params_clustered_method(self):
        """Test collecting synthetic parameters for clustered method."""
        wizard = DatasetWizard()
        inputs = [
            "3",
            "120",
            "80",
            "ATCG",
            "123",
            "5",
            "0.15",
        ]  # method=clustered, n, length, alphabet, seed, num_clusters, cluster_noise

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

                assert params["method"] == "clustered"
                assert params["n"] == 120
                assert params["length"] == 80
                assert params["alphabet"] == "ATCG"
                assert params["seed"] == 123
                assert params["num_clusters"] == 5
                assert params["cluster_noise"] == 0.15

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_collect_synthetic_params_mutations_method(self):
        """Test collecting synthetic parameters for mutations method."""
        wizard = DatasetWizard()
        inputs = [
            "4",
            "180",
            "90",
            "ATCG",
            "456",
            "0.08",
            "ATCGATCGATCG",
        ]  # method=mutations, n, length, alphabet, seed, mutation_rate, base_sequence

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

                assert params["method"] == "mutations"
                assert params["n"] == 180
                assert params["length"] == 90
                assert params["alphabet"] == "ATCG"
                assert params["seed"] == 456
                assert params["mutation_rate"] == 0.08
                assert params["base_sequence"] == "ATCGATCGATCG"

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_collect_synthetic_params_mutations_method_no_base(self):
        """Test collecting synthetic parameters for mutations method without base sequence."""
        wizard = DatasetWizard()
        inputs = [
            "4",
            "180",
            "90",
            "ATCG",
            "456",
            "0.08",
            "",
        ]  # method=mutations, n, length, alphabet, seed, mutation_rate, empty base

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

                assert params["method"] == "mutations"
                assert params["mutation_rate"] == 0.08
                assert (
                    "base_sequence" not in params
                )  # empty base should not be included

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_collect_synthetic_params_invalid_method_choice(self):
        """Test collecting synthetic parameters with invalid method choice defaults to random."""
        wizard = DatasetWizard()
        inputs = [
            "invalid",
            "100",
            "50",
            "ATCG",
            "",
        ]  # invalid method choice, n, length, alphabet, empty seed

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()

                assert params["method"] == "random"  # should default to random
                assert params["n"] == 100
                assert params["length"] == 50
                assert params["alphabet"] == "ATCG"
                assert "seed" not in params

    @patch(
        "src.infrastructure.external.dataset_entrez.ENTREZ_DEFAULTS",
        {
            "db": "nucleotide",
            "query": "*",
            "max_sequences": 1000,
            "min_length": None,
            "max_length": None,
        },
    )
    def test_collect_real_params_all_values(self):
        """Test collecting real dataset parameters with all values specified."""
        wizard = DatasetWizard()
        inputs = [
            "protein",
            "gene[All Fields]",
            "500",
            "50",
            "200",
        ]  # db, query, max_seq, min_len, max_len

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_real_params()

                assert params["db"] == "protein"
                assert params["query"] == "gene[All Fields]"
                assert params["max_sequences"] == 500
                assert params["min_length"] == 50
                assert params["max_length"] == 200

    @patch(
        "src.infrastructure.external.dataset_entrez.ENTREZ_DEFAULTS",
        {
            "db": "nucleotide",
            "query": "*",
            "max_sequences": 1000,
            "min_length": None,
            "max_length": None,
        },
    )
    def test_collect_real_params_with_defaults(self):
        """Test collecting real dataset parameters using defaults."""
        wizard = DatasetWizard()
        inputs = ["", "", "", "", ""]  # all empty, should use defaults

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print"):
                params = wizard.collect_real_params()

                assert params["db"] == "nucleotide"  # default value
                assert params["query"] == "*"  # fallback for empty query
                assert params["max_sequences"] == 1000  # default value
                assert params["min_length"] is None  # default value
                assert params["max_length"] is None  # default value

    @patch(
        "src.infrastructure.external.dataset_entrez.ENTREZ_DEFAULTS",
        {
            "db": "nucleotide",
            "query": "*",
            "max_sequences": 1000,
            "min_length": None,
            "max_length": None,
        },
    )
    def test_collect_real_params_min_max_validation(self):
        """Test collecting real dataset parameters with min/max length validation."""
        wizard = DatasetWizard()
        inputs = [
            "nucleotide",
            "bacteria[organism]",
            "300",
            "100",
            "50",
            "150",
        ]  # db, query, max_seq, min_len, invalid_max_len, corrected_max_len

        with patch("builtins.input", side_effect=inputs):
            with patch("builtins.print") as mock_print:
                params = wizard.collect_real_params()

                assert params["db"] == "nucleotide"
                assert params["query"] == "bacteria[organism]"
                assert params["max_sequences"] == 300
                assert params["min_length"] == 100
                assert params["max_length"] == 150  # corrected value
                # Should print validation error
                mock_print.assert_any_call(
                    "⚠️  'Maximum length' must be >= 'Minimum length'."
                )

    def test_generate_default_filename_synthetic(self):
        """Test generating default filename for synthetic dataset."""
        wizard = DatasetWizard()
        params = {"method": "random", "n": 100, "length": 50}

        filename = wizard.generate_default_filename("synthetic", params)
        assert filename == "synthetic_random_n100_L50.fasta"

    def test_generate_default_filename_real(self):
        """Test generating default filename for real dataset."""
        wizard = DatasetWizard()
        params = {"query": "bacteria[organism] AND complete genome"}

        filename = wizard.generate_default_filename("real", params)
        assert filename.startswith("real_bacteria_organism_AND_")
        assert filename.endswith(".fasta")

    def test_generate_default_filename_real_long_query(self):
        """Test generating default filename for real dataset with long query."""
        wizard = DatasetWizard()
        long_query = "a" * 100  # Very long query
        params = {"query": long_query}

        filename = wizard.generate_default_filename("real", params)
        # Should be truncated to 30 chars plus prefix and suffix
        assert len(filename) <= len("real_" + "a" * 30 + ".fasta")
        assert filename.endswith(".fasta")

    def test_generate_default_filename_real_special_chars(self):
        """Test generating default filename for real dataset with special characters in query."""
        wizard = DatasetWizard()
        params = {"query": "gene[All] AND species:Homo sapiens & complete!"}

        filename = wizard.generate_default_filename("real", params)
        # Special characters should be replaced with underscores and truncated
        assert "gene_All_AND_species_Homo_sapi" in filename
        assert filename.endswith(".fasta")

    def test_generate_default_filename_real_empty_query(self):
        """Test generating default filename for real dataset with empty query."""
        wizard = DatasetWizard()
        params = {"query": ""}

        filename = wizard.generate_default_filename("real", params)
        assert filename == "real_query.fasta"  # fallback to 'query'

    def test_generate_default_filename_unknown_kind(self):
        """Test generating default filename for unknown dataset kind."""
        wizard = DatasetWizard()
        params = {}

        filename = wizard.generate_default_filename("unknown", params)
        assert filename == "dataset.fasta"

    def test_get_output_filename_with_extension(self):
        """Test getting output filename when user provides name with .fasta extension."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="my_dataset.fasta"):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default.fasta")
                assert filename == "my_dataset.fasta"

    def test_get_output_filename_without_extension(self):
        """Test getting output filename when user provides name without extension."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="my_dataset"):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default.fasta")
                assert filename == "my_dataset.fasta"

    def test_get_output_filename_empty_uses_default(self):
        """Test getting output filename when user provides empty input."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value=""):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("suggested.fasta")
                assert filename == "suggested.fasta"

    def test_get_output_filename_case_insensitive_extension(self):
        """Test getting output filename with case-insensitive extension checking."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="my_dataset.FASTA"):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default.fasta")
                assert filename == "my_dataset.FASTA"  # preserves original case

    def test_get_output_filename_adds_extension_to_uppercase(self):
        """Test getting output filename adds extension when uppercase extension missing."""
        wizard = DatasetWizard()

        with patch("builtins.input", return_value="my_dataset.txt"):
            with patch("builtins.print"):
                filename = wizard.get_output_filename("default.fasta")
                assert filename == "my_dataset.txt.fasta"


class TestIntegrationScenarios:
    """Test integration scenarios combining multiple methods."""

    @patch(
        "src.application.services.dataset_generator.SYNTHETIC_DEFAULTS",
        {
            "n": 100,
            "length": 50,
            "alphabet": "ATCG",
            "noise_rate": 0.1,
            "num_clusters": 3,
            "cluster_distance": 0.2,
            "mutation_rate": 0.05,
        },
    )
    def test_full_synthetic_workflow_with_defaults(self):
        """Test complete synthetic dataset workflow with default values."""
        wizard = DatasetWizard()

        # Menu choice: synthetic (1), then all defaults for random method
        menu_inputs = ["1"]  # choose synthetic
        param_inputs = [
            "",
            "",
            "",
            "",
            "",
        ]  # all defaults: method, n, length, alphabet, seed
        filename_inputs = [""]  # use default filename

        with patch("builtins.input", side_effect=menu_inputs):
            with patch("builtins.print"):
                choice = wizard.show_main_menu()
                assert choice == "synthetic"

        with patch("builtins.input", side_effect=param_inputs):
            with patch("builtins.print"):
                params = wizard.collect_synthetic_params()
                assert params["method"] == "random"  # default method choice (1)
                assert params["n"] == 100  # default value
                assert params["length"] == 50  # default value
                assert params["alphabet"] == "ATCG"  # default value

        default_filename = wizard.generate_default_filename("synthetic", params)
        assert default_filename == "synthetic_random_n100_L50.fasta"

        with patch("builtins.input", side_effect=filename_inputs):
            with patch("builtins.print"):
                final_filename = wizard.get_output_filename(default_filename)
                assert final_filename == default_filename

    @patch(
        "src.infrastructure.external.dataset_entrez.ENTREZ_DEFAULTS",
        {
            "db": "nucleotide",
            "query": "*",
            "max_sequences": 1000,
            "min_length": None,
            "max_length": None,
        },
    )
    def test_full_real_workflow_with_custom_values(self):
        """Test complete real dataset workflow with custom values."""
        wizard = DatasetWizard()

        # Menu choice: real (2)
        menu_inputs = ["2"]
        param_inputs = [
            "protein",
            "bacteria[organism]",
            "500",
            "100",
            "300",
        ]  # custom values
        filename_inputs = ["bacteria_proteins"]  # custom filename

        with patch("builtins.input", side_effect=menu_inputs):
            with patch("builtins.print"):
                choice = wizard.show_main_menu()
                assert choice == "real"

        with patch("builtins.input", side_effect=param_inputs):
            with patch("builtins.print"):
                params = wizard.collect_real_params()
                assert params["db"] == "protein"
                assert params["query"] == "bacteria[organism]"
                assert params["max_sequences"] == 500
                assert params["min_length"] == 100
                assert params["max_length"] == 300

        default_filename = wizard.generate_default_filename("real", params)
        assert "bacteria_organism" in default_filename

        with patch("builtins.input", side_effect=filename_inputs):
            with patch("builtins.print"):
                final_filename = wizard.get_output_filename(default_filename)
                assert final_filename == "bacteria_proteins.fasta"
