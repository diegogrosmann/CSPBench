"""
Unit tests for DatasetGenerationOrchestrator.

Tests the high-level orchestrator that coordinates dataset
generation workflows using the DatasetWizard and various
repository/generator services.
"""

from unittest.mock import Mock, patch, MagicMock
from pathlib import Path

import pytest

from src.domain.dataset import Dataset
from infrastructure.orchestration.dataset_generation_orchestrator import (
    DatasetGenerationOrchestrator,
)


class TestDatasetGenerationOrchestratorInitialization:
    """Tests for orchestrator initialization."""

    def test_init_default_path(self):
        """Test initialization with default base path."""
        orchestrator = DatasetGenerationOrchestrator()

        assert orchestrator.base_path == "datasets"
        assert orchestrator.wizard is not None
        assert orchestrator.repository is not None

    def test_init_custom_path(self):
        """Test initialization with custom base path."""
        custom_path = "/custom/datasets"
        orchestrator = DatasetGenerationOrchestrator(base_path=custom_path)

        assert orchestrator.base_path == custom_path

    def test_init_wizard_creation(self):
        """Test that wizard is properly initialized."""
        orchestrator = DatasetGenerationOrchestrator()

        # Should have a DatasetWizard instance
        from src.presentation.cli.dataset_wizard import DatasetWizard

        assert isinstance(orchestrator.wizard, DatasetWizard)

    def test_init_repository_assignment(self):
        """Test that repository is properly assigned."""
        orchestrator = DatasetGenerationOrchestrator()

        # Should have FileDatasetRepository class
        from src.infrastructure.persistence.dataset_repository import (
            FileDatasetRepository,
        )

        assert orchestrator.repository is FileDatasetRepository


class TestDatasetGenerationOrchestratorInteractiveGeneration:
    """Tests for interactive generation workflow."""

    @patch(
        "src.infrastructure.orchestrators.dataset_generation_orchestrator.DatasetWizard"
    )
    def test_run_interactive_generation_exit(self, mock_wizard_class):
        """Test interactive generation with user exit."""
        mock_wizard = Mock()
        mock_wizard.show_main_menu.return_value = "exit"
        mock_wizard_class.return_value = mock_wizard

        orchestrator = DatasetGenerationOrchestrator()

        result = orchestrator.run_interactive_generation()

        assert result is None
        mock_wizard.show_main_menu.assert_called_once()

    @patch(
        "src.infrastructure.orchestrators.dataset_generation_orchestrator.DatasetWizard"
    )
    def test_run_interactive_generation_keyboard_interrupt(self, mock_wizard_class):
        """Test handling of keyboard interrupt."""
        mock_wizard = Mock()
        mock_wizard.show_main_menu.side_effect = KeyboardInterrupt()
        mock_wizard_class.return_value = mock_wizard

        orchestrator = DatasetGenerationOrchestrator()

        with patch("builtins.print"):  # Suppress output
            result = orchestrator.run_interactive_generation()

        assert result is None

    @patch(
        "src.infrastructure.orchestrators.dataset_generation_orchestrator.DatasetWizard"
    )
    def test_run_interactive_generation_eof_error(self, mock_wizard_class):
        """Test handling of EOF error."""
        mock_wizard = Mock()
        mock_wizard.show_main_menu.side_effect = EOFError()
        mock_wizard_class.return_value = mock_wizard

        orchestrator = DatasetGenerationOrchestrator()

        with patch("builtins.print"):  # Suppress output
            result = orchestrator.run_interactive_generation()

        assert result is None

    @patch(
        "src.infrastructure.orchestrators.dataset_generation_orchestrator.DatasetWizard"
    )
    def test_run_interactive_generation_file_error(self, mock_wizard_class):
        """Test handling of file-related errors."""
        mock_wizard = Mock()
        mock_wizard.show_main_menu.side_effect = FileNotFoundError("File not found")
        mock_wizard_class.return_value = mock_wizard

        orchestrator = DatasetGenerationOrchestrator()

        with patch("builtins.print"):  # Suppress output
            result = orchestrator.run_interactive_generation()

        assert result is None


class TestDatasetGenerationOrchestratorSyntheticGeneration:
    """Tests for synthetic dataset generation workflow."""

    def test_handle_synthetic_generation_success(self):
        """Test successful synthetic dataset generation."""
        orchestrator = DatasetGenerationOrchestrator()

        # Mock wizard responses
        mock_params = {
            "method": "random",
            "n": 10,
            "length": 8,
            "alphabet": "ACGT",
            "seed": 42,
        }

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_synthetic_params.return_value = mock_params
        orchestrator.wizard.generate_default_filename.return_value = (
            "synthetic_random_n10_L8.fasta"
        )
        orchestrator.wizard.get_output_filename.return_value = "test_dataset.fasta"

        # Mock generator
        mock_dataset = Dataset(["ACGT"] * 10)
        mock_gen_params = {"method": "random", "n": 10}

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetGenerator"
        ) as mock_generator:
            mock_generator.generate_from_config.return_value = (
                mock_dataset,
                mock_gen_params,
            )

            # Mock repository
            with patch.object(orchestrator.repository, "save") as mock_save:
                mock_save.return_value = "/path/to/test_dataset.fasta"

                with patch("builtins.print"):  # Suppress output
                    result = orchestrator._handle_synthetic_generation()

        assert result == "/path/to/test_dataset.fasta"
        orchestrator.wizard.collect_synthetic_params.assert_called_once()
        mock_generator.generate_from_config.assert_called_once()
        mock_save.assert_called_once()

    def test_handle_synthetic_generation_with_config_creation(self):
        """Test synthetic generation with proper config object creation."""
        orchestrator = DatasetGenerationOrchestrator()

        mock_params = {
            "method": "noise",
            "n": 5,
            "length": 6,
            "alphabet": "ACGT",
            "noise_rate": 0.2,
            "seed": 123,
        }

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_synthetic_params.return_value = mock_params
        orchestrator.wizard.generate_default_filename.return_value = (
            "synthetic_noise_n5_L6.fasta"
        )
        orchestrator.wizard.get_output_filename.return_value = "noise_dataset.fasta"

        mock_dataset = Dataset(["ACGT", "ACTT", "ACGT", "ATGT", "GCGT"])

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetGenerator"
        ) as mock_generator:
            with patch(
                "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetConfig"
            ) as mock_config_class:
                mock_config = Mock()
                mock_config_class.return_value = mock_config
                mock_generator.generate_from_config.return_value = (
                    mock_dataset,
                    mock_params,
                )

                with patch.object(orchestrator.repository, "save") as mock_save:
                    mock_save.return_value = "/path/to/noise_dataset.fasta"

                    with patch("builtins.print"):
                        result = orchestrator._handle_synthetic_generation()

        # Verify config was created with correct parameters
        mock_config_class.assert_called_once_with(**mock_params)
        assert result == "/path/to/noise_dataset.fasta"

    def test_handle_synthetic_generation_error(self):
        """Test error handling in synthetic generation."""
        orchestrator = DatasetGenerationOrchestrator()

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_synthetic_params.return_value = {
            "method": "invalid"
        }
        orchestrator.wizard.generate_default_filename.return_value = "test.fasta"
        orchestrator.wizard.get_output_filename.return_value = "test.fasta"

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetGenerator"
        ) as mock_generator:
            mock_generator.generate_from_config.side_effect = ValueError(
                "Invalid method"
            )

            with patch("builtins.print"):  # Suppress error output
                result = orchestrator._handle_synthetic_generation()

        assert result is None


class TestDatasetGenerationOrchestratorRealGeneration:
    """Tests for real dataset generation workflow."""

    def test_handle_real_generation_ncbi_success(self):
        """Test successful real dataset generation from NCBI."""
        orchestrator = DatasetGenerationOrchestrator()

        mock_params = {
            "source": "ncbi",
            "query": "COX1 AND mitochondrion",
            "max_sequences": 10,
            "min_length": 100,
            "max_length": 200,
        }

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_real_params.return_value = mock_params
        orchestrator.wizard.generate_default_filename.return_value = (
            "real_COX1_AND_mitochondrion.fasta"
        )
        orchestrator.wizard.get_output_filename.return_value = "cox1_dataset.fasta"

        mock_dataset = Dataset(["ACGT" * 50, "ATGT" * 50])  # Simulated sequences
        mock_download_params = {"term": "COX1 AND mitochondrion", "n": 10}

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.EntrezDatasetDownloader"
        ) as mock_downloader:
            mock_downloader.download.return_value = (mock_dataset, mock_download_params)

            with patch.object(orchestrator.repository, "save") as mock_save:
                mock_save.return_value = "/path/to/cox1_dataset.fasta"

                with patch("builtins.print"):
                    result = orchestrator._handle_real_generation()

        assert result == "/path/to/cox1_dataset.fasta"
        mock_downloader.download.assert_called_once()
        mock_save.assert_called_once()

    def test_handle_real_generation_file_success(self):
        """Test successful real dataset generation from file."""
        orchestrator = DatasetGenerationOrchestrator()

        mock_params = {"source": "file", "file_path": "/path/to/existing.fasta"}

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_real_params.return_value = mock_params
        orchestrator.wizard.generate_default_filename.return_value = "real_import.fasta"
        orchestrator.wizard.get_output_filename.return_value = "imported_dataset.fasta"

        mock_dataset = Dataset(["ACGTACGT", "ATGTACGT"])
        mock_load_params = {"file_path": "/path/to/existing.fasta"}

        with patch.object(orchestrator.repository, "load") as mock_load:
            mock_load.return_value = (mock_dataset, mock_load_params)

            with patch.object(orchestrator.repository, "save") as mock_save:
                mock_save.return_value = "/path/to/imported_dataset.fasta"

                with patch("builtins.print"):
                    result = orchestrator._handle_real_generation()

        assert result == "/path/to/imported_dataset.fasta"
        mock_load.assert_called_once_with("/path/to/existing.fasta")
        mock_save.assert_called_once()

    def test_handle_real_generation_ncbi_error(self):
        """Test error handling in NCBI download."""
        orchestrator = DatasetGenerationOrchestrator()

        mock_params = {"source": "ncbi", "query": "invalid query", "max_sequences": 5}

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_real_params.return_value = mock_params
        orchestrator.wizard.generate_default_filename.return_value = "real_query.fasta"
        orchestrator.wizard.get_output_filename.return_value = "query_dataset.fasta"

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.EntrezDatasetDownloader"
        ) as mock_downloader:
            mock_downloader.download.side_effect = ValueError("NCBI error")

            with patch("builtins.print"):  # Suppress error output
                result = orchestrator._handle_real_generation()

        assert result is None

    def test_handle_real_generation_file_not_found(self):
        """Test error handling when source file is not found."""
        orchestrator = DatasetGenerationOrchestrator()

        mock_params = {"source": "file", "file_path": "/nonexistent/file.fasta"}

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_real_params.return_value = mock_params
        orchestrator.wizard.generate_default_filename.return_value = "real_import.fasta"
        orchestrator.wizard.get_output_filename.return_value = "imported.fasta"

        with patch.object(orchestrator.repository, "load") as mock_load:
            mock_load.side_effect = FileNotFoundError("File not found")

            with patch("builtins.print"):  # Suppress error output
                result = orchestrator._handle_real_generation()

        assert result is None

    def test_handle_real_generation_unsupported_source(self):
        """Test handling of unsupported source type."""
        orchestrator = DatasetGenerationOrchestrator()

        mock_params = {"source": "unsupported_source", "some_param": "value"}

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_real_params.return_value = mock_params
        orchestrator.wizard.generate_default_filename.return_value = (
            "real_unknown.fasta"
        )
        orchestrator.wizard.get_output_filename.return_value = "unknown.fasta"

        with patch("builtins.print"):  # Suppress error output
            result = orchestrator._handle_real_generation()

        assert result is None


class TestDatasetGenerationOrchestratorConfigCreation:
    """Tests for configuration object creation."""

    def test_create_synthetic_config_random(self):
        """Test creation of synthetic config for random method."""
        orchestrator = DatasetGenerationOrchestrator()

        params = {
            "method": "random",
            "n": 10,
            "length": 8,
            "alphabet": "ACGT",
            "seed": 42,
        }

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetConfig"
        ) as mock_config:
            orchestrator._create_synthetic_config(params)
            mock_config.assert_called_once_with(**params)

    def test_create_synthetic_config_with_optional_params(self):
        """Test creation of synthetic config with optional parameters."""
        orchestrator = DatasetGenerationOrchestrator()

        params = {
            "method": "noise",
            "n": 5,
            "length": 6,
            "alphabet": "ACGT",
            "noise_rate": 0.3,
            "center_sequence": "ACGTAC",
        }

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetConfig"
        ) as mock_config:
            orchestrator._create_synthetic_config(params)
            mock_config.assert_called_once_with(**params)

    def test_create_entrez_config(self):
        """Test creation of Entrez config."""
        orchestrator = DatasetGenerationOrchestrator()

        params = {
            "query": "COX1 AND mitochondrion",
            "max_sequences": 20,
            "min_length": 100,
            "max_length": 200,
        }

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.EntrezDatasetConfig"
        ) as mock_config:
            orchestrator._create_entrez_config(params)

            # Should map max_sequences to retmax
            expected_params = {
                "query": "COX1 AND mitochondrion",
                "retmax": 20,
                "min_length": 100,
                "max_length": 200,
            }
            mock_config.assert_called_once_with(**expected_params)


class TestDatasetGenerationOrchestratorEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_empty_dataset_handling(self):
        """Test handling of empty datasets."""
        orchestrator = DatasetGenerationOrchestrator()

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_synthetic_params.return_value = {
            "method": "random",
            "n": 0,
            "length": 4,
            "alphabet": "ACGT",
        }
        orchestrator.wizard.generate_default_filename.return_value = "empty.fasta"
        orchestrator.wizard.get_output_filename.return_value = "empty.fasta"

        empty_dataset = Dataset([])

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetGenerator"
        ) as mock_generator:
            mock_generator.generate_from_config.return_value = (empty_dataset, {})

            with patch.object(orchestrator.repository, "save") as mock_save:
                mock_save.return_value = "/path/to/empty.fasta"

                with patch("builtins.print"):
                    result = orchestrator._handle_synthetic_generation()

        assert result == "/path/to/empty.fasta"
        mock_save.assert_called_once()

    def test_very_large_dataset_handling(self):
        """Test handling of very large datasets."""
        orchestrator = DatasetGenerationOrchestrator()

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_synthetic_params.return_value = {
            "method": "random",
            "n": 1000,
            "length": 100,
            "alphabet": "ACGT",
        }
        orchestrator.wizard.generate_default_filename.return_value = "large.fasta"
        orchestrator.wizard.get_output_filename.return_value = "large.fasta"

        # Simulate large dataset
        large_dataset = Dataset(["A" * 100] * 1000)

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetGenerator"
        ) as mock_generator:
            mock_generator.generate_from_config.return_value = (large_dataset, {})

            with patch.object(orchestrator.repository, "save") as mock_save:
                mock_save.return_value = "/path/to/large.fasta"

                with patch("builtins.print"):
                    result = orchestrator._handle_synthetic_generation()

        assert result == "/path/to/large.fasta"

    def test_special_characters_in_filename(self):
        """Test handling of special characters in filenames."""
        orchestrator = DatasetGenerationOrchestrator()

        orchestrator.wizard = Mock()
        orchestrator.wizard.collect_synthetic_params.return_value = {
            "method": "random",
            "n": 5,
            "length": 4,
            "alphabet": "ACGT",
        }
        orchestrator.wizard.generate_default_filename.return_value = "test*file?.fasta"
        orchestrator.wizard.get_output_filename.return_value = "test*file?.fasta"

        mock_dataset = Dataset(["ACGT"] * 5)

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetGenerator"
        ) as mock_generator:
            mock_generator.generate_from_config.return_value = (mock_dataset, {})

            with patch.object(orchestrator.repository, "save") as mock_save:
                # Repository should handle special characters appropriately
                mock_save.return_value = "/path/to/test_file_.fasta"

                with patch("builtins.print"):
                    result = orchestrator._handle_synthetic_generation()

        assert result == "/path/to/test_file_.fasta"


class TestDatasetGenerationOrchestratorIntegration:
    """Integration-style tests for complete workflows."""

    def test_complete_synthetic_workflow_integration(self):
        """Test complete synthetic workflow integration."""
        orchestrator = DatasetGenerationOrchestrator()

        # Mock complete workflow
        orchestrator.wizard = Mock()

        # Sequence of wizard calls
        orchestrator.wizard.show_main_menu.return_value = "synthetic"
        orchestrator.wizard.collect_synthetic_params.return_value = {
            "method": "random",
            "n": 5,
            "length": 8,
            "alphabet": "ACGT",
            "seed": 42,
        }
        orchestrator.wizard.generate_default_filename.return_value = (
            "synthetic_random_n5_L8.fasta"
        )
        orchestrator.wizard.get_output_filename.return_value = "my_dataset.fasta"

        # Mock generation
        mock_dataset = Dataset(
            ["ACGTACGT", "ATGTACGT", "GCGTACGT", "ACGTCCGT", "ATGTCCGT"]
        )

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.SyntheticDatasetGenerator"
        ) as mock_generator:
            mock_generator.generate_from_config.return_value = (mock_dataset, {})

            with patch.object(orchestrator.repository, "save") as mock_save:
                mock_save.return_value = "/datasets/my_dataset.fasta"

                with patch("builtins.print"):
                    result = orchestrator.run_interactive_generation()

        assert result == "/datasets/my_dataset.fasta"

        # Verify complete workflow was executed
        orchestrator.wizard.show_main_menu.assert_called_once()
        orchestrator.wizard.collect_synthetic_params.assert_called_once()
        orchestrator.wizard.generate_default_filename.assert_called_once()
        orchestrator.wizard.get_output_filename.assert_called_once()
        mock_generator.generate_from_config.assert_called_once()
        mock_save.assert_called_once()

    def test_complete_real_workflow_integration(self):
        """Test complete real dataset workflow integration."""
        orchestrator = DatasetGenerationOrchestrator()

        orchestrator.wizard = Mock()

        # Sequence of wizard calls
        orchestrator.wizard.show_main_menu.return_value = "real"
        orchestrator.wizard.collect_real_params.return_value = {
            "source": "ncbi",
            "query": "COX1",
            "max_sequences": 10,
            "min_length": 50,
            "max_length": 100,
        }
        orchestrator.wizard.generate_default_filename.return_value = "real_COX1.fasta"
        orchestrator.wizard.get_output_filename.return_value = "cox1_seqs.fasta"

        # Mock download
        mock_dataset = Dataset(["A" * 75, "T" * 75, "G" * 75])

        with patch(
            "src.infrastructure.orchestrators.dataset_generation_orchestrator.EntrezDatasetDownloader"
        ) as mock_downloader:
            mock_downloader.download.return_value = (mock_dataset, {})

            with patch.object(orchestrator.repository, "save") as mock_save:
                mock_save.return_value = "/datasets/cox1_seqs.fasta"

                with patch("builtins.print"):
                    result = orchestrator.run_interactive_generation()

        assert result == "/datasets/cox1_seqs.fasta"

        # Verify workflow
        orchestrator.wizard.show_main_menu.assert_called_once()
        orchestrator.wizard.collect_real_params.assert_called_once()
        mock_downloader.download.assert_called_once()
        mock_save.assert_called_once()
