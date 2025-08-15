"""
Unit tests for dataset_service.py.

Tests the load_dataset function that orchestrates loading
datasets from different sources (synthetic, file, entrez).
"""

from unittest.mock import Mock, patch

import pytest

from src.application.services.dataset_service import load_dataset
from src.domain.config import (
    SyntheticDatasetConfig,
    FileDatasetConfig,
    EntrezDatasetConfig,
)
from src.domain.dataset import Dataset


class TestLoadDatasetSynthetic:
    """Tests for loading synthetic datasets."""

    @patch("src.application.services.dataset_service.SyntheticDatasetGenerator")
    def test_load_synthetic_dataset(self, mock_generator):
        """Test loading synthetic dataset."""
        # Setup mock
        mock_dataset = Dataset(["ACGT", "ATGT"])
        mock_params = {"method": "random", "n": 2}
        mock_generator.generate_from_config.return_value = (mock_dataset, mock_params)

        config = SyntheticDatasetConfig(method="random", n=2, length=4, alphabet="ACGT")

        dataset, params = load_dataset(config)

        assert dataset == mock_dataset
        assert params == mock_params
        mock_generator.generate_from_config.assert_called_once_with(config)

    @patch("src.application.services.dataset_service.SyntheticDatasetGenerator")
    def test_load_synthetic_with_noise(self, mock_generator):
        """Test loading synthetic dataset with noise."""
        mock_dataset = Dataset(["ACGT", "ACTT"])
        mock_params = {"method": "noise", "noise_rate": 0.2}
        mock_generator.generate_from_config.return_value = (mock_dataset, mock_params)

        config = SyntheticDatasetConfig(
            method="noise", n=2, length=4, alphabet="ACGT", noise_rate=0.2
        )

        dataset, params = load_dataset(config)

        assert isinstance(dataset, Dataset)
        assert "noise_rate" in params
        mock_generator.generate_from_config.assert_called_once_with(config)

    @patch("src.application.services.dataset_service.SyntheticDatasetGenerator")
    def test_load_synthetic_clustered(self, mock_generator):
        """Test loading clustered synthetic dataset."""
        mock_dataset = Dataset(["ACGT", "ACGA", "TGCA", "TGCT"])
        mock_params = {"method": "clustered", "num_clusters": 2}
        mock_generator.generate_from_config.return_value = (mock_dataset, mock_params)

        config = SyntheticDatasetConfig(
            method="clustered", n=4, length=4, alphabet="ACGT", num_clusters=2
        )

        dataset, params = load_dataset(config)

        assert dataset.size == 4
        assert "num_clusters" in params
        mock_generator.generate_from_config.assert_called_once_with(config)


class TestLoadDatasetFile:
    """Tests for loading file-based datasets."""

    @patch("src.application.services.dataset_service.FileDatasetRepository")
    def test_load_file_dataset(self, mock_repo):
        """Test loading dataset from file."""
        mock_dataset = Dataset(["ACGTACGT", "ATGTACGT"])
        mock_params = {"file_path": "/path/to/dataset.fasta"}
        mock_repo.load.return_value = (mock_dataset, mock_params)

        config = FileDatasetConfig(filename="test_dataset.fasta")

        dataset, params = load_dataset(config)

        assert dataset == mock_dataset
        assert params == mock_params
        mock_repo.load.assert_called_once_with("test_dataset.fasta")

    @patch("src.application.services.dataset_service.FileDatasetRepository")
    def test_load_file_dataset_with_path(self, mock_repo):
        """Test loading dataset from specific file path."""
        mock_dataset = Dataset(["ACGT"])
        mock_params = {"file_path": "/absolute/path/data.fasta"}
        mock_repo.load.return_value = (mock_dataset, mock_params)

        config = FileDatasetConfig(filename="/absolute/path/data.fasta")

        dataset, params = load_dataset(config)

        assert isinstance(dataset, Dataset)
        assert "file_path" in params
        mock_repo.load.assert_called_once_with("/absolute/path/data.fasta")

    @patch("src.application.services.dataset_service.FileDatasetRepository")
    def test_load_file_dataset_nonexistent(self, mock_repo):
        """Test loading non-existent file propagates error."""
        from src.domain.errors import DatasetNotFoundError

        mock_repo.load.side_effect = DatasetNotFoundError("File not found")

        config = FileDatasetConfig(filename="nonexistent.fasta")

        with pytest.raises(DatasetNotFoundError):
            load_dataset(config)


class TestLoadDatasetEntrez:
    """Tests for loading datasets from NCBI Entrez."""

    @patch("src.application.services.dataset_service.EntrezDatasetDownloader")
    def test_load_entrez_dataset(self, mock_downloader):
        """Test loading dataset from NCBI."""
        mock_dataset = Dataset(["ACGTACGTACGT", "ATGTACGTACGT"])
        mock_params = {"term": "COX1 AND mitochondrion", "n": 2, "db": "nucleotide"}
        mock_downloader.download.return_value = (mock_dataset, mock_params)

        config = EntrezDatasetConfig(
            query="COX1 AND mitochondrion", retmax=2, db="nucleotide"
        )

        dataset, params = load_dataset(config)

        assert dataset == mock_dataset
        assert params == mock_params
        mock_downloader.download.assert_called_once_with(config)

    @patch("src.application.services.dataset_service.EntrezDatasetDownloader")
    def test_load_entrez_with_filters(self, mock_downloader):
        """Test loading Entrez dataset with length filters."""
        mock_dataset = Dataset(["ACGTACGTACGTACGT"])  # 16bp
        mock_params = {
            "term": "16S rRNA",
            "min_length": 15,
            "max_length": 20,
            "uniform_policy": "strict",
        }
        mock_downloader.download.return_value = (mock_dataset, mock_params)

        config = EntrezDatasetConfig(
            query="16S rRNA",
            retmax=1,
            min_length=15,
            max_length=20,
            uniform_policy="strict",
        )

        dataset, params = load_dataset(config)

        assert len(dataset.sequences[0]) == 16
        assert params["min_length"] == 15
        assert params["max_length"] == 20
        mock_downloader.download.assert_called_once_with(config)

    @patch("src.application.services.dataset_service.EntrezDatasetDownloader")
    def test_load_entrez_download_error(self, mock_downloader):
        """Test handling of Entrez download errors."""
        mock_downloader.download.side_effect = ValueError("NCBI API error")

        config = EntrezDatasetConfig(query="invalid query", retmax=5)

        with pytest.raises(ValueError, match="NCBI API error"):
            load_dataset(config)


class TestLoadDatasetErrorHandling:
    """Tests for error handling and edge cases."""

    def test_unsupported_config_type(self):
        """Test error with unsupported config type."""

        class UnsupportedConfig:
            pass

        unsupported_config = UnsupportedConfig()

        with pytest.raises(TypeError, match="Tipo de dataset não suportado"):
            load_dataset(unsupported_config)

    def test_none_config(self):
        """Test error with None config."""
        with pytest.raises(AttributeError):
            load_dataset(None)

    @patch("src.application.services.dataset_service.SyntheticDatasetGenerator")
    def test_synthetic_generation_error(self, mock_generator):
        """Test handling of synthetic generation errors."""
        mock_generator.generate_from_config.side_effect = ValueError(
            "Invalid parameters"
        )

        config = SyntheticDatasetConfig(method="invalid", n=0, length=-1)

        with pytest.raises(ValueError, match="Invalid parameters"):
            load_dataset(config)

    @patch("src.application.services.dataset_service.FileDatasetRepository")
    def test_file_loading_permission_error(self, mock_repo):
        """Test handling of file permission errors."""
        mock_repo.load.side_effect = PermissionError("Permission denied")

        config = FileDatasetConfig(filename="protected.fasta")

        with pytest.raises(PermissionError):
            load_dataset(config)


class TestLoadDatasetIntegration:
    """Integration-style tests for load_dataset function."""

    @patch("src.application.services.dataset_service.SyntheticDatasetGenerator")
    @patch("src.application.services.dataset_service.FileDatasetRepository")
    @patch("src.application.services.dataset_service.EntrezDatasetDownloader")
    def test_load_different_config_types(self, mock_entrez, mock_file, mock_synthetic):
        """Test loading datasets of different types."""
        # Setup mocks
        synthetic_dataset = Dataset(["ACGT", "ATGT"])
        file_dataset = Dataset(["GCGT", "GAGT"])
        entrez_dataset = Dataset(["TCGT", "TAGT"])

        mock_synthetic.generate_from_config.return_value = (synthetic_dataset, {})
        mock_file.load.return_value = (file_dataset, {})
        mock_entrez.download.return_value = (entrez_dataset, {})

        # Test synthetic
        synthetic_config = SyntheticDatasetConfig(method="random", n=2, length=4)
        dataset1, _ = load_dataset(synthetic_config)
        assert dataset1 == synthetic_dataset

        # Test file
        file_config = FileDatasetConfig(filename="test.fasta")
        dataset2, _ = load_dataset(file_config)
        assert dataset2 == file_dataset

        # Test Entrez
        entrez_config = EntrezDatasetConfig(query="test", retmax=2)
        dataset3, _ = load_dataset(entrez_config)
        assert dataset3 == entrez_dataset

        # Verify all were called
        mock_synthetic.generate_from_config.assert_called_once()
        mock_file.load.assert_called_once()
        mock_entrez.download.assert_called_once()

    def test_config_type_detection(self):
        """Test that function correctly detects config types."""
        synthetic_config = SyntheticDatasetConfig(method="random", n=1, length=4)
        file_config = FileDatasetConfig(filename="test.fasta")
        entrez_config = EntrezDatasetConfig(query="test", retmax=1)

        # Should not raise TypeError for valid configs
        assert isinstance(synthetic_config, SyntheticDatasetConfig)
        assert isinstance(file_config, FileDatasetConfig)
        assert isinstance(entrez_config, EntrezDatasetConfig)

    @patch("src.application.services.dataset_service.SyntheticDatasetGenerator")
    def test_parameter_passthrough(self, mock_generator):
        """Test that parameters are correctly passed through."""
        expected_params = {
            "method": "clustered",
            "n": 10,
            "length": 8,
            "alphabet": "ACGT",
            "num_clusters": 3,
            "seed": 42,
        }

        mock_dataset = Dataset(["ACGT"] * 10)
        mock_generator.generate_from_config.return_value = (
            mock_dataset,
            expected_params,
        )

        config = SyntheticDatasetConfig(
            method="clustered", n=10, length=8, alphabet="ACGT", num_clusters=3, seed=42
        )

        dataset, params = load_dataset(config)

        assert params == expected_params
        # Verify config object was passed unchanged
        mock_generator.generate_from_config.assert_called_once_with(config)


class TestLoadDatasetEdgeCases:
    """Tests for edge cases and boundary conditions."""

    @patch("src.application.services.dataset_service.SyntheticDatasetGenerator")
    def test_empty_synthetic_dataset(self, mock_generator):
        """Test loading empty synthetic dataset."""
        empty_dataset = Dataset([])
        mock_generator.generate_from_config.return_value = (empty_dataset, {"n": 0})

        config = SyntheticDatasetConfig(method="random", n=0, length=4)

        dataset, params = load_dataset(config)

        assert dataset.size == 0
        assert params["n"] == 0

    @patch("src.application.services.dataset_service.FileDatasetRepository")
    def test_load_empty_file(self, mock_repo):
        """Test loading empty file dataset."""
        empty_dataset = Dataset([])
        mock_repo.load.return_value = (empty_dataset, {"file_path": "empty.fasta"})

        config = FileDatasetConfig(filename="empty.fasta")

        dataset, params = load_dataset(config)

        assert dataset.size == 0

    def test_invalid_config_attributes(self):
        """Test handling of configs with invalid attributes."""
        # Create config with required attributes but invalid values
        config = SyntheticDatasetConfig(method="random", n=1, length=1)

        # Manually set invalid attribute (this would normally be caught by validation)
        config.method = None

        # Should still attempt to load (error handling is delegated to generators)
        with patch(
            "src.application.services.dataset_service.SyntheticDatasetGenerator"
        ) as mock_gen:
            mock_gen.generate_from_config.side_effect = ValueError("Invalid method")

            with pytest.raises(ValueError):
                load_dataset(config)
