"""Tests for DatasetFactory module."""

import os
import tempfile
from unittest.mock import patch

import pytest

from src.core.data.dataset_factory import DatasetError, DatasetFactory, DatasetType


class TestDatasetFactory:
    """Test DatasetFactory functionality."""

    def setup_method(self):
        """Setup test fixtures."""
        self.factory = DatasetFactory()

    def test_init(self):
        """Test DatasetFactory initialization."""
        factory = DatasetFactory()
        assert factory is not None
        assert hasattr(factory, "create_dataset")
        assert hasattr(factory, "validate_config")
        assert factory.cache_enabled is True

    def test_init_with_cache_disabled(self):
        """Test DatasetFactory initialization with cache disabled."""
        factory = DatasetFactory(cache_enabled=False)
        assert factory.cache_enabled is False

    def test_create_dataset_invalid_config_type(self):
        """Test creating dataset with invalid config type."""
        with pytest.raises(DatasetError, match="deve ter campo 'tipo'"):
            self.factory.create_dataset({"invalid": "config"})

    def test_create_dataset_missing_tipo(self):
        """Test creating dataset with missing tipo field."""
        config = {"n": 5, "L": 4}
        with pytest.raises(DatasetError, match="Configuração de dataset deve ter campo 'tipo'"):
            self.factory.create_dataset(config)

    def test_create_dataset_unknown_type(self):
        """Test creating dataset with unknown type."""
        config = {"tipo": "unknown"}
        with pytest.raises(DatasetError, match="Tipo de dataset não suportado"):
            self.factory.create_dataset(config)

    @patch("src.core.data.dataset_factory.generate_dataset_with_params")
    def test_create_dataset_synthetic(self, mock_generate):
        """Test creating synthetic dataset."""
        mock_generate.return_value = (
            ["ATCG", "GCTA"],
            {"n": 2, "L": 4, "alphabet": "ATCG"},
        )

        config = {"tipo": "synthetic", "n": 2, "L": 4, "alphabet": "ATCG"}
        sequences, metadata = self.factory.create_dataset(config)

        assert len(sequences) == 2
        assert "ATCG" in sequences
        assert "GCTA" in sequences
        assert metadata["n"] == 2
        assert metadata["L"] == 4
        assert metadata["dataset_type"] == "synthetic"
        mock_generate.assert_called_once()

    @patch("src.core.data.dataset_factory.load_dataset_with_params")
    def test_create_dataset_file(self, mock_load):
        """Test creating dataset from file."""
        mock_load.return_value = (
            ["ATCG", "GCTA"],
            {"n": 2, "L": 4, "filepath": "test.fasta"},
        )

        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\nATCG\n>seq2\nGCTA\n")
            temp_path = f.name

        try:
            config = {"tipo": "file", "filepath": temp_path}
            sequences, metadata = self.factory.create_dataset(config)

            assert len(sequences) == 2
            assert "ATCG" in sequences
            assert "GCTA" in sequences
            assert metadata["dataset_type"] == "file"
            mock_load.assert_called_once()
        finally:
            os.unlink(temp_path)

    def test_create_dataset_file_missing_filepath(self):
        """Test creating file dataset with missing filepath."""
        config = {"tipo": "file"}
        with pytest.raises(DatasetError, match="Dataset tipo 'file' deve ter campo 'filepath'"):
            self.factory.create_dataset(config)

    def test_create_dataset_file_nonexistent(self):
        """Test creating dataset from nonexistent file."""
        config = {"tipo": "file", "filepath": "/nonexistent/path.fasta"}
        with pytest.raises(DatasetError, match="Arquivo não encontrado"):
            self.factory.create_dataset(config)

    @patch("src.core.data.dataset_factory.fetch_dataset_silent")
    def test_create_dataset_entrez(self, mock_fetch):
        """Test creating dataset from Entrez."""
        mock_fetch.return_value = (
            ["ATCG", "GCTA"],
            {"n": 2, "L": 4, "db": "nucleotide"},
        )

        config = {
            "tipo": "entrez",
            "email": "test@example.com",
            "db": "nucleotide",
            "term": "bacteria",
        }
        sequences, metadata = self.factory.create_dataset(config)

        assert len(sequences) == 2
        assert "ATCG" in sequences
        assert "GCTA" in sequences
        assert metadata["dataset_type"] == "entrez"
        mock_fetch.assert_called_once()

    def test_create_dataset_entrez_missing_fields(self):
        """Test creating entrez dataset with missing fields."""
        config = {"tipo": "entrez", "email": "test@example.com"}
        with pytest.raises(DatasetError, match="Dataset tipo 'entrez' deve ter campos"):
            self.factory.create_dataset(config)

    def test_create_dataset_synthetic_missing_fields(self):
        """Test creating synthetic dataset with missing fields."""
        config = {"tipo": "synthetic", "n": 10}
        with pytest.raises(DatasetError, match="Dataset tipo 'synthetic' deve ter campos"):
            self.factory.create_dataset(config)

    @patch("src.core.data.dataset_factory.generate_dataset_with_params")
    def test_create_dataset_with_caching(self, mock_generate):
        """Test dataset creation with caching."""
        mock_generate.return_value = (
            ["ATCG", "GCTA"],
            {"n": 2, "L": 4, "alphabet": "ATCG"},
        )

        config = {"tipo": "synthetic", "n": 2, "L": 4, "alphabet": "ATCG", "seed": 42}

        # First call should create and cache
        sequences1, metadata1 = self.factory.create_dataset(config)
        assert len(sequences1) == 2
        assert mock_generate.call_count == 1

        # Second call should return cached result
        sequences2, metadata2 = self.factory.create_dataset(config)
        assert sequences1 == sequences2
        assert metadata1 == metadata2
        # Still only one call to generate function
        assert mock_generate.call_count == 1

    @patch("src.core.data.dataset_factory.generate_dataset_with_params")
    def test_create_dataset_with_base_index(self, mock_generate):
        """Test dataset creation with base index."""
        mock_generate.return_value = (
            ["ATCG"],
            {"n": 1, "L": 4, "alphabet": "ATCG", "seed": 44},
        )

        config = {"tipo": "synthetic", "n": 1, "L": 4, "alphabet": "ATCG", "seed": 42}
        sequences, metadata = self.factory.create_dataset(config, base_index=2)

        assert len(sequences) == 1
        assert metadata["base_index"] == 2
        # Verify seed was modified by base index
        call_args = mock_generate.call_args[0][0]
        assert call_args["seed"] == 44  # 42 + 2

    def test_clear_cache(self):
        """Test clearing the cache."""
        # Create a dataset to populate cache
        with patch("src.core.data.dataset_factory.generate_dataset_with_params") as mock_generate:
            mock_generate.return_value = (
                ["ATCG"],
                {"n": 1, "L": 4, "alphabet": "ATCG"},
            )

            config = {
                "tipo": "synthetic",
                "n": 1,
                "L": 4,
                "alphabet": "ATCG",
                "seed": 42,
            }
            self.factory.create_dataset(config)

            # Cache should have entries
            assert len(self.factory._cache) > 0

            # Clear cache
            self.factory.clear_cache()

            # Cache should be empty
            assert len(self.factory._cache) == 0

    def test_get_cache_info(self):
        """Test getting cache information."""
        # Initially empty
        cache_info = self.factory.get_cache_info()
        assert cache_info["enabled"] is True
        assert cache_info["size"] == 0
        assert cache_info["keys"] == []

        # Create dataset to populate cache
        with patch("src.core.data.dataset_factory.generate_dataset_with_params") as mock_generate:
            mock_generate.return_value = (
                ["ATCG"],
                {"n": 1, "L": 4, "alphabet": "ATCG"},
            )

            config = {"tipo": "synthetic", "n": 1, "L": 4, "alphabet": "ATCG"}
            self.factory.create_dataset(config)

            # Cache should have entries
            cache_info = self.factory.get_cache_info()
            assert cache_info["enabled"] is True
            assert cache_info["size"] == 1
            assert len(cache_info["keys"]) == 1

    def test_validate_config_valid(self):
        """Test validating valid configs."""
        valid_configs = [
            {"tipo": "synthetic", "n": 5, "L": 4, "alphabet": "ATCG"},
            {"tipo": "file", "filepath": "/tmp/test.fasta"},
            {
                "tipo": "entrez",
                "email": "test@example.com",
                "db": "nucleotide",
                "term": "bacteria",
            },
        ]

        for config in valid_configs:
            # Should not raise exception
            self.factory.validate_config(config)

    def test_validate_config_invalid_type(self):
        """Test validating config with invalid type."""
        with pytest.raises(DatasetError, match="deve ter campo 'tipo'"):
            self.factory.validate_config({"invalid": "config"})

    def test_validate_config_missing_tipo(self):
        """Test validating config with missing tipo."""
        config = {"n": 5, "L": 4}
        with pytest.raises(DatasetError, match="Configuração de dataset deve ter campo 'tipo'"):
            self.factory.validate_config(config)

    def test_validate_config_unknown_tipo(self):
        """Test validating config with unknown tipo."""
        config = {"tipo": "unknown"}
        with pytest.raises(DatasetError, match="Tipo de dataset não suportado"):
            self.factory.validate_config(config)

    def test_validate_config_file_missing_filepath(self):
        """Test validating file config with missing filepath."""
        config = {"tipo": "file"}
        with pytest.raises(DatasetError, match="Dataset tipo 'file' deve ter campo 'filepath'"):
            self.factory.validate_config(config)

    def test_validate_config_entrez_missing_fields(self):
        """Test validating entrez config with missing fields."""
        config = {"tipo": "entrez", "email": "test@example.com"}
        with pytest.raises(DatasetError, match="Dataset tipo 'entrez' deve ter campos"):
            self.factory.validate_config(config)

    def test_validate_config_synthetic_missing_fields(self):
        """Test validating synthetic config with missing fields."""
        config = {"tipo": "synthetic", "n": 5}
        with pytest.raises(DatasetError, match="Dataset tipo 'synthetic' deve ter campos"):
            self.factory.validate_config(config)

    @patch("src.core.data.dataset_factory.generate_dataset_with_params")
    def test_create_multiple_datasets(self, mock_generate):
        """Test creating multiple datasets."""

        # Mock should return different base_index for each call
        def mock_side_effect(config, base_index=0):
            return (
                ["ATCG"],
                {"n": 1, "L": 4, "alphabet": "ATCG", "base_index": base_index},
            )

        mock_generate.side_effect = mock_side_effect

        config = {"tipo": "synthetic", "n": 1, "L": 4, "alphabet": "ATCG"}
        datasets = self.factory.create_multiple_datasets(config, num_bases=3)

        assert len(datasets) == 3
        for sequences, metadata in datasets:
            assert len(sequences) == 1
            assert metadata["dataset_type"] == "synthetic"
        # Each base should have different base_index
        assert datasets[0][1]["base_index"] == 0
        assert datasets[1][1]["base_index"] == 1
        assert datasets[2][1]["base_index"] == 2

    def test_create_multiple_datasets_zero_bases(self):
        """Test creating multiple datasets with zero bases."""
        config = {"tipo": "synthetic", "n": 1, "L": 4, "alphabet": "ATCG"}
        datasets = self.factory.create_multiple_datasets(config, num_bases=0)
        assert len(datasets) == 0

    def test_get_dataset_info_valid_file(self):
        """Test getting dataset info for valid file config."""
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\nATCG\n")
            temp_path = f.name

        try:
            config = {"tipo": "file", "filepath": temp_path}
            info = self.factory.get_dataset_info(config)

            assert info["tipo"] == "file"
            assert info["valid"] is True
            assert info["error"] is None
            assert info["filepath"] == temp_path
            assert info["exists"] is True
            assert info["size"] > 0
        finally:
            os.unlink(temp_path)

    def test_get_dataset_info_nonexistent_file(self):
        """Test getting dataset info for nonexistent file."""
        config = {"tipo": "file", "filepath": "/nonexistent/file.fasta"}
        info = self.factory.get_dataset_info(config)

        assert info["tipo"] == "file"
        assert info["valid"] is True  # Config is valid, but file doesn't exist
        assert info["filepath"] == "/nonexistent/file.fasta"
        assert info["exists"] is False
        assert info["size"] == 0

    def test_get_dataset_info_entrez(self):
        """Test getting dataset info for entrez config."""
        config = {
            "tipo": "entrez",
            "email": "test@example.com",
            "db": "nucleotide",
            "term": "bacteria",
            "n": 10,
        }
        info = self.factory.get_dataset_info(config)

        assert info["tipo"] == "entrez"
        assert info["valid"] is True
        assert info["error"] is None
        assert info["db"] == "nucleotide"
        assert info["term"] == "bacteria"
        assert info["n"] == 10

    def test_get_dataset_info_synthetic(self):
        """Test getting dataset info for synthetic config."""
        config = {"tipo": "synthetic", "n": 5, "L": 4, "alphabet": "ATCG", "noise": 0.1}
        info = self.factory.get_dataset_info(config)

        assert info["tipo"] == "synthetic"
        assert info["valid"] is True
        assert info["error"] is None
        assert info["n"] == 5
        assert info["L"] == 4
        assert info["alphabet"] == "ATCG"
        assert info["noise"] == 0.1

    def test_get_dataset_info_invalid_config(self):
        """Test getting dataset info for invalid config."""
        config = {"tipo": "synthetic", "n": 5}  # Missing L and alphabet
        info = self.factory.get_dataset_info(config)

        assert info["tipo"] == "synthetic"
        assert info["valid"] is False
        assert info["error"] is not None
        assert "deve ter campos" in info["error"]

    def test_get_cache_key_synthetic(self):
        """Test cache key generation for synthetic dataset."""
        config = {"tipo": "synthetic", "n": 5, "L": 4, "alphabet": "ATCG", "seed": 42}
        key = self.factory._get_cache_key(config, base_index=1)

        assert isinstance(key, str)
        assert "synthetic" in key
        assert "1" in key  # base_index
        assert "5" in key  # n
        assert "4" in key  # L
        assert "ATCG" in key  # alphabet
        assert "42" in key  # seed

    def test_get_cache_key_file(self):
        """Test cache key generation for file dataset."""
        config = {"tipo": "file", "filepath": "/tmp/test.fasta"}
        key = self.factory._get_cache_key(config, base_index=0)

        assert isinstance(key, str)
        assert "file" in key
        assert "0" in key  # base_index
        assert "/tmp/test.fasta" in key

    def test_get_cache_key_entrez(self):
        """Test cache key generation for entrez dataset."""
        config = {
            "tipo": "entrez",
            "email": "test@example.com",
            "db": "nucleotide",
            "term": "bacteria",
            "n": 10,
            "seed": 42,
        }
        key = self.factory._get_cache_key(config, base_index=2)

        assert isinstance(key, str)
        assert "entrez" in key
        assert "2" in key  # base_index
        assert "nucleotide" in key
        assert "bacteria" in key
        assert "10" in key  # n
        assert "42" in key  # seed

    def test_dataset_type_enum(self):
        """Test DatasetType enum values."""
        assert DatasetType.SYNTHETIC.value == "synthetic"
        assert DatasetType.FILE.value == "file"
        assert DatasetType.ENTREZ.value == "entrez"
