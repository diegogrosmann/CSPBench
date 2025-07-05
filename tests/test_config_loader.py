"""Tests for ConfigLoader module."""

import os
import tempfile
from pathlib import Path

import pytest
import yaml

from src.core.config.config_loader import ConfigError, ConfigLoader


class TestConfigLoader:
    """Test ConfigLoader functionality."""

    def setup_method(self):
        """Setup test fixtures."""
        self.config_loader = ConfigLoader()

        # Sample valid config aligned with real interface
        self.valid_config = {
            "nome": "teste",
            "algoritmos": ["Baseline", "BLF-GA"],
            "dataset": {"tipo": "file", "filepath": "test.fasta"},
            "execucoes_por_algoritmo_por_base": 3,
            "num_bases": 1,
            "timeout": 300,
        }

        # Sample invalid config
        self.invalid_config = {"nome": "teste", "algoritmos": []}

    def test_init(self):
        """Test ConfigLoader initialization."""
        loader = ConfigLoader()
        assert loader is not None
        assert hasattr(loader, "load_config")
        assert hasattr(loader, "_validate_config")  # Private method

    def test_load_config_valid_file(self):
        """Test loading a valid config file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(self.valid_config, f)
            temp_path = f.name

        try:
            config = self.config_loader.load_config(temp_path)
            assert config is not None
            assert config["nome"] == "teste"
            assert config["algoritmos"] == ["Baseline", "BLF-GA"]
            assert config["dataset"]["tipo"] == "file"
            # The filepath gets resolved to absolute path
            assert config["dataset"]["filepath"].endswith("test.fasta")
            assert config["execucoes_por_algoritmo_por_base"] == 3
            assert config["num_bases"] == 1
        finally:
            os.unlink(temp_path)

    def test_load_config_nonexistent_file(self):
        """Test loading a nonexistent config file."""
        with pytest.raises(ConfigError):
            self.config_loader.load_config("/nonexistent/path.yaml")

    def test_load_config_invalid_yaml(self):
        """Test loading invalid YAML."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("invalid: yaml: content: [")
            temp_path = f.name

        try:
            with pytest.raises(ConfigError):
                self.config_loader.load_config(temp_path)
        finally:
            os.unlink(temp_path)

    def test_load_multiple_configs(self):
        """Test loading multiple config files."""
        config1 = {
            "nome": "config1",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "file", "filepath": "test1.fasta"},
            "execucoes_por_algoritmo_por_base": 1,
            "num_bases": 1,
        }
        config2 = {
            "nome": "config2",
            "algoritmos": ["BLF-GA"],
            "dataset": {"tipo": "file", "filepath": "test2.fasta"},
            "execucoes_por_algoritmo_por_base": 2,
            "num_bases": 1,
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f1:
            yaml.dump(config1, f1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f2:
            yaml.dump(config2, f2)
            path2 = f2.name

        try:
            configs = self.config_loader.load_multiple_configs([path1, path2])
            assert len(configs) == 2
            assert configs[0]["nome"] == "config1"
            assert configs[0]["algoritmos"] == ["Baseline"]
            assert configs[1]["nome"] == "config2"
            assert configs[1]["algoritmos"] == ["BLF-GA"]
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_merge_configs(self):
        """Test merging multiple configs."""
        config1 = {
            "nome": "base",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "file", "filepath": "base.fasta"},
            "execucoes_por_algoritmo_por_base": 1,
            "num_bases": 1,
            "timeout": 300,
        }
        config2 = {"nome": "override", "timeout": 600, "max_workers": 8}

        merged = self.config_loader.merge_configs(config1, config2)
        assert merged["nome"] == "override"
        assert merged["timeout"] == 600
        assert merged["max_workers"] == 8
        assert merged["algoritmos"] == ["Baseline"]
        assert merged["dataset"]["tipo"] == "file"

    def test_validate_config_valid(self):
        """Test validating a valid config."""
        # Should not raise any exception
        self.config_loader._validate_config(self.valid_config)

    def test_validate_config_missing_dataset(self):
        """Test validating config with missing dataset."""
        invalid_config = {
            "nome": "teste",
            "algoritmos": ["Baseline"],
            "execucoes_por_algoritmo_por_base": 1,
            "num_bases": 1,
        }
        with pytest.raises(ConfigError, match="Campos obrigatórios ausentes"):
            self.config_loader._validate_config(invalid_config)

    def test_validate_config_missing_algorithms(self):
        """Test validating config with missing algorithms."""
        invalid_config = {
            "nome": "teste",
            "dataset": {"tipo": "file", "filepath": "test.fasta"},
            "execucoes_por_algoritmo_por_base": 1,
            "num_bases": 1,
        }
        with pytest.raises(ConfigError, match="Campos obrigatórios ausentes"):
            self.config_loader._validate_config(invalid_config)

    def test_validate_algorithms_empty(self):
        """Test validating empty algorithms list."""
        with pytest.raises(ConfigError, match="Lista de algoritmos não pode estar vazia"):
            self.config_loader._validate_algorithms([])

    def test_validate_algorithms_invalid_type(self):
        """Test validating non-string algorithm."""
        with pytest.raises(ConfigError, match="Nomes de algoritmos devem ser strings"):
            self.config_loader._validate_algorithms(["Baseline", 123])

    def test_validate_dataset_missing_type(self):
        """Test validating dataset with missing type."""
        dataset = {"filepath": "test.fasta"}
        with pytest.raises(ConfigError, match="Dataset deve ter campo 'tipo'"):
            self.config_loader._validate_dataset(dataset)

    def test_validate_dataset_unknown_type(self):
        """Test validating dataset with unknown type."""
        dataset = {"tipo": "unknown"}
        with pytest.raises(ConfigError, match="Tipo de dataset deve ser um de"):
            self.config_loader._validate_dataset(dataset)

    def test_validate_dataset_synthetic(self):
        """Test validating synthetic dataset."""
        dataset = {"tipo": "synthetic", "n": 10, "L": 4, "alphabet": "ATCG"}
        # Should not raise any exception
        self.config_loader._validate_dataset(dataset)

    def test_validate_dataset_file(self):
        """Test validating file dataset."""
        dataset = {"tipo": "file", "filepath": "test.fasta"}
        # Should not raise any exception
        self.config_loader._validate_dataset(dataset)

    def test_validate_dataset_file_missing_filepath(self):
        """Test validating file dataset with missing filepath."""
        dataset = {"tipo": "file"}
        with pytest.raises(ConfigError, match="Dataset tipo 'file' deve ter campo 'filepath'"):
            self.config_loader._validate_dataset(dataset)

    def test_validate_dataset_entrez(self):
        """Test validating entrez dataset."""
        dataset = {"tipo": "entrez", "email": "test@example.com", "db": "nucleotide", "term": "bacteria"}
        # Should not raise any exception
        self.config_loader._validate_dataset(dataset)

    def test_validate_dataset_entrez_missing_fields(self):
        """Test validating entrez dataset with missing fields."""
        dataset = {"tipo": "entrez", "email": "test@example.com"}
        with pytest.raises(ConfigError, match="Dataset tipo 'entrez' deve ter campos"):
            self.config_loader._validate_dataset(dataset)

    def test_validate_numeric_fields(self):
        """Test validating numeric fields."""
        config = {
            "nome": "teste",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "file", "filepath": "test.fasta"},
            "execucoes_por_algoritmo_por_base": 5,
            "num_bases": 3,
            "timeout": 600,
            "max_workers": 8,
            "seed": 42,
        }
        # Should not raise any exception
        self.config_loader._validate_numeric_fields(config)

    def test_validate_numeric_fields_invalid_type(self):
        """Test validating invalid numeric fields."""
        config = {
            "nome": "teste",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "file", "filepath": "test.fasta"},
            "execucoes_por_algoritmo_por_base": "not_a_number",
            "num_bases": 1,
        }
        with pytest.raises(ConfigError, match="deve ser do tipo int"):
            self.config_loader._validate_numeric_fields(config)

    def test_validate_numeric_fields_invalid_value(self):
        """Test validating invalid numeric values."""
        config = {
            "nome": "teste",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "file", "filepath": "test.fasta"},
            "execucoes_por_algoritmo_por_base": 0,  # Must be > 0
            "num_bases": 1,
        }
        with pytest.raises(ConfigError, match="tem valor inválido"):
            self.config_loader._validate_numeric_fields(config)

    def test_create_default_config(self):
        """Test creating default config."""
        default_config = self.config_loader.create_default_config()
        assert default_config is not None
        assert default_config["nome"] == "config_padrao"
        assert default_config["algoritmos"] == ["Baseline", "BLF-GA"]
        assert default_config["dataset"]["tipo"] == "file"
        assert default_config["dataset"]["filepath"] == "saved_datasets/dataset_custom.fasta"
        assert default_config["execucoes_por_algoritmo_por_base"] == 3
        assert default_config["num_bases"] == 1
        assert default_config["timeout"] == 300
        assert default_config["max_workers"] == 4
        assert default_config["seed"] == 42
        assert default_config["export_format"] == "csv"
        assert default_config["output_dir"] == "results"
        assert default_config["verbose"] is False

    def test_save_config(self):
        """Test saving config to file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            temp_path = f.name

        try:
            # Remove the file so save_config can create it
            os.unlink(temp_path)

            self.config_loader.save_config(self.valid_config, temp_path)

            # Verify file was created and contains correct data
            assert os.path.exists(temp_path)
            with open(temp_path) as f:
                saved_config = yaml.safe_load(f)

            assert saved_config["nome"] == "teste"
            assert saved_config["algoritmos"] == ["Baseline", "BLF-GA"]
            assert saved_config["dataset"]["tipo"] == "file"
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_save_config_permission_error(self):
        """Test saving config with permission error."""
        with pytest.raises(ConfigError, match="Erro ao salvar configuração"):
            self.config_loader.save_config(self.valid_config, "/root/readonly.yaml")

    def test_resolve_paths(self):
        """Test resolving relative paths in config."""
        config = {
            "nome": "teste",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "file", "filepath": "test.fasta"},
            "execucoes_por_algoritmo_por_base": 1,
            "num_bases": 1,
            "output_dir": "results",
        }

        config_dir = Path("/config/dir")
        resolved_config = self.config_loader._resolve_paths(config, config_dir)

        # Path should be resolved to absolute path
        assert resolved_config["output_dir"] == "/config/dir/results"
        assert resolved_config["dataset"]["filepath"] == "/config/dir/test.fasta"

    def test_resolve_paths_no_dataset(self):
        """Test resolving paths when no dataset section exists."""
        config = {"nome": "teste", "algoritmos": ["Baseline"], "execucoes_por_algoritmo_por_base": 1, "num_bases": 1}
        config_dir = Path("/config/dir")
        resolved = self.config_loader._resolve_paths(config, config_dir)
        assert resolved == config

    def test_resolve_paths_absolute_paths(self):
        """Test resolving paths when paths are already absolute."""
        config = {
            "nome": "teste",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "file", "filepath": "/absolute/path/test.fasta"},
            "execucoes_por_algoritmo_por_base": 1,
            "num_bases": 1,
            "output_dir": "/absolute/results",
        }
        config_dir = Path("/config/dir")
        resolved = self.config_loader._resolve_paths(config, config_dir)
        # Absolute paths should not change
        assert resolved["output_dir"] == "/absolute/results"
        assert resolved["dataset"]["filepath"] == "/absolute/path/test.fasta"
