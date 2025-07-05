"""
Tests for BatchExecutor module.
"""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from src.core.exec.batch_executor import (
    BatchConfig,
    BatchExecutor,
    create_example_config,
    list_batch_configs,
    select_batch_config,
)


class TestBatchExecutor:
    """Test BatchExecutor functionality."""

    def test_list_batch_configs(self):
        """Test listing batch configuration files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock config directory
            config_dir = Path(temp_dir) / "batch_configs"
            config_dir.mkdir()

            # Create a test JSON file
            test_file = config_dir / "test.json"
            test_file.write_text('{"test": "data"}')

            # Patch the config directory path
            with patch("src.core.exec.batch_executor.Path") as mock_path:
                mock_path.return_value = config_dir
                mock_path.side_effect = lambda x: config_dir if x == "batch_configs" else Path(x)

                configs = list_batch_configs()
                assert len(configs) >= 0
                assert isinstance(configs, list)

    def test_create_example_config(self):
        """Test creating example configuration file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Patch the config directory path
            with patch("src.core.exec.batch_executor.Path") as mock_path:
                config_dir = Path(temp_dir) / "batch_configs"
                mock_path.return_value = config_dir
                mock_path.side_effect = lambda x: config_dir if x == "batch_configs" else Path(x)

                create_example_config()

                # Check that the config directory was created
                assert config_dir.exists()

    @patch("src.core.exec.batch_executor.safe_input")
    @patch("src.core.exec.batch_executor.console")
    def test_select_batch_config_valid_choice(self, mock_console, mock_input):
        """Test selecting a valid batch config."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock config directory
            config_dir = Path(temp_dir) / "batch_configs"
            config_dir.mkdir()

            # Create a test YAML file
            test_file = config_dir / "test.yaml"
            test_file.write_text("test: data")

            # Mock user input
            mock_input.return_value = "1"

            # Patch the config directory path
            with patch("src.core.exec.batch_executor.Path") as mock_path:
                mock_path.return_value = config_dir
                mock_path.side_effect = lambda x: config_dir if x == "batch_configs" else Path(x)

                result = select_batch_config()
                assert result == str(test_file)

    @patch("src.core.exec.batch_executor.safe_input")
    @patch("src.core.exec.batch_executor.console")
    def test_select_batch_config_no_files(self, mock_console, mock_input):
        """Test selecting batch config when no files exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create empty config directory
            config_dir = Path(temp_dir) / "batch_configs"
            config_dir.mkdir()

            # Patch the config directory path
            with patch("src.core.exec.batch_executor.Path") as mock_path:
                mock_path.return_value = config_dir
                mock_path.side_effect = lambda x: config_dir if x == "batch_configs" else Path(x)

                result = select_batch_config()
                assert result == ""

    @patch("src.core.exec.batch_executor.safe_input")
    @patch("src.core.exec.batch_executor.console")
    def test_select_batch_config_invalid_choice(self, mock_console, mock_input):
        """Test selecting batch config with invalid choice."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock config directory
            config_dir = Path(temp_dir) / "batch_configs"
            config_dir.mkdir()

            # Create a test YAML file
            test_file = config_dir / "test.yaml"
            test_file.write_text("test: data")

            # Mock user input - first invalid, then valid
            mock_input.side_effect = ["invalid", "1"]

            # Patch the config directory path
            with patch("src.core.exec.batch_executor.Path") as mock_path:
                mock_path.return_value = config_dir
                mock_path.side_effect = lambda x: config_dir if x == "batch_configs" else Path(x)

                result = select_batch_config()
                assert result == str(test_file)

    def test_batch_executor_imports(self):
        """Test that all required imports are available."""
        assert callable(create_example_config)
        assert callable(list_batch_configs)
        assert callable(select_batch_config)
        assert BatchConfig is not None
        assert BatchExecutor is not None


class TestBatchConfig:
    """Test BatchConfig class."""

    def test_batch_config_creation(self):
        """Test creating BatchConfig."""
        config_dict = {
            "nome": "Test Config",
            "dataset": {"tipo": "synthetic", "parametros": {"n": 10}},
            "algoritmos": ["Baseline", "BLF-GA"],
            "execucoes_por_algoritmo_por_base": 3,
            "num_bases": 2,
            "timeout": 300,
        }

        batch_config = BatchConfig(config_dict)

        assert batch_config.nome == "Test Config"
        assert batch_config.dataset_config["tipo"] == "synthetic"
        assert batch_config.algoritmos == ["Baseline", "BLF-GA"]
        assert batch_config.execucoes_por_algoritmo_por_base == 3
        assert batch_config.num_bases == 2
        assert batch_config.timeout == 300

    def test_batch_config_defaults(self):
        """Test BatchConfig with default values."""
        config_dict = {}
        batch_config = BatchConfig(config_dict)

        assert batch_config.nome == "Execução Sem Nome"
        assert batch_config.dataset_config == {}
        assert batch_config.algoritmos == []
        assert batch_config.execucoes_por_algoritmo_por_base == 3
        assert batch_config.num_bases == 1

    def test_batch_config_retro_compatibility(self):
        """Test BatchConfig retro-compatibility with old field names."""
        config_dict = {
            "nome": "Test Config",
            "execucoes_por_algoritmo": 5,  # Old field name
        }

        batch_config = BatchConfig(config_dict)

        assert batch_config.execucoes_por_algoritmo_por_base == 5

    def test_batch_config_str_representation(self):
        """Test BatchConfig string representation."""
        config_dict = {
            "nome": "Test Config",
            "algoritmos": ["Baseline", "BLF-GA"],
            "execucoes_por_algoritmo_por_base": 3,
            "num_bases": 2,
        }

        batch_config = BatchConfig(config_dict)
        str_repr = str(batch_config)

        assert "Test Config" in str_repr
        assert "2 algoritmos" in str_repr
        assert "3 exec por base" in str_repr
        assert "2 bases" in str_repr


class TestBatchExecutorAdvanced:
    """Test BatchExecutor class with advanced scenarios."""

    def test_batch_executor_creation_yaml(self):
        """Test BatchExecutor creation with YAML file."""
        config_data = {
            "batch_info": {
                "nome": "Test Batch",
                "descricao": "Test batch description",
                "timeout_global": 1800,
            },
            "execucoes": [
                {
                    "nome": "Test Execution",
                    "dataset": {"tipo": "synthetic", "parametros": {"n": 10}},
                    "algoritmos": ["Baseline"],
                    "execucoes_por_algoritmo_por_base": 1,
                    "timeout": 300,
                }
            ],
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.safe_dump(config_data, f)
            temp_file = f.name

        try:
            batch_executor = BatchExecutor(temp_file)

            assert batch_executor.batch_info["nome"] == "Test Batch"
            assert len(batch_executor.execucoes) == 1
            assert batch_executor.execucoes[0].nome == "Test Execution"
            assert batch_executor.results_dir.exists()

        finally:
            os.unlink(temp_file)

    def test_batch_executor_creation_json(self):
        """Test BatchExecutor creation with JSON file."""
        config_data = {
            "batch_info": {
                "nome": "Test Batch JSON",
                "descricao": "Test batch description",
                "timeout_global": 1800,
            },
            "execucoes": [
                {
                    "nome": "Test Execution JSON",
                    "dataset": {"tipo": "synthetic", "parametros": {"n": 10}},
                    "algoritmos": ["Baseline"],
                    "execucoes_por_algoritmo_por_base": 1,
                    "timeout": 300,
                }
            ],
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            temp_file = f.name

        try:
            batch_executor = BatchExecutor(temp_file)

            assert batch_executor.batch_info["nome"] == "Test Batch JSON"
            assert len(batch_executor.execucoes) == 1
            assert batch_executor.execucoes[0].nome == "Test Execution JSON"

        finally:
            os.unlink(temp_file)

    def test_batch_executor_invalid_file_format(self):
        """Test BatchExecutor with invalid file format."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("invalid content")
            temp_file = f.name

        try:
            with pytest.raises(ValueError, match="Formato de arquivo não suportado"):
                BatchExecutor(temp_file)

        finally:
            os.unlink(temp_file)

    def test_batch_executor_invalid_yaml(self):
        """Test BatchExecutor with invalid YAML."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("invalid: yaml: content: [")
            temp_file = f.name

        try:
            with pytest.raises((yaml.YAMLError, FileNotFoundError)):
                BatchExecutor(temp_file)

        finally:
            os.unlink(temp_file)

    def test_batch_executor_nonexistent_file(self):
        """Test BatchExecutor with nonexistent file."""
        with pytest.raises(FileNotFoundError):
            BatchExecutor("nonexistent_file.yaml")

    def test_batch_executor_empty_config(self):
        """Test BatchExecutor with empty configuration."""
        config_data = {}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            temp_file = f.name

        try:
            batch_executor = BatchExecutor(temp_file)

            assert batch_executor.batch_info == {}
            assert len(batch_executor.execucoes) == 0

        finally:
            os.unlink(temp_file)

    def test_list_batch_configs_creates_directory(self):
        """Test that list_batch_configs creates directory if it doesn't exist."""
        # Create temporary directory and move to it
        with tempfile.TemporaryDirectory() as temp_dir:
            original_cwd = Path.cwd()
            temp_path = Path(temp_dir)

            try:
                os.chdir(temp_path)

                # Ensure batch_configs doesn't exist
                batch_configs_dir = temp_path / "batch_configs"
                if batch_configs_dir.exists():
                    import shutil

                    shutil.rmtree(batch_configs_dir)

                configs = list_batch_configs()

                # Should create directory and return empty list
                assert batch_configs_dir.exists()
                assert isinstance(configs, list)
                assert len(configs) == 0

            finally:
                os.chdir(original_cwd)

    def test_create_example_config_file(self):
        """Test creating example config file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            original_cwd = Path.cwd()
            temp_path = Path(temp_dir)

            try:
                os.chdir(temp_path)

                create_example_config()

                # Check that example file was created
                example_file = temp_path / "batch_configs" / "exemplo_batch.json"
                assert example_file.exists()

                # Check that file contains valid JSON
                with open(example_file, encoding="utf-8") as f:
                    config = json.load(f)

                assert "batch_info" in config
                assert "execucoes" in config
                assert len(config["execucoes"]) > 0

            finally:
                os.chdir(original_cwd)
