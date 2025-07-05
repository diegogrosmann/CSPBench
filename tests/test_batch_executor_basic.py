"""
Testes básicos para BatchExecutor.
"""

import os
import tempfile

import pytest
import yaml

from src.core.exec.batch_executor import BatchConfig, BatchExecutor


class TestBatchConfig:
    """Testes para a classe BatchConfig."""

    def test_batch_config_creation(self):
        """Testa criação básica de BatchConfig."""
        config = {
            "nome": "teste",
            "algoritmos": ["Baseline", "BLF-GA"],
            "dataset": {"tipo": "synthetic", "n": 10, "L": 4},
            "execucoes_por_algoritmo_por_base": 5,
            "num_bases": 2,
        }

        batch_config = BatchConfig(config)

        assert batch_config.nome == "teste"
        assert batch_config.algoritmos == ["Baseline", "BLF-GA"]
        assert batch_config.execucoes_por_algoritmo_por_base == 5
        assert batch_config.num_bases == 2

    def test_batch_config_defaults(self):
        """Testa valores padrão de BatchConfig."""
        config = {"nome": "teste", "algoritmos": ["Baseline"], "dataset": {"tipo": "synthetic", "n": 10, "L": 4}}

        batch_config = BatchConfig(config)

        assert batch_config.execucoes_por_algoritmo_por_base == 3
        assert batch_config.num_bases == 1
        assert batch_config.timeout == 300

    def test_batch_config_str_representation(self):
        """Testa representação string de BatchConfig."""
        config = {"nome": "teste", "algoritmos": ["Baseline"], "dataset": {"tipo": "synthetic", "n": 10, "L": 4}}

        batch_config = BatchConfig(config)
        str_repr = str(batch_config)

        assert "teste" in str_repr
        assert "1 algoritmos" in str_repr
        assert "3 exec por base" in str_repr
        assert "1 bases" in str_repr


class TestBatchExecutor:
    """Testes para a classe BatchExecutor."""

    def setup_method(self):
        """Configuração inicial para cada teste."""
        self.config = {
            "nome": "teste",
            "algoritmos": ["Baseline"],
            "dataset": {"tipo": "synthetic", "n": 10, "L": 4},
            "execucoes_por_algoritmo_por_base": 1,
            "num_bases": 1,
        }

    def test_batch_executor_creation(self):
        """Testa criação de BatchExecutor."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(self.config, f)
            temp_path = f.name

        try:
            executor = BatchExecutor(temp_path)
            assert executor.config_file.name.endswith(".yaml")
        finally:
            os.unlink(temp_path)

    def test_batch_executor_invalid_file(self):
        """Testa erro com arquivo inexistente."""
        with pytest.raises(FileNotFoundError):
            BatchExecutor("/nonexistent/file.yaml")

    def test_batch_executor_with_valid_config_file(self):
        """Testa BatchExecutor com arquivo de configuração válido."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(self.config, f)
            temp_path = f.name

        try:
            executor = BatchExecutor(temp_path)
            assert executor.config_file.exists()
            assert str(executor.config_file).endswith(".yaml")
        finally:
            os.unlink(temp_path)
