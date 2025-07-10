# Testes Unitários para Módulo de Configuração - CSPBench
# Testa as funcionalidades do módulo src.utils.config

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from src.utils.config import (
    get_config_value,
    merge_configs,
    safe_input,
    validate_config,
)


class TestValidateConfig:
    """Testes para função validate_config."""

    def test_validate_valid_config(self):
        """Testa validação de configuração válida."""
        config = {
            "algorithm": "BLF-GA",
            "dataset": {"type": "synthetic"},
            "optimization": {"n_trials": 100},
        }

        # Deve retornar True para configuração válida
        assert validate_config(config) is True

    def test_validate_empty_config(self):
        """Testa validação de configuração vazia."""
        config = {}

        # A função validate_config pode aceitar configuração vazia como válida
        result = validate_config(config)
        assert isinstance(result, bool)

    def test_validate_invalid_config(self):
        """Testa validação de configuração inválida."""
        config = {"algorithm": None, "dataset": {}}

        # Deve retornar False para configuração inválida quando chaves obrigatórias faltam
        assert (
            validate_config(config, required_keys=["algorithm", "optimization"])
            is False
        )


class TestMergeConfigs:
    """Testes para função merge_configs."""

    def test_merge_two_configs(self):
        """Testa merge de duas configurações."""
        config1 = {"algorithm": "BLF-GA", "dataset": {"type": "synthetic"}}
        config2 = {"optimization": {"n_trials": 100}}

        merged = merge_configs(config1, config2)

        assert merged["algorithm"] == "BLF-GA"
        assert merged["dataset"]["type"] == "synthetic"
        assert merged["optimization"]["n_trials"] == 100

    def test_merge_overlapping_configs(self):
        """Testa merge de configurações com sobreposição."""
        config1 = {"algorithm": "BLF-GA", "dataset": {"type": "synthetic"}}
        config2 = {"algorithm": "CSC", "dataset": {"n": 10}}

        merged = merge_configs(config1, config2)

        # Segundo config deve sobrescrever primeiro
        assert merged["algorithm"] == "CSC"
        assert merged["dataset"]["n"] == 10

    def test_merge_empty_configs(self):
        """Testa merge de configurações vazias."""
        config1 = {}
        config2 = {}

        merged = merge_configs(config1, config2)

        assert merged == {}


class TestGetConfigValue:
    """Testes para função get_config_value."""

    def test_get_config_value_with_default(self):
        """Testa obtenção de valor com default."""
        config = {"algorithm": "BLF-GA"}

        value = get_config_value("algorithm", "default", config)
        assert value == "BLF-GA"

        value = get_config_value("nonexistent", "default", config)
        assert value == "default"

    def test_get_config_value_nested(self):
        """Testa obtenção de valor com configuração aninhada."""
        config = {"dataset": {"type": "synthetic", "n": 10}}

        # Função não suporta chaves aninhadas, testamos apenas chaves de primeiro nível
        value = get_config_value("dataset", None, config)
        assert value == {"type": "synthetic", "n": 10}

        value = get_config_value("nonexistent", "default", config)
        assert value == "default"


class TestSafeInput:
    """Testes para função safe_input."""

    @patch("builtins.input")
    def test_safe_input_with_value(self, mock_input):
        """Testa entrada com valor fornecido."""
        mock_input.return_value = "test_value"

        result = safe_input("Enter value: ", "default")
        assert result == "test_value"

    @patch("builtins.input")
    def test_safe_input_with_default(self, mock_input):
        """Testa entrada com valor padrão."""
        mock_input.return_value = ""

        result = safe_input("Enter value: ", "default")
        assert result == "default"

    @patch("builtins.input")
    def test_safe_input_keyboard_interrupt(self, mock_input):
        """Testa entrada com interrupção do teclado."""
        mock_input.side_effect = KeyboardInterrupt()

        with pytest.raises(SystemExit):
            safe_input("Enter value: ", "default")

    @patch("builtins.input")
    def test_safe_input_eof_error(self, mock_input):
        """Testa entrada com erro EOF."""
        mock_input.side_effect = EOFError()

        with pytest.raises(SystemExit):
            safe_input("Enter value: ", "default")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
