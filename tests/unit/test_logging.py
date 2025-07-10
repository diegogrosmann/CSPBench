"""
CSPBench - Testes para utilitários de logging.

Este módulo contém testes para as funcionalidades de logging
incluindo configuração, níveis e gerenciamento de arquivos.
"""

import logging
import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from src.utils.logging import (
    cleanup_old_logs,
    configure_algorithm_logging,
    get_logger,
    log_execution_context,
    set_module_level,
    setup_logging,
)


class TestSetupLogging(unittest.TestCase):
    """Testes para setup_logging."""

    def test_setup_logging_default(self):
        """Testa configuração padrão do logging."""
        with tempfile.TemporaryDirectory() as temp_dir:
            logger = setup_logging(base_name="test", log_dir=temp_dir)
            assert logger is not None
            assert isinstance(logger, logging.Logger)
            # A função setup_logging retorna o logger root configurado

    def test_setup_logging_with_debug(self):
        """Testa configuração com debug habilitado."""
        with tempfile.TemporaryDirectory() as temp_dir:
            logger = setup_logging(base_name="test", log_dir=temp_dir, debug=True)
            assert logger is not None
            assert isinstance(logger, logging.Logger)

    def test_setup_logging_with_console(self):
        """Testa configuração com console."""
        with tempfile.TemporaryDirectory() as temp_dir:
            logger = setup_logging(base_name="test", log_dir=temp_dir, console=True)
            assert logger is not None
            assert isinstance(logger, logging.Logger)


class TestGetLogger(unittest.TestCase):
    """Testes para get_logger."""

    def test_get_logger_default(self):
        """Testa obtenção de logger padrão."""
        logger = get_logger("test_module")
        assert logger is not None
        assert isinstance(logger, logging.Logger)
        assert logger.name == "test_module"

    def test_get_logger_with_level(self):
        """Testa obtenção de logger com nível específico."""
        logger = get_logger("test_module", level="DEBUG")
        assert logger is not None
        assert isinstance(logger, logging.Logger)

    def test_get_logger_multiple_calls(self):
        """Testa múltiplas chamadas para o mesmo logger."""
        logger1 = get_logger("test_module")
        logger2 = get_logger("test_module")
        assert logger1 is logger2


class TestSetModuleLevel(unittest.TestCase):
    """Testes para set_module_level."""

    def test_set_module_level_debug(self):
        """Testa configuração de nível DEBUG."""
        set_module_level("test_module", "DEBUG")
        # Verifica se não houve exceção

    def test_set_module_level_info(self):
        """Testa configuração de nível INFO."""
        set_module_level("test_module", "INFO")
        # Verifica se não houve exceção

    def test_set_module_level_invalid(self):
        """Testa configuração com nível inválido."""
        with self.assertRaises(ValueError):
            set_module_level("test_module", "INVALID")


class TestConfigureAlgorithmLogging(unittest.TestCase):
    """Testes para configure_algorithm_logging."""

    def test_configure_algorithm_logging_default(self):
        """Testa configuração padrão para algoritmo."""
        logger = configure_algorithm_logging("BLF-GA")
        assert logger is not None
        assert isinstance(logger, logging.Logger)

    def test_configure_algorithm_logging_with_debug(self):
        """Testa configuração com debug para algoritmo."""
        logger = configure_algorithm_logging("BLF-GA", verbose=True)
        assert logger is not None
        assert isinstance(logger, logging.Logger)


class TestLogExecutionContext(unittest.TestCase):
    """Testes para log_execution_context."""

    def test_log_execution_context_basic(self):
        """Testa logging de contexto básico."""
        context = {"algorithm": "BLF-GA", "dataset": "test"}
        logger = logging.getLogger("test")
        # Não deve gerar exceção
        log_execution_context(logger, context)

    def test_log_execution_context_empty(self):
        """Testa logging de contexto vazio."""
        context = {}
        logger = logging.getLogger("test")
        # Não deve gerar exceção
        log_execution_context(logger, context)

    def test_log_execution_context_complex(self):
        """Testa logging de contexto complexo."""
        context = {
            "algorithm": "BLF-GA",
            "dataset": {"type": "synthetic", "n": 10},
            "config": {"param1": "value1", "param2": 42},
        }
        logger = logging.getLogger("test")
        # Não deve gerar exceção
        log_execution_context(logger, context)


class TestCleanupOldLogs(unittest.TestCase):
    """Testes para cleanup_old_logs."""

    def test_cleanup_old_logs_empty_dir(self):
        """Testa limpeza em diretório vazio."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Não deve gerar exceção
            cleanup_old_logs(temp_dir)

    def test_cleanup_old_logs_with_files(self):
        """Testa limpeza com arquivos existentes."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Cria alguns arquivos de log
            for i in range(5):
                log_file = os.path.join(temp_dir, f"test_{i}.log")
                with open(log_file, "w") as f:
                    f.write("test log content")

            # Não deve gerar exceção
            cleanup_old_logs(temp_dir)

    def test_cleanup_old_logs_nonexistent_dir(self):
        """Testa limpeza em diretório inexistente."""
        # Não deve gerar exceção
        cleanup_old_logs("/path/that/does/not/exist")


if __name__ == "__main__":
    unittest.main()
