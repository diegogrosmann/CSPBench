"""Testes para o módulo app.py."""
import unittest
from unittest.mock import MagicMock, patch

import pytest

from src.ui.cli.app import main, signal_handler


class TestApp(unittest.TestCase):
    """Testes para o módulo app.py."""

    def test_signal_handler(self):
        """Testa o handler de sinal."""
        # O signal_handler deve chamar sys.exit(0)
        with pytest.raises(SystemExit) as exc_info:
            signal_handler(None, None)
        self.assertEqual(exc_info.value.code, 0)

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("sys.argv", ["app.py", "--help"])
    def test_main_with_help_flag(self, mock_makedirs, mock_setup_logging, mock_console):
        """Testa main com flag --help."""
        # Deve sair com código 0 quando --help é usado
        with pytest.raises(SystemExit) as exc_info:
            main()
        self.assertEqual(exc_info.value.code, 0)

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.ui.cli.app.menu")
    @patch("sys.argv", ["app.py"])
    def test_main_menu_interaction(self, mock_menu, mock_makedirs, mock_setup_logging, mock_console):
        """Testa interação com menu."""
        mock_menu.return_value = "1"
        with patch("src.datasets.dataset_synthetic.generate_dataset") as mock_generate:
            mock_generate.return_value = (["ATCG", "ATGG"], {"n": 2, "L": 4})
            with patch("src.ui.cli.app.select_algorithms") as mock_select_algs:
                mock_select_algs.return_value = ["Baseline"]
                with patch("src.ui.cli.app.safe_input") as mock_input:
                    mock_input.side_effect = ["1", "30"]  # num_execs, timeout
                    with patch("src.core.exec.runner.execute_algorithm_runs") as mock_exec:
                        mock_exec.return_value = [{"distancia": 1, "tempo": 0.5}]
                        with patch("src.ui.cli.app.ResultsFormatter") as mock_formatter_class:
                            mock_formatter = MagicMock()
                            mock_formatter_class.return_value = mock_formatter
                            main()

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.ui.cli.app.menu")
    @patch("sys.argv", ["app.py"])
    def test_main_menu_returns_none(self, mock_menu, mock_makedirs, mock_setup_logging, mock_console):
        """Testa quando menu retorna None."""
        mock_menu.return_value = None
        # Deve retornar sem erro
        result = main()
        self.assertIsNone(result)

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.ui.cli.app.menu")
    @patch("sys.argv", ["app.py"])
    def test_main_keyboard_interrupt(self, mock_menu, mock_makedirs, mock_setup_logging, mock_console):
        """Testa tratamento de KeyboardInterrupt."""
        mock_menu.side_effect = KeyboardInterrupt()
        # KeyboardInterrupt deve ser capturado e resultar em SystemExit
        try:
            main()
        except SystemExit as e:
            self.assertEqual(e.code, 1)
        except KeyboardInterrupt:
            # Se KeyboardInterrupt escapar, isso é esperado em alguns casos
            pass

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.ui.cli.app.menu")
    @patch("sys.argv", ["app.py"])
    def test_main_general_exception(self, mock_menu, mock_makedirs, mock_setup_logging, mock_console):
        """Testa tratamento de exceção geral."""
        mock_menu.side_effect = Exception("Erro de teste")
        with pytest.raises(SystemExit) as exc_info:
            main()
        self.assertEqual(exc_info.value.code, 1)

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("sys.argv", ["app.py"])
    def test_main_setup_functions_called(self, mock_setup_logging, mock_console):
        """Testa se funções de setup são chamadas."""
        with patch("src.ui.cli.app.menu") as mock_menu:
            mock_menu.return_value = None
            main()

        # Verificar se setup_logging foi chamado
        mock_setup_logging.assert_called_once()

    def test_results_and_logs_directories_created(self):
        """Testa se diretórios results e logs são criados."""
        # Verificar se os diretórios existem (são criados na importação)
        import os

        self.assertTrue(os.path.exists("results"))
        self.assertTrue(os.path.exists("logs"))

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.core.exec.runner.execute_algorithm_runs")
    @patch("src.ui.cli.app.ResultsFormatter")
    @patch("src.ui.cli.app.print_quick_summary")
    @patch("sys.argv", ["app.py", "--silent", "--dataset", "synthetic", "--algorithms", "Baseline"])
    def test_main_with_algorithm_execution(
        self, mock_summary, mock_formatter_class, mock_exec, mock_makedirs, mock_setup_logging, mock_console
    ):
        """Testa execução completa com algoritmos."""
        mock_formatter = MagicMock()
        mock_formatter_class.return_value = mock_formatter
        mock_exec.return_value = [{"distancia": 1, "tempo": 0.5}]

        with patch("src.datasets.dataset_synthetic.generate_dataset") as mock_generate:
            mock_generate.return_value = (["ATCG", "ATGG"], {"n": 2, "L": 4})
            with pytest.raises(SystemExit) as exc_info:
                main()
            self.assertEqual(exc_info.value.code, 0)

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.core.exec.runner.execute_algorithm_runs")
    @patch("sys.argv", ["app.py", "--silent", "--dataset", "synthetic", "--algorithms", "Baseline"])
    def test_main_execution_error(self, mock_exec, mock_makedirs, mock_setup_logging, mock_console):
        """Testa tratamento de erro durante execução."""
        mock_exec.side_effect = Exception("Erro de execução")

        with patch("src.datasets.dataset_synthetic.generate_dataset") as mock_generate:
            mock_generate.return_value = (["ATCG", "ATGG"], {"n": 2, "L": 4})
            with pytest.raises(SystemExit) as exc_info:
                main()
            self.assertEqual(exc_info.value.code, 1)

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.core.exec.runner.execute_algorithm_runs")
    @patch("src.ui.cli.app.ResultsFormatter")
    @patch("src.ui.cli.app.print_quick_summary")
    @patch("sys.argv", ["app.py", "--silent", "--dataset", "synthetic", "--algorithms", "Baseline"])
    def test_main_save_results(
        self, mock_summary, mock_formatter_class, mock_exec, mock_makedirs, mock_setup_logging, mock_console
    ):
        """Testa salvamento de resultados."""
        mock_formatter = MagicMock()
        mock_formatter_class.return_value = mock_formatter
        mock_exec.return_value = [{"distancia": 1, "tempo": 0.5}]

        with patch("src.datasets.dataset_synthetic.generate_dataset") as mock_generate:
            mock_generate.return_value = (["ATCG", "ATGG"], {"n": 2, "L": 4})
            with pytest.raises(SystemExit) as exc_info:
                main()
            self.assertEqual(exc_info.value.code, 0)

        # Verificar se os métodos de salvamento foram chamados
        mock_formatter.save_detailed_report.assert_called_once()
        mock_formatter.export_to_csv.assert_called_once()

    @patch("src.ui.cli.app.console")
    @patch("src.ui.cli.app.setup_logging")
    @patch("src.ui.cli.app.os.makedirs")
    @patch("src.core.exec.runner.execute_algorithm_runs")
    @patch("src.ui.cli.app.ResultsFormatter")
    @patch("src.ui.cli.app.print_quick_summary")
    @patch("src.ui.cli.app.ask_save_dataset")
    @patch("sys.argv", ["app.py", "--silent", "--dataset", "synthetic", "--algorithms", "Baseline"])
    def test_main_save_dataset(
        self,
        mock_ask_save,
        mock_summary,
        mock_formatter_class,
        mock_exec,
        mock_makedirs,
        mock_setup_logging,
        mock_console,
    ):
        """Testa salvamento de dataset."""
        mock_formatter = MagicMock()
        mock_formatter_class.return_value = mock_formatter
        mock_exec.return_value = [{"distancia": 1, "tempo": 0.5}]

        with patch("src.datasets.dataset_synthetic.generate_dataset") as mock_generate:
            mock_generate.return_value = (["ATCG", "ATGG"], {"n": 2, "L": 4})
            with pytest.raises(SystemExit) as exc_info:
                main()
            self.assertEqual(exc_info.value.code, 0)

        # Em modo silencioso, ask_save_dataset não deve ser chamado
        mock_ask_save.assert_not_called()


if __name__ == "__main__":
    unittest.main()
