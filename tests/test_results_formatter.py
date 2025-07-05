"""Testes para o módulo results_formatter.py."""

import os
import tempfile
import unittest
from unittest.mock import patch

from src.core.io.results_formatter import ResultsFormatter


class TestResultsFormatter(unittest.TestCase):
    """Testes para a classe ResultsFormatter."""

    def setUp(self):
        """Configuração antes de cada teste."""
        self.formatter = ResultsFormatter()

    def test_results_formatter_creation(self):
        """Testa criação do ResultsFormatter."""
        formatter = ResultsFormatter()

        assert formatter.results == {}
        assert formatter.extra_info == {}
        assert formatter.exporter is not None

    def test_add_algorithm_results(self):
        """Testa adição de resultados de algoritmo."""
        executions = [
            {"distancia": 2, "tempo": 1.5, "melhor_string": "ACGT"},
            {"distancia": 3, "tempo": 1.2, "melhor_string": "ATGT"},
        ]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)

        assert "TestAlgorithm" in self.formatter.results
        assert self.formatter.results["TestAlgorithm"] == executions

    def test_add_multiple_algorithm_results(self):
        """Testa adição de resultados de múltiplos algoritmos."""
        executions1 = [{"distancia": 2, "tempo": 1.5}]
        executions2 = [{"distancia": 3, "tempo": 2.0}]

        self.formatter.add_algorithm_results("Algorithm1", executions1)
        self.formatter.add_algorithm_results("Algorithm2", executions2)

        assert len(self.formatter.results) == 2
        assert "Algorithm1" in self.formatter.results
        assert "Algorithm2" in self.formatter.results

    def test_format_detailed_results_empty(self):
        """Testa formatação de resultados detalhados vazios."""
        result = self.formatter.format_detailed_results()

        assert isinstance(result, str)
        assert "Nenhum resultado disponível" in result

    def test_format_detailed_results_with_data(self):
        """Testa formatação de resultados detalhados com dados."""
        executions = [
            {"distancia": 2, "tempo": 1.5, "melhor_string": "ACGT", "iteracoes": 100},
            {"distancia": 3, "tempo": 1.2, "melhor_string": "ATGT", "iteracoes": 150},
        ]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        result = self.formatter.format_detailed_results()

        assert isinstance(result, str)
        assert "TestAlgorithm" in result
        assert "ACGT" in result
        assert "ATGT" in result

    def test_format_quick_summary_empty(self):
        """Testa formatação de resumo rápido vazio."""
        result = self.formatter.format_quick_summary()

        assert isinstance(result, str)
        assert "Nenhum resultado disponível" in result

    def test_format_quick_summary_with_data(self):
        """Testa formatação de resumo rápido com dados."""
        executions = [{"distancia": 2, "tempo": 1.5}, {"distancia": 3, "tempo": 1.2}]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        result = self.formatter.format_quick_summary()

        assert isinstance(result, str)
        assert "TestAlgorithm" in result

    def test_save_detailed_report(self):
        """Testa salvamento de relatório detalhado."""
        executions = [
            {"distancia": 2, "tempo": 1.5, "melhor_string": "ACGT"},
        ]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            try:
                self.formatter.save_detailed_report(f.name)

                # Verificar se arquivo foi criado
                assert os.path.exists(f.name)

                # Verificar conteúdo do arquivo
                with open(f.name) as file:
                    content = file.read()
                    assert "TestAlgorithm" in content
                    assert "ACGT" in content

            finally:
                # Limpar arquivo temporário
                os.unlink(f.name)

    def test_save_detailed_report_with_extra_info(self):
        """Testa salvamento de relatório com informações extras."""
        executions = [
            {"distancia": 2, "tempo": 1.5, "melhor_string": "ACGT"},
        ]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        self.formatter.extra_info = {"seed": 42, "dataset_strings": ["ACGT", "ATGT"], "params": {"n": 2, "L": 4}}

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            try:
                self.formatter.save_detailed_report(f.name)

                # Verificar se arquivo foi criado
                assert os.path.exists(f.name)

                # Verificar conteúdo do arquivo
                with open(f.name) as file:
                    content = file.read()
                    # Verificar que informações extras estão presentes
                    assert "params" in content or "Parâmetros" in content

            finally:
                # Limpar arquivo temporário
                os.unlink(f.name)

    def test_export_to_csv(self):
        """Testa exportação para CSV."""
        executions = [
            {"distancia": 2, "tempo": 1.5, "melhor_string": "ACGT", "iteracoes": 100},
            {"distancia": 3, "tempo": 1.2, "melhor_string": "ATGT", "iteracoes": 150},
        ]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv") as f:
            try:
                with patch.object(self.formatter.exporter, "export_to_csv") as mock_export:
                    self.formatter.export_to_csv(f.name)

                    # Verificar se o método de exportação foi chamado
                    mock_export.assert_called_once()

            finally:
                # Limpar arquivo temporário
                if os.path.exists(f.name):
                    os.unlink(f.name)

    def test_export_to_csv_empty_results(self):
        """Testa exportação para CSV com resultados vazios."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv") as f:
            try:
                with patch.object(self.formatter.exporter, "export_to_csv") as mock_export:
                    self.formatter.export_to_csv(f.name)

                    # Deve tentar exportar mesmo com resultados vazios
                    mock_export.assert_called_once()

            finally:
                # Limpar arquivo temporário
                if os.path.exists(f.name):
                    os.unlink(f.name)

    def test_format_execution_summary_single_execution(self):
        """Testa formatação de resumo de execuções simples."""
        executions = [{"distancia": 2, "tempo": 1.5}]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        result = self.formatter.format_quick_summary()

        assert "TestAlgorithm" in result
        assert "2" in result  # distancia
        assert "1.5" in result  # tempo

    def test_format_execution_summary_multiple_executions(self):
        """Testa formatação de resumo de múltiplas execuções."""
        executions = [{"distancia": 2, "tempo": 1.5}, {"distancia": 3, "tempo": 1.2}, {"distancia": 1, "tempo": 2.0}]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        result = self.formatter.format_quick_summary()

        assert "TestAlgorithm" in result
        # Deve mostrar estatísticas (melhor, média, etc.)

    def test_format_execution_summary_with_errors(self):
        """Testa formatação de resumo com execuções com erro."""
        executions = [
            {"distancia": 2, "tempo": 1.5},
            {"erro": "Timeout", "tempo": 30.0},
            {"distancia": 3, "tempo": 1.2},
        ]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        result = self.formatter.format_quick_summary()

        assert "TestAlgorithm" in result

    def test_get_algorithm_stats_valid_results(self):
        """Testa obtenção de estatísticas de algoritmo com resultados válidos."""
        executions = [{"distancia": 2, "tempo": 1.5}, {"distancia": 3, "tempo": 1.2}, {"distancia": 1, "tempo": 2.0}]

        # Testar se o formatter consegue processar estatísticas
        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        result = self.formatter.format_detailed_results()

        assert isinstance(result, str)
        assert "TestAlgorithm" in result

    def test_get_algorithm_stats_with_infinite_distance(self):
        """Testa obtenção de estatísticas com distâncias infinitas."""
        executions = [
            {"distancia": 2, "tempo": 1.5},
            {"distancia": float("inf"), "tempo": 30.0},
            {"distancia": 3, "tempo": 1.2},
        ]

        self.formatter.add_algorithm_results("TestAlgorithm", executions)
        result = self.formatter.format_detailed_results()

        assert isinstance(result, str)
        assert "TestAlgorithm" in result

    def test_results_formatter_with_complex_data(self):
        """Testa formatter com dados complexos."""
        executions = [
            {
                "distancia": 2,
                "tempo": 1.5,
                "melhor_string": "ACGT",
                "iteracoes": 100,
                "seed": 42,
                "parametros": {"mutation_rate": 0.1},
            },
            {"distancia": 3, "tempo": 1.2, "melhor_string": "ATGT", "iteracoes": 150, "seed": 42},
        ]

        self.formatter.add_algorithm_results("ComplexAlgorithm", executions)
        self.formatter.extra_info = {
            "dataset_strings": ["ACGT", "ATGT", "ACCT"],
            "params": {"n": 3, "L": 4},
            "seed": 42,
        }

        # Teste de formatação detalhada
        detailed = self.formatter.format_detailed_results()
        assert "ComplexAlgorithm" in detailed
        assert "ACGT" in detailed

        # Teste de resumo rápido
        summary = self.formatter.format_quick_summary()
        assert "ComplexAlgorithm" in summary

    def test_results_formatter_str_representation(self):
        """Testa representação string do formatter."""
        self.formatter.add_algorithm_results("TestAlgorithm", [{"distancia": 2}])

        str_repr = str(self.formatter)
        assert isinstance(str_repr, str)

    def test_results_formatter_with_no_valid_results(self):
        """Testa formatter com execuções sem resultados válidos."""
        executions = [{"erro": "Timeout", "tempo": 30.0}, {"erro": "Memory error", "tempo": 25.0}]

        self.formatter.add_algorithm_results("FailingAlgorithm", executions)

        detailed = self.formatter.format_detailed_results()
        summary = self.formatter.format_quick_summary()

        assert "FailingAlgorithm" in detailed
        assert "FailingAlgorithm" in summary


if __name__ == "__main__":
    unittest.main()
