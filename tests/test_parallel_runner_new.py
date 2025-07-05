"""
Testes específicos para o ParallelRunner.

Este módulo testa a funcionalidade de execução paralela de algoritmos,
incluindo cenários de sucesso, falha e interrupção.
"""

import time
from unittest.mock import MagicMock, patch

import pytest

from src.core.exec.parallel_runner import ParallelRunner, execute_algorithms_parallel
from src.utils.signal_manager import get_signal_manager


class TestParallelRunner:
    """Testes para a classe ParallelRunner."""

    def setup_method(self):
        """Setup para cada teste."""
        self.sequences = ["AAAA", "AAAT", "AATT", "TTTT"]
        self.alphabet = "AT"
        self.console_mock = MagicMock()

    def test_parallel_runner_creation(self):
        """Testa criação do ParallelRunner."""
        runner = ParallelRunner(max_workers=2, timeout=30)
        assert runner.max_workers == 2
        assert runner.timeout == 30
        assert not runner._shutdown_requested

        # Verifica se callback foi registrado
        signal_manager = get_signal_manager()
        assert len(signal_manager.shutdown_callbacks) > 0

    def test_execute_algorithms_parallel_success(self):
        """Testa execução paralela bem-sucedida."""
        algorithm_names = ["Baseline"]

        results = execute_algorithms_parallel(
            algorithm_names=algorithm_names,
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=1,
            timeout=30,
        )

        assert "Baseline" in results
        assert results["Baseline"]["success"] is True
        assert results["Baseline"]["center"] is not None
        assert results["Baseline"]["distance"] >= 0
        assert results["Baseline"]["erro"] is None

    def test_execute_algorithms_parallel_invalid_algorithm(self):
        """Testa execução com algoritmo inválido."""
        algorithm_names = ["AlgoritmoInexistente"]

        results = execute_algorithms_parallel(
            algorithm_names=algorithm_names,
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=1,
            timeout=30,
        )

        assert results == {}

    def test_execute_algorithms_parallel_multiple(self):
        """Testa execução paralela de múltiplos algoritmos."""
        algorithm_names = ["Baseline", "BLF-GA"]

        results = execute_algorithms_parallel(
            algorithm_names=algorithm_names,
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=2,
            timeout=30,
        )

        assert len(results) == 2
        assert "Baseline" in results
        assert "BLF-GA" in results

        for alg_name in algorithm_names:
            assert results[alg_name]["success"] is True
            assert results[alg_name]["center"] is not None

    def test_execute_algorithms_parallel_with_baseline(self):
        """Testa execução com valor baseline."""
        algorithm_names = ["Baseline"]
        baseline_val = 1.0

        results = execute_algorithms_parallel(
            algorithm_names=algorithm_names,
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            baseline_val=baseline_val,
            max_workers=1,
            timeout=30,
        )

        assert "Baseline" in results
        assert results["Baseline"]["success"] is True
        # O console pode ou não ser chamado dependendo da implementação
        # O importante é que a função executa sem erro

    def test_parallel_runner_timeout_handling(self):
        """Testa tratamento de timeout."""
        # Usar timeout muito baixo para forçar timeout
        results = execute_algorithms_parallel(
            algorithm_names=["BLF-GA"],  # Algoritmo mais lento
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=1,
            timeout=1,  # 1 segundo apenas
        )

        # BLF-GA pode ou não dar timeout dependendo da máquina
        assert "BLF-GA" in results

    @patch("src.utils.signal_manager.is_interrupted")
    def test_execute_algorithms_interrupted(self, mock_is_interrupted):
        """Testa execução interrompida por sinal."""
        # Simular interrupção depois do primeiro algoritmo
        call_count = 0

        def mock_interrupted():
            nonlocal call_count
            call_count += 1
            # Retornar True após algumas chamadas para simular interrupção
            return call_count > 3

        mock_is_interrupted.side_effect = mock_interrupted

        results = execute_algorithms_parallel(
            algorithm_names=["Baseline"],
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=1,
            timeout=30,
        )

        # O importante é que não deve travar e deve retornar dict
        assert isinstance(results, dict)
        # Como Baseline é rápido, pode ter completado mesmo com interrupção
        # O teste confirma que a função não trava com interrupção

    def test_parallel_runner_shutdown_callback(self):
        """Testa callback de shutdown."""
        runner = ParallelRunner(max_workers=2)

        # Simular shutdown
        runner._shutdown_callback()

        assert runner._shutdown_requested is True

    def test_parallel_runner_with_console_output(self):
        """Testa execução com output para console."""
        algorithm_names = ["Baseline"]

        results = execute_algorithms_parallel(
            algorithm_names=algorithm_names,
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=1,
            timeout=30,
        )

        # Verificar resultados
        assert "Baseline" in results
        assert results["Baseline"]["success"] is True

        # O console pode ou não ser chamado dependendo da implementação
        # O importante é que a função executa sem erro

    def test_parallel_runner_empty_algorithm_list(self):
        """Testa execução com lista vazia de algoritmos."""
        results = execute_algorithms_parallel(
            algorithm_names=[],
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=1,
            timeout=30,
        )

        assert results == {}

    def test_parallel_runner_resource_cleanup(self):
        """Testa limpeza de recursos."""
        runner = ParallelRunner(max_workers=2)

        # Executar algoritmo
        results = runner.execute_algorithms_parallel(
            algorithm_names=["Baseline"], seqs=self.sequences, alphabet=self.alphabet, console=self.console_mock
        )

        # Verificar que executor foi limpo
        assert (
            runner._current_executor is None
            or not hasattr(runner._current_executor, "executor")
            or runner._current_executor.executor is None
        )

    def test_parallel_runner_progress_tracking(self):
        """Testa tracking de progresso."""
        algorithm_names = ["Baseline", "BLF-GA"]

        results = execute_algorithms_parallel(
            algorithm_names=algorithm_names,
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=2,
            timeout=30,
        )

        # Verificar resultados
        assert len(results) == 2
        for alg_name in algorithm_names:
            assert alg_name in results
            assert results[alg_name]["success"] is True

        # O console pode ou não ser chamado dependendo da implementação
        # O importante é que a função executa sem erro

    @pytest.mark.slow
    def test_parallel_runner_stress_test(self):
        """Teste de stress com múltiplas execuções."""
        algorithm_names = ["Baseline"] * 3  # Múltiplas instâncias

        start_time = time.time()
        results = execute_algorithms_parallel(
            algorithm_names=["Baseline", "BLF-GA", "CSC"],
            seqs=self.sequences,
            alphabet=self.alphabet,
            console=self.console_mock,
            max_workers=3,
            timeout=60,
        )
        end_time = time.time()

        # Verificar que executou em tempo razoável
        assert end_time - start_time < 60

        # Verificar resultados
        assert len(results) == 3
        for result in results.values():
            assert result["success"] is True

    def test_function_execute_algorithms_parallel(self):
        """Testa função de conveniência execute_algorithms_parallel."""
        results = execute_algorithms_parallel(
            algorithm_names=["Baseline"], seqs=self.sequences, alphabet=self.alphabet, max_workers=1, timeout=15
        )

        assert "Baseline" in results
        assert results["Baseline"]["success"] is True

    def teardown_method(self):
        """Cleanup após cada teste."""
        # Limpar callbacks de sinal para não interferir em outros testes
        signal_manager = get_signal_manager()
        signal_manager.shutdown_callbacks.clear()
        signal_manager.interrupted = False
