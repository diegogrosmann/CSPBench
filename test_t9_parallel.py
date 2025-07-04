#!/usr/bin/env python3
"""
Exemplo de uso do sistema paralelo modernizado (T9).

Demonstra o uso de:
- CSPAlgorithm ABC
- ModernParallelExecutor
- ParallelRunner
"""

import logging
import time
from pprint import pprint

# Configurar logging
logging.basicConfig(level=logging.INFO)

from algorithms.base import CSPAlgorithm

# Imports usando a nova estrutura
from algorithms.baseline.algorithm import BaselineAlg
from algorithms.blf_ga.algorithm import BLFGAAlgorithm
from csp_blfga.core.exec.algorithm_executor import ModernParallelExecutor
from csp_blfga.core.exec.parallel_runner import ParallelRunner
from csp_blfga.ui.cli.console_manager import ConsoleManager


def test_modern_parallel_system():
    """Testa o sistema paralelo modernizado."""
    print("🚀 Testando Sistema Paralelo Modernizado (T9)")
    print("=" * 50)

    # Dados de teste
    seqs = ["ACGTACGTACGT", "ACGTACGTACGT", "ACGTACGTACGT", "ACGTACGTACGT"]
    alphabet = "ACGT"

    # Verificar que algoritmos são CSPAlgorithm
    print("📋 Verificando interfaces dos algoritmos:")
    print(f"  BaselineAlg é CSPAlgorithm: {issubclass(BaselineAlg, CSPAlgorithm)}")
    print(
        f"  BLFGAAlgorithm é CSPAlgorithm: {issubclass(BLFGAAlgorithm, CSPAlgorithm)}"
    )
    print()

    # Teste 1: Execução sequencial de algoritmo individual
    print("🔧 Teste 1: Execução sequencial")
    baseline = BaselineAlg(seqs, alphabet)

    start_time = time.time()
    result = baseline.run()
    exec_time = time.time() - start_time

    print(f"  Resultado: {len(result)} elementos")
    print(f"  Centro: {result[0]}")
    print(f"  Distância: {result[1]}")
    print(f"  Metadata: {result[2]}")
    print(f"  Tempo: {exec_time:.3f}s")
    print()

    # Teste 2: Execução paralela com ModernParallelExecutor
    print("⚡ Teste 2: Execução paralela com ModernParallelExecutor")

    executor = ModernParallelExecutor(max_workers=2, timeout=60)

    tasks = [
        {
            "alg_class": BaselineAlg,
            "strings": seqs,
            "alphabet": alphabet,
            "params": {},
            "name": "Baseline_1",
        },
        {
            "alg_class": BaselineAlg,
            "strings": seqs,
            "alphabet": alphabet,
            "params": {},
            "name": "Baseline_2",
        },
    ]

    start_time = time.time()
    results = executor.execute_algorithm_batch(tasks)
    exec_time = time.time() - start_time

    print(f"  Executadas {len(tasks)} tarefas em {exec_time:.3f}s")
    print("  Resultados:")
    for i, result in enumerate(results):
        print(
            f"    Tarefa {i+1}: sucesso={result['success']}, dist={result['distance']}"
        )
    print()

    # Teste 3: ParallelRunner com múltiplos algoritmos
    print("🎯 Teste 3: ParallelRunner com múltiplos algoritmos")

    console = ConsoleManager()
    runner = ParallelRunner(max_workers=2, timeout=60)

    algorithm_names = ["Baseline"]  # Apenas Baseline para teste rápido

    start_time = time.time()
    results = runner.execute_algorithms_parallel(
        algorithm_names=algorithm_names,
        seqs=seqs,
        alphabet=alphabet,
        console=console,
        baseline_val=None,
    )
    exec_time = time.time() - start_time

    print(f"  Executados {len(algorithm_names)} algoritmos em {exec_time:.3f}s")
    print("  Resultados finais:")
    pprint(results)
    print()

    print("✅ Todos os testes do sistema paralelo modernizado foram concluídos!")
    print("🎉 Tarefa T9 implementada com sucesso!")


if __name__ == "__main__":
    test_modern_parallel_system()
