#!/usr/bin/env python3
"""
Exemplo de uso do sistema paralelo T9 - CSPAlgorithm + ParallelRunner.

Este script demonstra como usar:
- CSPAlgorithm (ABC moderna)
- ParallelRunner (execução paralela)
- ProcessPoolExecutor (paralelização avançada)
"""

import sys
from pathlib import Path

# Adicionar diretório raiz ao path
sys.path.insert(0, str(Path(__file__).parent.parent))

from csp_blfga.core.exec.parallel_runner import execute_algorithms_parallel
from csp_blfga.ui.cli.console_manager import ConsoleManager
from datasets.dataset_synthetic import generate_synthetic_dataset


def main():
    """Demonstra o uso do sistema paralelo T9."""

    console = ConsoleManager()

    console.print("🚀 Exemplo T9 - Sistema Paralelo CSPAlgorithm")
    console.print("=" * 60)

    # Carregar dataset de exemplo
    try:
        # Gerar dataset sintético simples
        seqs = generate_synthetic_dataset(n=3, L=10, alphabet="ACGT", noise=0.1)
        alphabet = "ACGT"

        console.print(
            f"📊 Dataset carregado: {len(seqs)} sequências de tamanho {len(seqs[0])}"
        )
        console.print(f"🔤 Alfabeto: {alphabet}")
        console.print()

    except Exception as e:
        console.print(f"❌ Erro ao carregar dataset: {e}")
        return

    # Lista de algoritmos a executar
    algorithms = ["Baseline", "BLF-GA", "DP-CSP", "CSC", "H³-CSP"]

    console.print(f"🔄 Executando {len(algorithms)} algoritmos em paralelo...")
    console.print(f"📋 Algoritmos: {', '.join(algorithms)}")
    console.print()

    # Executar algoritmos em paralelo
    results = execute_algorithms_parallel(
        algorithm_names=algorithms,
        seqs=seqs,
        alphabet=alphabet,
        console=console,
        baseline_val=None,
        max_workers=3,  # Usar 3 workers
        timeout=60,  # Timeout de 60 segundos
    )

    console.print("\n" + "=" * 60)
    console.print("📊 RESULTADOS DA EXECUÇÃO PARALELA")
    console.print("=" * 60)

    # Mostrar resultados
    for alg_name, result in results.items():
        if result["success"]:
            console.print(f"✅ {alg_name}:")
            console.print(f"   Distância: {result['distance']}")
            console.print(f"   Tempo: {result['tempo']:.2f}s")
            console.print(f"   Iterações: {result['iteracoes']}")
            console.print(f"   Centro: {result['center'][:20]}...")
        else:
            console.print(f"❌ {alg_name}: {result['erro']}")
        console.print()

    # Estatísticas finais
    successful = [r for r in results.values() if r["success"]]
    if successful:
        best_distance = min(r["distance"] for r in successful)
        best_algorithm = [
            name
            for name, r in results.items()
            if r["success"] and r["distance"] == best_distance
        ][0]

        console.print(
            f"🏆 Melhor resultado: {best_algorithm} (distância {best_distance})"
        )
        console.print(f"✅ Sucesso: {len(successful)}/{len(results)} algoritmos")
    else:
        console.print("❌ Nenhum algoritmo executado com sucesso")

    console.print("\n🎯 Teste T9 completo!")


if __name__ == "__main__":
    main()
