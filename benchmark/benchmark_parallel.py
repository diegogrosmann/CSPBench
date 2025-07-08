#!/usr/bin/env python3
"""
Benchmark rápido para testar paralelização do sistema de otimização e análise de sensibilidade.

Este script executa testes pequenos com e sem paralelismo para medir o ganho de performance.
Meta: aceleração ≥ 2x com paralelismo.

Uso:
    python benchmark_parallel.py [--verbose] [--skip-optuna] [--skip-salib]
"""

import argparse
import multiprocessing
import os
import sys
import time
from typing import Any, Dict, List, Tuple

# Adicionar src ao path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.insert(0, project_root)

from src.datasets.dataset_synthetic import generate_dataset_from_params
from src.optimization.optuna_optimizer import optimize_algorithm
from src.optimization.sensitivity_analyzer import analyze_algorithm_sensitivity
from src.utils.worker_calculator import get_cpu_count


def create_test_dataset(size: str = "small") -> Tuple[List[str], str]:
    """Cria dataset de teste baseado no tamanho."""
    if size == "small":
        # Dataset pequeno para testes rápidos
        sequences, _ = generate_dataset_from_params(
            n=8, L=25, alphabet="ACGT", noise=0.1, seed=42
        )
    elif size == "medium":
        # Dataset médio para testes mais realistas
        sequences, _ = generate_dataset_from_params(
            n=15, L=50, alphabet="ACGT", noise=0.15, seed=42
        )
    else:
        raise ValueError(f"Tamanho desconhecido: {size}")

    return sequences, "ACGT"


def benchmark_optuna(
    sequences: List[str], alphabet: str, n_trials: int = 20, verbose: bool = False
) -> Dict[str, float]:
    """Benchmarks para otimização Optuna."""
    results = {}

    print(f"📊 Benchmark Optuna ({n_trials} trials)")

    # Teste serial
    print("  🔄 Executando modo serial...")
    start_time = time.time()
    try:
        result_serial = optimize_algorithm(
            algorithm_name="BLF-GA",
            sequences=sequences,
            alphabet=alphabet,
            n_trials=n_trials,
            timeout_per_trial=30,
            show_progress=verbose,
            yaml_config={"optimization_config": {"parallel": {"n_jobs": 1}}},
        )
        serial_time = time.time() - start_time
        results["serial"] = serial_time
        if verbose:
            print(
                f"    ✅ Serial: {serial_time:.2f}s, melhor={result_serial.best_value:.2f}"
            )
    except Exception as e:
        print(f"    ❌ Erro no modo serial: {e}")
        results["serial"] = float("inf")

    # Teste paralelo
    n_jobs = min(
        get_cpu_count(), 4
    )  # Usar no máximo número de CPUs disponíveis, máximo 4
    print(f"  🚀 Executando modo paralelo ({n_jobs} jobs)...")
    start_time = time.time()
    try:
        result_parallel = optimize_algorithm(
            algorithm_name="BLF-GA",
            sequences=sequences,
            alphabet=alphabet,
            n_trials=n_trials,
            timeout_per_trial=30,
            show_progress=False,  # Desabilitar progresso no paralelo
            yaml_config={"optimization_config": {"parallel": {"n_jobs": n_jobs}}},
        )
        parallel_time = time.time() - start_time
        results["parallel"] = parallel_time
        if verbose:
            print(
                f"    ✅ Paralelo: {parallel_time:.2f}s, melhor={result_parallel.best_value:.2f}"
            )
    except Exception as e:
        print(f"    ❌ Erro no modo paralelo: {e}")
        results["parallel"] = float("inf")

    return results


def benchmark_salib(
    sequences: List[str], alphabet: str, n_samples: int = 40, verbose: bool = False
) -> Dict[str, float]:
    """Benchmarks para análise de sensibilidade SALib."""
    results = {}

    print(f"🔬 Benchmark SALib ({n_samples} amostras)")

    # Teste serial
    print("  🔄 Executando modo serial...")
    start_time = time.time()
    try:
        result_serial = analyze_algorithm_sensitivity(
            algorithm_name="BLF-GA",
            sequences=sequences,
            alphabet=alphabet,
            n_samples=n_samples,
            method="morris",
            timeout_per_sample=15,
            show_progress=verbose,
            yaml_config={"sensitivity_config": {"parallel": {"n_jobs": 1}}},
        )
        serial_time = time.time() - start_time
        results["serial"] = serial_time
        if verbose:
            print(f"    ✅ Serial: {serial_time:.2f}s, análise={result_serial.method}")
    except Exception as e:
        print(f"    ❌ Erro no modo serial: {e}")
        results["serial"] = float("inf")

    # Teste paralelo
    n_jobs = min(
        get_cpu_count(), 4
    )  # Usar no máximo número de CPUs disponíveis, máximo 4
    print(f"  🚀 Executando modo paralelo ({n_jobs} jobs)...")
    start_time = time.time()
    try:
        result_parallel = analyze_algorithm_sensitivity(
            algorithm_name="BLF-GA",
            sequences=sequences,
            alphabet=alphabet,
            n_samples=n_samples,
            method="morris",
            timeout_per_sample=15,
            show_progress=verbose,
            yaml_config={"sensitivity_config": {"parallel": {"n_jobs": n_jobs}}},
        )
        parallel_time = time.time() - start_time
        results["parallel"] = parallel_time
        if verbose:
            print(
                f"    ✅ Paralelo: {parallel_time:.2f}s, análise={result_parallel.method}"
            )
    except Exception as e:
        print(f"    ❌ Erro no modo paralelo: {e}")
        results["parallel"] = float("inf")

    return results


def calculate_speedup(results: Dict[str, float]) -> float:
    """Calcula speedup baseado nos resultados."""
    if "serial" in results and "parallel" in results:
        serial_time = results["serial"]
        parallel_time = results["parallel"]

        if parallel_time > 0 and serial_time > 0:
            return serial_time / parallel_time

    return 0.0


def print_results(
    optuna_results: Dict[str, float],
    salib_results: Dict[str, float],
    system_info: Dict[str, Any],
) -> None:
    """Imprime resultados do benchmark."""
    print(f"\n{'='*60}")
    print("📈 RESULTADOS DO BENCHMARK")
    print(f"{'='*60}")

    print(f"🖥️  Sistema:")
    print(f"   CPUs: {system_info['cpu_count']}")
    print(f"   Python: {system_info['python_version']}")

    # Resultados Optuna
    if optuna_results:
        print(f"\n🔧 Optuna (Otimização):")
        if "serial" in optuna_results:
            print(f"   Serial:   {optuna_results['serial']:.2f}s")
        if "parallel" in optuna_results:
            print(f"   Paralelo: {optuna_results['parallel']:.2f}s")

        speedup = calculate_speedup(optuna_results)
        if speedup > 0:
            print(f"   Speedup:  {speedup:.2f}x")
            if speedup >= 2.0:
                print(f"   ✅ Meta alcançada (≥2x)")
            else:
                print(f"   ⚠️  Meta não alcançada (<2x)")

    # Resultados SALib
    if salib_results:
        print(f"\n🔬 SALib (Análise de Sensibilidade):")
        if "serial" in salib_results:
            print(f"   Serial:   {salib_results['serial']:.2f}s")
        if "parallel" in salib_results:
            print(f"   Paralelo: {salib_results['parallel']:.2f}s")

        speedup = calculate_speedup(salib_results)
        if speedup > 0:
            print(f"   Speedup:  {speedup:.2f}x")
            if speedup >= 2.0:
                print(f"   ✅ Meta alcançada (≥2x)")
            else:
                print(f"   ⚠️  Meta não alcançada (<2x)")

    # Resumo geral
    total_speedups = []
    if optuna_results:
        optuna_speedup = calculate_speedup(optuna_results)
        if optuna_speedup > 0:
            total_speedups.append(optuna_speedup)

    if salib_results:
        salib_speedup = calculate_speedup(salib_results)
        if salib_speedup > 0:
            total_speedups.append(salib_speedup)

    if total_speedups:
        avg_speedup = sum(total_speedups) / len(total_speedups)
        print(f"\n📊 Speedup médio: {avg_speedup:.2f}x")
        if avg_speedup >= 2.0:
            print("🎉 Paralelização bem-sucedida!")
        else:
            print("⚠️  Paralelização pode ser melhorada")

    print(f"\n{'='*60}")


def main():
    """Função principal do benchmark."""
    parser = argparse.ArgumentParser(description="Benchmark de paralelização")
    parser.add_argument("--verbose", "-v", action="store_true", help="Saída detalhada")
    parser.add_argument(
        "--skip-optuna", action="store_true", help="Pular benchmark Optuna"
    )
    parser.add_argument(
        "--skip-salib", action="store_true", help="Pular benchmark SALib"
    )
    parser.add_argument(
        "--dataset-size",
        choices=["small", "medium"],
        default="small",
        help="Tamanho do dataset de teste",
    )

    args = parser.parse_args()

    print("🚀 Iniciando benchmark de paralelização")
    print("=" * 60)

    # Informações do sistema
    cpu_count = get_cpu_count()
    system_info = {
        "cpu_count": cpu_count,
        "python_version": sys.version.split()[0],
    }

    print(f"🖥️  Sistema: {cpu_count} CPUs, Python {system_info['python_version']}")

    # Criar dataset de teste
    print(f"📊 Criando dataset de teste ({args.dataset_size})...")
    sequences, alphabet = create_test_dataset(args.dataset_size)
    print(f"   Dataset: {len(sequences)} sequências, tamanho {len(sequences[0])}")

    # Executar benchmarks
    optuna_results = {}
    salib_results = {}

    try:
        if not args.skip_optuna:
            optuna_results = benchmark_optuna(
                sequences,
                alphabet,
                n_trials=15 if args.dataset_size == "small" else 25,
                verbose=args.verbose,
            )

        if not args.skip_salib:
            salib_results = benchmark_salib(
                sequences,
                alphabet,
                n_samples=30 if args.dataset_size == "small" else 50,
                verbose=args.verbose,
            )

    except KeyboardInterrupt:
        print("\n❌ Benchmark interrompido pelo usuário")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Erro durante benchmark: {e}")
        sys.exit(1)

    # Mostrar resultados
    print_results(optuna_results, salib_results, system_info)


if __name__ == "__main__":
    main()
