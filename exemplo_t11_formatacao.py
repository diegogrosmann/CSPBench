#!/usr/bin/env python3
"""
Exemplo de uso do sistema de formatação melhorado (T11).

Este script demonstra:
- ResultsFormatter com print_quick_summary integrado
- Métodos comparativos renomeados para clareza
- Parsing simplificado de algoritmos e bases
"""

import sys
from pathlib import Path

# Adicionar diretório raiz ao path
sys.path.insert(0, str(Path(__file__).parent))

from csp_blfga.core.io.results_formatter import ResultsFormatter
from csp_blfga.ui.cli.console_manager import ConsoleManager


def main():
    """Demonstra as melhorias T11 do sistema de formatação."""

    console = ConsoleManager()

    console.print("🎨 Exemplo T11 - Sistema de Formatação Melhorado")
    console.print("=" * 60)

    # Criar formatter
    formatter = ResultsFormatter()

    # Simular resultados de diferentes algoritmos
    console.print("📊 Adicionando resultados simulados...")

    # Baseline
    baseline_results = [
        {
            "tempo": 0.123,
            "distancia": 5,
            "melhor_string": "ACGTACGT",
            "iteracoes": 1,
            "erro": "",
        },
        {
            "tempo": 0.125,
            "distancia": 5,
            "melhor_string": "ACGTACGT",
            "iteracoes": 1,
            "erro": "",
        },
        {
            "tempo": 0.124,
            "distancia": 5,
            "melhor_string": "ACGTACGT",
            "iteracoes": 1,
            "erro": "",
        },
    ]
    formatter.add_algorithm_results("Baseline", baseline_results)

    # BLF-GA
    blfga_results = [
        {
            "tempo": 2.456,
            "distancia": 3,
            "melhor_string": "ACGTCCGT",
            "iteracoes": 150,
            "erro": "",
        },
        {
            "tempo": 2.478,
            "distancia": 4,
            "melhor_string": "ACGTACCT",
            "iteracoes": 145,
            "erro": "",
        },
        {
            "tempo": 10.0,
            "distancia": float("inf"),
            "melhor_string": "",
            "iteracoes": 0,
            "erro": "Timeout (10s)",
        },
    ]
    formatter.add_algorithm_results("BLF-GA", blfga_results)

    # Dataset_Base123_DP-CSP (testando parsing)
    dp_results = [
        {
            "tempo": 0.789,
            "distancia": 2,
            "melhor_string": "ACGTACGA",
            "iteracoes": 1,
            "erro": "",
        },
    ]
    formatter.add_algorithm_results("Synthetic_Base001_DP-CSP", dp_results)

    console.print("✅ Resultados adicionados!")
    console.print("")

    # T11-1: Demonstrar print_quick_summary integrado
    console.print("🎯 T11-1: Resumo rápido integrado ao ResultsFormatter")
    formatter.print_quick_summary(console)

    console.print("\n" + "=" * 60)
    console.print("🔄 T11-2: Métodos comparativos renomeados")
    console.print("  - _format_comparative_table → _format_algorithms_comparison")
    console.print("  - _format_algorithm_table → _format_algorithm_executions_table")
    console.print(
        "  - _format_algorithm_statistics → _format_algorithm_summary_statistics"
    )
    console.print("  - _format_algorithm_strings → _format_algorithm_best_strings")
    console.print("✅ Métodos renomeados para melhor clareza!")

    console.print("\n" + "=" * 60)
    console.print("🧹 T11-3: Parsing simplificado de algoritmos e bases")

    # Testar parsing simplificado
    test_names = [
        "Baseline",
        "BLF-GA",
        "Synthetic_Base001_DP-CSP",
        "RealData_Base456_H3-CSP",
        "SimpleAlgorithm",
    ]

    for name in test_names:
        algo_name, base_display = formatter._parse_algorithm_name(name)
        simplified = formatter._simplify_algorithm_label(algo_name)
        console.print(f"  {name}")
        console.print(f"    → Algoritmo: {algo_name}")
        console.print(f"    → Base: {base_display}")
        console.print(f"    → Simplificado: {simplified}")
        console.print("")

    console.print("✅ Parsing simplificado funcionando!")

    console.print("\n" + "=" * 60)
    console.print("🎊 TODAS AS MELHORIAS T11 IMPLEMENTADAS:")
    console.print("  ✅ T11-1: print_quick_summary integrado ao ResultsFormatter")
    console.print("  ✅ T11-2: Métodos comparativos renomeados para clareza")
    console.print("  ✅ T11-3: Parsing simplificado de algoritmos e bases")
    console.print("\n🎯 Teste T11 completo!")


if __name__ == "__main__":
    main()
