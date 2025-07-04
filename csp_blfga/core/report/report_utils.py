"""
Utilitários para exibição de resumo de resultados.

Funções:
    print_quick_summary(results, console): Exibe resumo rápido dos resultados no console.
"""


def print_quick_summary(results, console):
    """
    Exibe um resumo rápido dos resultados dos algoritmos no console.

    Args:
        results (dict): Dicionário com resultados por algoritmo.
        console: Instância do gerenciador de console para saída.
    """
    console.print("\n" + "=" * 60)
    console.print("RESUMO RÁPIDO")
    console.print("=" * 60)
    console.print(
        f"{'Algoritmo':<12} {'Melhor Dist':<12} {'Dist. String Base':<16} {'Tempo Médio':<12}"
    )
    console.print("-" * 55)
    for alg_name, res in results.items():
        dist_str = str(res.get("dist", "-"))
        dist_base_str = str(res.get("dist_base", "-"))
        time_str = f"{res.get('time', 0):.3f}"
        console.print(
            f"{alg_name:<12} {dist_str:<12} {dist_base_str:<16} {time_str:<12}"
        )
