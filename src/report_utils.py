"""
Utilitários para exibição de resumo e salvamento de relatórios detalhados.

Funções:
    print_quick_summary(...): Exibe resumo rápido dos resultados no console.
    save_detailed_report(...): Salva relatório detalhado usando ResultsFormatter.
"""

def print_quick_summary(baseline_val, baseline_time, results, console):
    console.print("\n" + "="*50)
    console.print("RESUMO RÁPIDO")
    console.print("="*50)
    console.print(f"{'Algoritmo':<12} {'Melhor Dist':<12} {'Tempo Médio':<12} {'Gap%':<6}")
    console.print("-" * 45)
    console.print(f"{'Baseline':<12} {baseline_val:<12} {baseline_time:<12.3f} {'-':<6}")
    for alg_name, res in results.items():
        dist_str = str(res.get('dist', '-'))
        time_str = f"{res.get('time', 0):.3f}"
        gap_str = f"{res.get('gap', 0):.1f}" if isinstance(res.get('gap', 0), (int, float)) else "-"
        console.print(f"{alg_name:<12} {dist_str:<12} {time_str:<12} {gap_str:<6}")

def save_detailed_report(formatter, filename):
    formatter.save_detailed_report(filename)
