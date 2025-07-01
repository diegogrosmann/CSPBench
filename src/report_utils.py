"""
Utilitários para exibição de resumo e salvamento de relatórios detalhados.

Funções:
    print_quick_summary(results, console): Exibe resumo rápido dos resultados no console.
    save_detailed_report(formatter, filename): Salva relatório detalhado usando ResultsFormatter.
"""

def print_quick_summary(results, console):
    """
    Exibe um resumo rápido dos resultados dos algoritmos no console.

    Args:
        results (dict): Dicionário com resultados por algoritmo.
        console: Instância do gerenciador de console para saída.
    """
    console.print("\n" + "="*60)
    console.print("RESUMO RÁPIDO")
    console.print("="*60)
    console.print(f"{'Algoritmo':<12} {'Melhor Dist':<12} {'Dist. String Base':<16} {'Tempo Médio':<12}")
    console.print("-" * 55)
    for alg_name, res in results.items():
        dist_str = str(res.get('dist', '-'))
        dist_base_str = str(res.get('dist_base', '-'))
        time_str = f"{res.get('time', 0):.3f}"
        console.print(f"{alg_name:<12} {dist_str:<12} {dist_base_str:<16} {time_str:<12}")

def _convert_dict_keys_to_str(obj):
    """
    Recursivamente converte todas as chaves de dicionários para string.

    Args:
        obj: Objeto a ser convertido.
    Returns:
        Objeto com todas as chaves de dicionário como string.
    """
    if isinstance(obj, dict):
        return {str(k): _convert_dict_keys_to_str(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_convert_dict_keys_to_str(i) for i in obj]
    else:
        return obj

def save_detailed_report(formatter, filename):
    """
    Salva relatório detalhado dos resultados em arquivo JSON ou formato antigo.

    Args:
        formatter: Instância de ResultsFormatter.
        filename (str): Caminho do arquivo de saída.
    """
    # Antes de salvar, converter todas as chaves de dicionário para string
    data = formatter.get_detailed_report_data() if hasattr(formatter, "get_detailed_report_data") else None
    if data is not None:
        data = _convert_dict_keys_to_str(data)
        import json
        with open(filename, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
    else:
        # fallback para método antigo
        formatter.save_detailed_report(filename)
