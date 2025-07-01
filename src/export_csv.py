import csv
from pathlib import Path


def export_results_to_csv(formatter, filename):
    """
    Exporta todos os dados de execuções detalhadas para um arquivo CSV.
    Cada linha corresponde a uma execução de um algoritmo.
    """
    file_path = Path(filename)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Cabeçalhos padrão
    headers = [
        'algoritmo', 'execucao', 'melhor_string', 'distancia', 'tempo', 'erro', 'distancia_string_base', 'seed'
    ]
    with open(file_path, 'w', encoding='utf-8', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for alg_name, executions in formatter.results.items():
            for idx, exec_data in enumerate(executions, 1):
                row = {
                    'algoritmo': alg_name,
                    'execucao': idx,
                    'melhor_string': exec_data.get('melhor_string', ''),
                    'distancia': exec_data.get('distancia', exec_data.get('melhor_distancia', '')),
                    'tempo': exec_data.get('tempo', ''),
                    'erro': exec_data.get('erro', ''),
                    'distancia_string_base': exec_data.get('distancia_string_base', ''),
                    'seed': exec_data.get('seed', '')
                }
                writer.writerow(row)
