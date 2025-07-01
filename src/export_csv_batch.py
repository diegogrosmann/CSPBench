import csv
import json
from pathlib import Path

def export_batch_json_to_csv(json_path, csv_path):
    """
    Exporta os resultados detalhados do batch (JSON) para CSV.
    """
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    execucoes = data.get('execucoes', [])
    rows = []
    for execucao in execucoes:
        config_nome = execucao.get('config_nome', '')
        bases_info = execucao.get('bases_info', [])
        algoritmos_executados = execucao.get('algoritmos_executados', {})
        for alg, bases in algoritmos_executados.items():
            for base_key, result in bases.items():
                base_idx = base_key.replace('base_', '')
                base_params = {}
                for b in bases_info:
                    if str(b.get('base_idx')) == str(base_idx):
                        base_params = b.get('params', {})
                        break
                row = {
                    'config_nome': config_nome,
                    'algoritmo': alg,
                    'base_idx': base_idx,
                    'dist': result.get('dist', ''),
                    'dist_base': result.get('dist_base', ''),
                    'tempo': result.get('time', ''),
                    'status': result.get('status', ''),
                    'n': base_params.get('n', ''),
                    'L': base_params.get('L', ''),
                    'alphabet': base_params.get('alphabet', ''),
                    'noise': base_params.get('noise', ''),
                    'seed': base_params.get('seed', ''),
                }
                rows.append(row)
    headers = ['config_nome', 'algoritmo', 'base_idx', 'dist', 'dist_base', 'tempo', 'status', 'n', 'L', 'alphabet', 'noise', 'seed']
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
