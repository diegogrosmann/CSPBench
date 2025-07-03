import csv
import json
from pathlib import Path

def export_batch_json_to_csv(json_path, csv_path):
    """
    Exporta todas as execuções detalhadas do batch (JSON) para CSV.
    Cada linha corresponde a uma execução individual de um algoritmo em uma base.
    """
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    execucoes = data.get('execucoes', [])
    rows = []
    for execucao in execucoes:
        config_nome = execucao.get('config_nome', '')
        bases_info = execucao.get('bases_info', [])
        algoritmos_executados = execucao.get('algoritmos_executados', {})
        # Para cada algoritmo
        for alg, bases in algoritmos_executados.items():
            # Para cada base
            for base_key, result in bases.items():
                base_idx = base_key.replace('base_', '')
                base_params = {}
                for b in bases_info:
                    if str(b.get('base_idx')) == str(base_idx):
                        base_params = b.get('params', {})
                        break
                # Se houver execuções detalhadas, exportar todas
                execucoes_detalhadas = result.get('execucoes_detalhadas')
                if execucoes_detalhadas and isinstance(execucoes_detalhadas, list):
                    for exec_idx, exec_data in enumerate(execucoes_detalhadas, 1):
                        row = {
                            'config_nome': config_nome,
                            'algoritmo': alg,
                            'base_idx': base_idx,
                            'execucao_idx': exec_idx,
                            'dist': exec_data.get('distancia', exec_data.get('melhor_distancia', '')),
                            'dist_base': result.get('dist_base', ''),
                            'tempo': exec_data.get('tempo', ''),
                            'status': exec_data.get('status', result.get('status', '')),
                            'erro': exec_data.get('erro', ''),
                            'melhor_string': exec_data.get('melhor_string', ''),
                            'n': base_params.get('n', ''),
                            'L': base_params.get('L', ''),
                            'alphabet': base_params.get('alphabet', ''),
                            'noise': base_params.get('noise', ''),
                            'seed': base_params.get('seed', ''),
                            'distancia_string_base': base_params.get('distancia_string_base', ''),
                        }
                        rows.append(row)
                else:
                    # Caso antigo: só resultado agregado (manter para compatibilidade)
                    row = {
                        'config_nome': config_nome,
                        'algoritmo': alg,
                        'base_idx': base_idx,
                        'execucao_idx': 1,  # Assumir que é uma única execução
                        'dist': result.get('dist', ''),
                        'dist_base': result.get('dist_base', ''),
                        'tempo': result.get('time', ''),
                        'status': result.get('status', ''),
                        'erro': result.get('erro', ''),
                        'melhor_string': '',
                        'n': base_params.get('n', ''),
                        'L': base_params.get('L', ''),
                        'alphabet': base_params.get('alphabet', ''),
                        'noise': base_params.get('noise', ''),
                        'seed': base_params.get('seed', ''),
                        'distancia_string_base': base_params.get('distancia_string_base', ''),
                    }
                    rows.append(row)
    headers = ['config_nome', 'algoritmo', 'base_idx', 'execucao_idx', 'dist', 'dist_base', 'tempo', 'status', 'erro', 'melhor_string', 'n', 'L', 'alphabet', 'noise', 'seed', 'distancia_string_base']
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
