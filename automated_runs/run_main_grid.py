import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import yaml
import itertools
import subprocess
from datetime import datetime

CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'configs', 'dataset.yaml')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results', 'main_grid')
os.makedirs(RESULTS_DIR, exist_ok=True)

# Carrega o grid de parâmetros do dataset
def load_dataset_grid():
    with open(CONFIG_PATH, 'r') as f:
        grid = yaml.safe_load(f)
    keys = list(grid.keys())
    values = [grid[k] if isinstance(grid[k], list) else [grid[k]] for k in keys]
    return keys, list(itertools.product(*values))

def run_main_with_params(params, keys, idx):
    # Define variáveis de ambiente para o main.py ler
    env = os.environ.copy()
    for k, v in zip(keys, params):
        env[f'CSP_{k.upper()}'] = str(v)
    # Define nome do log/result
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_id = f"run_{idx+1}_{ts}"
    log_path = os.path.join(RESULTS_DIR, f"{run_id}.log")
    # Executa main.py em modo silencioso
    with open(log_path, 'w') as logf:
        proc = subprocess.run(['python3', 'main.py', '--silent'], cwd=os.path.dirname(os.path.dirname(__file__)), env=env, stdout=logf, stderr=subprocess.STDOUT)
    print(f"Execução {run_id} finalizada. Log: {log_path}")

def main():
    keys, combos = load_dataset_grid()
    print(f"Total de execuções: {len(combos)}")
    for idx, params in enumerate(combos):
        print(f"Executando combinação {idx+1}/{len(combos)}: {dict(zip(keys, params))}")
        run_main_with_params(params, keys, idx)

if __name__ == "__main__":
    main()
