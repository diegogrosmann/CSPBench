import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import yaml
import pandas as pd
from itertools import product
from datasets.dataset_synthetic import generate_dataset_with_params
from algorithms.blf_ga.config import BLF_GA_DEFAULTS
from algorithms.blf_ga.algorithm import BLFGAAlgorithm

CONFIG_DIR = os.path.join(os.path.dirname(__file__), 'configs')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

def load_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Carregar parâmetros do dataset
    dataset_params_config = load_yaml(os.path.join(CONFIG_DIR, 'dataset.yaml'))
    # Carregar grid de parâmetros do BLF-GA
    blfga_param_grid = load_yaml(os.path.join(CONFIG_DIR, 'blfga_grid.yaml'))
    if blfga_param_grid is None:
        blfga_param_grid = {}

    # Gerar todas as combinações de parâmetros de dataset
    dataset_param_names = list(dataset_params_config.keys())
    dataset_param_values = [
        val if isinstance(val, list) else [val] 
        for val in dataset_params_config.values()
    ]
    dataset_configs = [
        dict(zip(dataset_param_names, values)) 
        for values in product(*dataset_param_values)
    ]
    
    print(f"Parâmetros usados: n={dataset_params_config.get('n', 'N/A')}, L={dataset_params_config.get('L', 'N/A')}, alphabet='{dataset_params_config.get('alphabet', 'N/A')}', noise={dataset_params_config.get('noise', 'N/A')}, fully_random={dataset_params_config.get('fully_random', 'N/A')}")

    all_results = []

    for ds_idx, ds_params in enumerate(dataset_configs):
        # Gerar dataset sintético
        strings, params_usados = generate_dataset_with_params(ds_params)
        print(f"\n--- Dataset {ds_idx+1}/{len(dataset_configs)}: n={len(strings)}, L={len(strings[0])}, |Σ|={len(params_usados['alphabet'])} ---")

        # Gerar grid de experimentos para BLF-GA
        param_names = list(BLF_GA_DEFAULTS.keys())
        param_values = [
            blfga_param_grid.get(k, [BLF_GA_DEFAULTS[k]]) 
            if isinstance(blfga_param_grid.get(k, BLF_GA_DEFAULTS[k]), list) 
            else [blfga_param_grid.get(k, BLF_GA_DEFAULTS[k])] 
            for k in param_names
        ]
        blfga_experimentos = list(product(*param_values))
        print(f"Total de configurações BLF-GA: {len(blfga_experimentos)}")

        resultados_dataset = []
        for i, valores in enumerate(blfga_experimentos):
            params = {k: v for k, v in zip(param_names, valores)}
            alg = BLFGAAlgorithm(strings, params_usados['alphabet'], **params)
            import time
            t0 = time.time()
            center, dist, history = alg.run_with_history()
            t1 = time.time()
            
            result_row = {
                'dataset_n': len(strings),
                'dataset_L': len(strings[0]),
                **params,
                'dist': dist,
                'tempo': t1-t0,
                'history': history
            }
            resultados_dataset.append(result_row)
            print(f"Experimento {i+1}/{len(blfga_experimentos)}: dist={dist}, tempo={t1-t0:.2f}s, params={params}")
        
        all_results.extend(resultados_dataset)

    # Salvar resultados consolidados
    if all_results:
        df = pd.DataFrame(all_results)
        csv_path = os.path.join(RESULTS_DIR, 'blfga_experimentos.csv')
        df.to_csv(csv_path, index=False)
        print(f"\nResultados salvos em: {csv_path}")
    else:
        print("Nenhum resultado foi gerado.")

if __name__ == "__main__":
    main()
