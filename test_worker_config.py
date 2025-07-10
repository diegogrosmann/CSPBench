#!/usr/bin/env python3
"""
Teste rápido da configuração de workers
"""

import yaml
from src.utils.worker_calculator import get_worker_config

# Carregar configuração YAML
with open("batch_configs/otimizacao_T1.yaml", "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

# Testar configuração de workers
worker_config = get_worker_config(
    yaml_config=config, context="optuna", algorithm_name="BLF-GA"
)

print("Configuração de Workers:")
print(f"- Total CPUs: {worker_config['total_cpus']}")
print(f"- Workers Optuna (externos): {worker_config['optuna_workers']}")
print(f"- Workers Internos: {worker_config['internal_workers']}")

# Verificar se a configuração está sendo extraída corretamente
print("\nConfiguração do YAML:")
if "advanced" in config:
    advanced = config["advanced"]
    if "parallel" in advanced:
        parallel = advanced["parallel"]
        print(f"- n_jobs: {parallel.get('n_jobs', 'não definido')}")
        print(f"- internal_workers: {parallel.get('internal_workers', 'não definido')}")
    else:
        print("- Seção 'parallel' não encontrada em 'advanced'")
else:
    print("- Seção 'advanced' não encontrada")
