#!/usr/bin/env python3
"""
Teste rápido dos logs de paralelismo
"""

import logging
import yaml
from src.utils.worker_calculator import get_worker_config

# Configurar logging para DEBUG
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)

# Carregar configuração YAML
with open("batch_configs/otimizacao_T1.yaml", "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

print("=== TESTE DE CONFIGURAÇÃO DE WORKERS ===")

# Testar configuração de workers Optuna
worker_config_optuna = get_worker_config(
    yaml_config=config, context="optuna", algorithm_name="BLF-GA"
)

print(f"\n1. Optuna Workers:")
print(f"   - Total CPUs: {worker_config_optuna['total_cpus']}")
print(f"   - Workers Optuna (externos): {worker_config_optuna['optuna_workers']}")
print(f"   - Workers Internos: {worker_config_optuna['internal_workers']}")

# Testar configuração de workers SALib
worker_config_salib = get_worker_config(
    yaml_config=config, context="salib", algorithm_name="BLF-GA", n_samples=1000
)

print(f"\n2. SALib Workers:")
print(f"   - Total CPUs: {worker_config_salib['total_cpus']}")
print(f"   - Workers SALib: {worker_config_salib['salib_workers']}")

print(f"\n3. Configuração YAML extraída:")
if "advanced" in config:
    advanced = config["advanced"]
    if "parallel" in advanced:
        parallel = advanced["parallel"]
        print(f"   - n_jobs: {parallel.get('n_jobs', 'não definido')}")
        print(
            f"   - internal_workers: {parallel.get('internal_workers', 'não definido')}"
        )
        print(f"   - log_level: {advanced.get('log_level', 'não definido')}")

print("\n✅ Teste concluído - configuração funcionando!")
