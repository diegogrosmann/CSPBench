#!/usr/bin/env python3
"""
Teste rápido para verificar paralelismo na análise de sensibilidade
"""

import os
import sys

sys.path.append("/home/diego_grosmann/CSPBench")

from src.utils.worker_calculator import calculate_salib_workers

# Testar com configuração YAML nova
yaml_config = {
    "advanced": {"parallel": {"n_jobs": -1, "internal_workers": 2}}  # Todos os CPUs
}

print("=== TESTE DE CÁLCULO DE WORKERS SALib ===")
print()

# Teste 1: Poucas amostras (5) - antes dava 1 worker
workers_5 = calculate_salib_workers(yaml_config=yaml_config, n_samples=5)
print(f"5 amostras: {workers_5} workers")

# Teste 2: Muitas amostras (20) - deve dar mais workers
workers_20 = calculate_salib_workers(yaml_config=yaml_config, n_samples=20)
print(f"20 amostras: {workers_20} workers")

# Teste 3: Muitas amostras (100) - deve dar máximo de CPUs
workers_100 = calculate_salib_workers(yaml_config=yaml_config, n_samples=100)
print(f"100 amostras: {workers_100} workers")

print()
print("=== RESULTADO ===")
if workers_20 > 1:
    print("✅ Paralelismo SALib FUNCIONANDO!")
else:
    print("❌ Paralelismo SALib ainda com problema")
