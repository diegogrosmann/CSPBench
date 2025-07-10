#!/usr/bin/env python3
"""
Teste rápido da paralelização corrigida
"""

import yaml
from src.optimization.optuna_optimizer import optimize_algorithm

# Dados de teste pequenos
sequences = ["ACGT", "AAGT", "ACTT"]
alphabet = "ACGT"

# Carregar configuração YAML
with open("batch_configs/otimizacao_T1.yaml", "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

print("=== TESTE DE PARALELIZAÇÃO CORRIGIDA ===")
print(f"Configuração YAML carregada:")
print(f"  - n_jobs: {config['advanced']['parallel']['n_jobs']}")
print(f"  - internal_workers: {config['advanced']['parallel']['internal_workers']}")

print(f"\nTeste rápido com 2 trials, timeout 30s...")

try:
    result = optimize_algorithm(
        algorithm_name="BLF-GA",
        sequences=sequences,
        alphabet=alphabet,
        n_trials=2,  # Apenas 2 trials para teste rápido
        timeout_per_trial=30,  # 30 segundos por trial
        show_progress=True,
        yaml_config=config,  # Incluir configuração YAML
    )

    print(f"\n✅ Teste concluído com sucesso!")
    print(f"  - Melhor valor: {result.best_value}")
    print(f"  - Trials executados: {result.n_trials}")
    print(f"  - Tempo total: {result.optimization_time:.2f}s")

except Exception as e:
    print(f"\n❌ Erro no teste: {e}")
