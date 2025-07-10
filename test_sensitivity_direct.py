#!/usr/bin/env python3
"""
Teste direto da análise de sensibilidade com paralelismo
"""

import os
import sys

sys.path.append("/home/diego_grosmann/CSPBench")

from algorithms.blf_ga.algorithm import BLFGA
from src.optimization.sensitivity_analyzer import analyze_algorithm_sensitivity

print("=== TESTE DIRETO ANÁLISE DE SENSIBILIDADE ===")

# Criar sequências simples para teste
sequences = ["ACGTACGT", "ACGTACCG", "ACGTAAGT", "ACGGACGT"]

print(f"Dataset teste: {len(sequences)} sequências")

# Configurar análise de sensibilidade com paralelismo
yaml_config = {
    "advanced": {
        "parallel": {"n_jobs": 4, "internal_workers": 2}  # 4 workers explícitos
    }
}

sensitivity_config = {
    "n_samples": 8,  # Poucas amostras para teste rápido
    "method": "morris",
    "timeout_per_sample": 10,
}

print(
    f"Configurado para usar paralelismo com {yaml_config['advanced']['parallel']['n_jobs']} workers"
)

# Executar análise
try:
    result = analyze_algorithm_sensitivity(
        algorithm_name="BLF-GA",
        algorithm_class=BLFGA,
        sequences=sequences,
        alphabet="ACGT",
        analysis_config=sensitivity_config,
        yaml_config=yaml_config,
    )
    print("✅ Análise concluída com sucesso!")
    print(f"Resultados: {type(result)}")
except Exception as e:
    print(f"❌ Erro na análise: {e}")
    import traceback

    traceback.print_exc()
