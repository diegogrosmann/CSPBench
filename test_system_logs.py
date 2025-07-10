#!/usr/bin/env python3
"""
Teste simplificado via linha de comando
"""

import subprocess
import time

# Executar o sistema principal e capturar os logs
print("=== EXECUTANDO SISTEMA PRINCIPAL ===")
print("Selecionando: Opção 4 (Batch Unificado) -> Opção 4 (sensibilidade_exemplo.yaml)")
print()

# Simular entrada: 4 (batch unificado) + 4 (sensibilidade)
input_data = "4\n4\n"

# Executar com timeout de 60 segundos
try:
    result = subprocess.run(
        ["/home/diego_grosmann/CSPBench/.venv/bin/python", "main.py"],
        input=input_data,
        text=True,
        timeout=60,
        capture_output=True,
        cwd="/home/diego_grosmann/CSPBench",
    )

    print("STDOUT:")
    print(result.stdout)
    print("\nSTDERR:")
    print(result.stderr)
    print(f"\nReturn code: {result.returncode}")

except subprocess.TimeoutExpired:
    print("✅ Teste executou por 60s - verificar logs para análise de paralelismo")

except Exception as e:
    print(f"❌ Erro na execução: {e}")

# Listar logs recentes
print("\n=== LOGS RECENTES ===")
import os

logs_dir = "/home/diego_grosmann/CSPBench/outputs/logs"
if os.path.exists(logs_dir):
    files = sorted(os.listdir(logs_dir))
    print(f"Log mais recente: {files[-1] if files else 'Nenhum'}")
