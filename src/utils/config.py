"""
Parâmetros padrão e configurações centrais do projeto CSP.

Atributos:
    DEBUG_DEFAULT (str): Modo debug padrão ('n' = desabilitado).
    ALGORITHM_TIMEOUT (int): Timeout padrão para execução de algoritmos.
    SYNTHETIC_DEFAULTS (dict): Parâmetros padrão para datasets sintéticos.
    FILE_DEFAULTS (dict): Parâmetros padrão para datasets de arquivo.
    ENTREZ_DEFAULTS (dict): Parâmetros padrão para datasets NCBI.

Funções:
    safe_input(prompt, default): Entrada segura, compatível com modo automatizado.
"""

# config.py
# =========

import os
import sys

# --------------------------------------------------
# Geral
# --------------------------------------------------
DEBUG_DEFAULT = "n"  # Modo debug desabilitado por padrão
ALGORITHM_TIMEOUT = 300  # Timeout padrão de 5 minutos para algoritmos


def safe_input(prompt: str, default: str = "") -> str:
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        # Retorna sempre o default ou vazio para todos os prompts
        return default
    try:
        return input(prompt).strip()
    except (KeyboardInterrupt, EOFError):
        print("\nOperação cancelada pelo usuário.")
        sys.exit(0)
    return default  # Garante retorno em todos os caminhos


# --------------------------------------------------
# Dataset sintético
# --------------------------------------------------
SYNTHETIC_DEFAULTS = {
    "n": 20,
    "L": 100,
    "alphabet": "ACGT",
    "noise": 0.10,
}

# --------------------------------------------------
# Dataset de arquivo
# --------------------------------------------------
FILE_DEFAULTS = {
    "filepath": "saved_datasets/sequences.fasta",
}

# --------------------------------------------------
# Dataset NCBI (Bio.Entrez)
# --------------------------------------------------
ENTREZ_DEFAULTS = {
    "email": "diegogrosmann@gmail.com",
    "db": "nucleotide",
    "term": "COI[Gene] AND 600:650[SLEN]",
    "n": 20,
    "api_key": "40aff1a7f51fc7e0711203ea3b2f5ae37c09",
}

# --------------------------------------------------
# Execução em lote
# --------------------------------------------------
BATCH_DEFAULTS = {
    "timeout_global": 1800,  # 30 minutos por execução
    "max_concurrent": 1,  # Execuções paralelas (futuro)
    "save_individual_reports": True,
    "save_consolidated_report": True,
}
