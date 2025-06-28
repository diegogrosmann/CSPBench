"""
config.py
=========

Parâmetros padrão e configurações centrais do projeto.

Atributos:
    DEBUG_DEFAULT (str): Modo debug padrão ('n' = desabilitado).
    ALGORITHM_TIMEOUT (int): Tempo limite em segundos para execução de algoritmos.
    SYNTHETIC_DEFAULTS (dict): Parâmetros padrão para datasets sintéticos.
    FILE_DEFAULTS (dict): Parâmetros padrão para datasets de arquivo.
    ENTREZ_DEFAULTS (dict): Parâmetros padrão para datasets NCBI.
"""

import sys

# --------------------------------------------------
# Geral
# --------------------------------------------------
DEBUG_DEFAULT = 'n'            # Modo debug desabilitado por padrão
ALGORITHM_TIMEOUT = 300        # Timeout padrão de 5 minutos para algoritmos

def safe_input(prompt: str, default: str = "") -> str:
    """Input seguro que trata KeyboardInterrupt de forma consistente."""
    try:
        return input(prompt).strip()
    except (KeyboardInterrupt, EOFError):
        print("\nOperação cancelada pelo usuário.")
        sys.exit(0)

# --------------------------------------------------
# Dataset sintético
# --------------------------------------------------
SYNTHETIC_DEFAULTS = {
    'n': 20,
    'L': 100,
    'alphabet': 'ACGT',
    'noise': 0.10,
}

# --------------------------------------------------
# Dataset de arquivo
# --------------------------------------------------
FILE_DEFAULTS = {
    'filepath': 'saved_datasets/sequences.fasta',
}

# --------------------------------------------------
# Dataset NCBI (Bio.Entrez)
# --------------------------------------------------
ENTREZ_DEFAULTS = {
    'email':  'diegogrosmann@gmail.com',
    'db':     'nucleotide',
    'term':   'COI[Gene] AND 600:650[SLEN]',
    'n':      20,
    'api_key': '40aff1a7f51fc7e0711203ea3b2f5ae37c09',
}
