"""
config.py
=========

Parâmetros padrão e configurações centrais do projeto.
"""

# --------------------------------------------------
# Geral
# --------------------------------------------------
DEBUG_DEFAULT = 'n'            # Modo debug desabilitado por padrão

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
