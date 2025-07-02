"""
Configurações padrão para o algoritmo BLF-GA.

Atributos:
    BLF_GA_DEFAULTS (dict): Parâmetros padrão do BLF-GA.
"""

# BLF-GA Configuration
BLF_GA_DEFAULTS = {
    'pop_size': 100,
    'initial_blocks': 20,
    'min_block_len': 3,
    'cross_prob': 2,
    'mut_prob': 0.9,
    'elite_rate': 0.05,
    'rediv_freq': 10,
    'max_gens': 30,
    'max_time': 60.0,
    'seed': None,
}
