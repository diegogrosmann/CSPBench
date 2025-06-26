"""
config.py
=========

Armazena os parâmetros padrão para o algoritmo BLF-GA e para
as entradas de usuário. Este arquivo centraliza todas as configurações
padrão que podem ser referenciadas pelo programa principal.
"""

# Parâmetros gerais
DEBUG_DEFAULT = 'n'  # Modo debug desabilitado por padrão

# Parâmetros para dataset sintético
SYNTHETIC_DEFAULTS = {
    'n': 20,              # Número de strings
    'L': 100,             # Comprimento das strings
    'alphabet': 'ACGT',   # DNA como padrão
    'noise': 0.1          # Taxa de ruído padrão (10%)
}

# Parâmetros para dataset de arquivo
FILE_DEFAULTS = {
    'filepath': 'datasets/sequences.fasta'  # Caminho padrão atualizado
}

# Parâmetros para dataset do NCBI
ENTREZ_DEFAULTS = {
    'email': 'diegogrosmann@gmail.com',                 # Email precisa ser fornecido pelo usuário
    'db': 'nucleotide',                                 # Base padrão
    'term': 'COI[Gene] AND 600:650[SLEN]',              # Termo de busca padrão
    'n': 20,                                            # Número de registros padrão
    'api_key': '40aff1a7f51fc7e0711203ea3b2f5ae37c09'   # API key do NCBI (adicione aqui se disponível)
}

# Parâmetros do algoritmo BLF-GA
BLFGA_DEFAULTS = {
    'pop_size': 50,        # Tamanho da população
    'initial_blocks': 5,   # Número inicial de blocos
    'min_block_len': 3,    # Tamanho mínimo dos blocos
    'cross_prob': 0.9,     # Probabilidade de crossover
    'mut_prob': 0.5,       # Probabilidade de mutação
    'elite_rate': 0.05,    # Taxa de elitismo
    'rediv_freq': 10,      # Frequência de redivisão
    'max_gens': 300,       # Número máximo de gerações
    'max_time': 60.0,      # Tempo máximo de execução (segundos)
    'seed': 0              # Semente para o gerador de números aleatórios
}
