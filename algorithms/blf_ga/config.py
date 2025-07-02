"""
Configurações padrão para o algoritmo BLF-GA.

Atributos:
    BLF_GA_DEFAULTS (dict): Parâmetros padrão do BLF-GA.
"""

# BLF-GA Configuration
BLF_GA_DEFAULTS = {
    'pop_size': 150, # Tamanho da população - Valor típico: Entre 50 e 200.
    'initial_blocks': 20, # Número de blocos iniciais
    'min_block_len': 1, # Tamanho mínimo do bloco
    'cross_prob': 0.9, # Probabilidade de crossover - Valor típico: Entre 0.6 e 0.95.
    'mut_prob': 0.1, # Probabilidade de mutação - Valor típico: Entre 0.01 e 0.2.
    'elite_rate': 0.05, # Taxa de elite - Valor típico: Entre 0.01 e 0.1.
    'rediv_freq': 10, # Frequência de redivisão - Valor típico: A cada 10 gerações.
    'max_gens': 400, # Número máximo de gerações - Valor típico: Entre 30 e 100.
    'max_time': 1200.0, # Tempo máximo em segundos - Valor típico: 60 segundos.
    'seed': None, # Semente para reprodutibilidade
    # Novos parâmetros para melhorias
    'immigrant_freq': 10, # Gera imigrantes a cada X gerações
    'immigrant_ratio': 0.2, # Proporção de imigrantes
    'diversity_threshold': 0.4, # Limite para diversidade
    'mutation_adapt_N': 10, # N gerações para detectar convergência
    'mutation_adapt_factor': 2.0, # Fator de aumento temporário da mutação
    'mutation_adapt_duration': 5, # Duração do aumento da mutação
    'mutation_type': 'multi', # multi, inversion, transposition
    'mutation_multi_n': 2, # Número de posições para mutação multi
    'tournament_k': 2, # Parâmetro externo do torneio
    'crossover_type': 'one_point', # one_point, uniform, blend_blocks
    'niching': False, # Ativa niching
    'niching_radius': 3, # Raio de nicho
    'refinement_type': 'greedy', # greedy, swap, insertion, 2opt
    'refine_elites': 'best', # all, best
    'refine_iter_limit': 100, # Limite de iterações por refinamento
    'restart_patience': 20, # Gerações sem melhoria para restart
    'restart_ratio': 0.3, # Proporção da população a reiniciar
    'disable_elitism_gens': 0, # Gerações sem elitismo
}
