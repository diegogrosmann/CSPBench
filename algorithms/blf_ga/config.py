"""
Configurações padrão para o algoritmo BLF-GA.

Atributos:
    BLF_GA_DEFAULTS (dict): Parâmetros padrão do BLF-GA.
"""

# BLF-GA Configuration
BLF_GA_DEFAULTS = {
    # --- Parâmetros de População ---
    "pop_size": 1.5,  # Tamanho da população (int fixo ou float >=1 para multiplicador de n)
    "seed": None,  # Semente para reprodutibilidade
    # --- Parâmetros de Blocos ---
    "initial_blocks": 0.2,  # Número de blocos iniciais (int fixo ou float 0-1 para proporção de L)
    "min_block_len": 1,  # Tamanho mínimo do bloco
    "rediv_freq": 10,  # Frequência de redivisão (a cada X gerações)
    # --- Operadores Genéticos ---
    "cross_prob": 0.9,  # Probabilidade de crossover (0.6 a 0.95)
    "crossover_type": "one_point",  # Tipo de crossover: one_point, uniform, blend_blocks
    "mut_prob": 0.1,  # Probabilidade de mutação (0.01 a 0.2)
    "mutation_type": "multi",  # Tipo de mutação: multi, inversion, transposition
    "mutation_multi_n": 2,  # Número de posições para mutação multi
    "elite_rate": 0.05,  # Taxa de elite (0.01 a 0.1)
    "tournament_k": 2,  # Parâmetro externo do torneio
    # --- Diversidade e Imigrantes ---
    "immigrant_freq": 10,  # Gera imigrantes a cada X gerações
    "immigrant_ratio": 0.2,  # Proporção de imigrantes
    "diversity_threshold": 0.4,  # Limite para diversidade
    "mutation_adapt_N": 10,  # N gerações para detectar convergência
    "mutation_adapt_factor": 2.0,  # Fator de aumento temporário da mutação
    "mutation_adapt_duration": 5,  # Duração do aumento da mutação
    # --- Niching (Nichos) ---
    "niching": False,  # Ativa niching
    "niching_radius": 3,  # Raio de nicho
    # --- Refinamento Local ---
    "refinement_type": "greedy",  # greedy, swap, insertion, 2opt
    "refine_elites": "best",  # all, best
    "refine_iter_limit": 100,  # Limite de iterações por refinamento
    # --- Critérios de Parada e Reinício ---
    "max_gens": 400,  # Número máximo de gerações (30 a 100+)
    "max_time": 1200.0,  # Tempo máximo em segundos
    "no_improve_patience": 0.2,  # Gerações sem melhoria para encerrar (int fixo ou float 0-1 para proporção de max_gens)
    "restart_patience": 20,  # Gerações sem melhoria para restart
    "restart_ratio": 0.3,  # Proporção da população a reiniciar
    "disable_elitism_gens": 5,  # Gerações sem elitismo
}
