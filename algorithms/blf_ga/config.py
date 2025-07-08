"""
Configurações padrão para o algoritmo BLF-GA.

O BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) é uma metaheurística híbrida
que combina aprendizado local por blocos com algoritmo genético global.

PARÂMETROS PRINCIPAIS:
- População: Controla tamanho e qualidade da busca
- Blocos: Define granularidade do aprendizado local
- Operadores: Configura evolução genética
- Adaptação: Mecanismos dinâmicos de ajuste
- Critérios: Controla quando parar/reiniciar

RECOMENDAÇÕES DE USO:
- Instâncias pequenas (n<10): pop_size=50-100, max_gens=50-100
- Instâncias médias (n=10-50): pop_size=100-200, max_gens=100-200
- Instâncias grandes (n>50): pop_size=200-500, max_gens=200-500

Atributos:
    BLF_GA_DEFAULTS (dict): Parâmetros padrão do BLF-GA.
"""

# BLF-GA Configuration
BLF_GA_DEFAULTS = {
    # --- PARÂMETROS DE POPULAÇÃO ---
    "pop_size": 1.5,  # Tamanho da população (int fixo ou float >=1 para multiplicador de n)
    # Exemplos: 100 (fixo), 1.5 (1.5*n), 2.0 (2*n)
    "min_pop_size": 20,  # Tamanho mínimo da população quando usar proporção de n
    # Evita populações muito pequenas em instâncias menores
    "seed": None,  # Semente para reprodutibilidade (None = aleatório)
    # --- PARÂMETROS DE BLOCOS (LEARNING) ---
    "initial_blocks": 0.2,  # Número de blocos iniciais (int fixo ou float 0-1 para proporção de L)
    # Exemplos: 5 (fixo), 0.2 (20% do comprimento)
    "min_block_len": 1,  # Tamanho mínimo do bloco (evita blocos muito pequenos)
    "rediv_freq": 10,  # Frequência de redivisão (a cada X gerações)
    # Permite adaptação dinâmica da estrutura de blocos
    # --- OPERADORES GENÉTICOS ---
    "cross_prob": 0.9,  # Probabilidade de crossover (0.6 a 0.95)
    # Maior valor = mais recombinação
    "crossover_type": "one_point",  # Tipo de crossover: one_point, uniform, blend_blocks
    "mut_prob": 0.1,  # Probabilidade de mutação (0.01 a 0.2)
    # Balança exploração vs exploração
    "mutation_type": "multi",  # Tipo de mutação: multi, inversion, transposition
    "mutation_multi_n": 2,  # Número de posições para mutação multi
    "elite_rate": 0.05,  # Taxa de elite (0.01 a 0.1) - % dos melhores preservados
    "tournament_k": 2,  # Tamanho do torneio para seleção
    # --- DIVERSIDADE E IMIGRANTES ---
    "immigrant_freq": 10,  # Gera imigrantes a cada X gerações
    # Injeta diversidade substituindo piores por aleatórios
    "immigrant_ratio": 0.2,  # Proporção de imigrantes (% da população)
    "diversity_threshold": 0.4,  # Limite para diversidade (0-1)
    # Abaixo deste valor, ativa mutação adaptativa
    # --- MUTAÇÃO ADAPTATIVA ---
    "mutation_adapt_N": 10,  # N gerações para detectar convergência
    "mutation_adapt_factor": 2.0,  # Fator de aumento temporário da mutação
    "mutation_adapt_duration": 5,  # Duração do aumento da mutação (gerações)
    # --- NICHING (NICHOS) ---
    "niching": False,  # Ativa niching (preserva diversidade local)
    "niching_radius": 3,  # Raio de nicho (distância mínima entre soluções)
    # --- REFINAMENTO LOCAL ---
    "refinement_type": "greedy",  # Tipo: greedy, swap, insertion, 2opt
    "refine_elites": "best",  # Quem refinar: all (todos elites), best (apenas melhor)
    "refine_iter_limit": 100,  # Limite de iterações por refinamento
    # --- CRITÉRIOS DE PARADA E REINÍCIO ---
    "max_gens": 100,  # Número máximo de gerações (30 a 100+)
    "max_time": 1200.0,  # Tempo máximo em segundos (20 minutos)
    "no_improve_patience": 0.2,  # Gerações sem melhoria para early stopping
    # (int fixo ou float 0-1 para proporção de max_gens)
    "restart_patience": 20,  # Gerações sem melhoria para restart parcial
    "restart_ratio": 0.3,  # Proporção da população a reiniciar
    "disable_elitism_gens": 5,  # Desabilita elitismo a cada X gerações
    # Evita convergência prematura
}
