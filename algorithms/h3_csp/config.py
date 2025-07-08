"""
Configurações padrão para o algoritmo H³-CSP.

Este módulo define todos os parâmetros padrão do algoritmo H³-CSP,
incluindo parâmetros de divisão de blocos, técnicas de busca,
refinamento e controle de execução.

Constantes:
    H3_CSP_DEFAULTS (dict): Dicionário com todos os parâmetros padrão
                           do algoritmo H³-CSP.
"""

# H³-CSP Configuration
H3_CSP_DEFAULTS = {
    # === Parâmetros de Divisão de Blocos ===
    "auto_blocks": True,  # Usa divisão automática por √L
    "min_block_size": 2,  # Tamanho mínimo de bloco
    "max_blocks": None,  # Máximo de blocos (None = automático)
    "block_size": 2,  # Tamanho base de bloco (se não automático)
    "block_strategy": None,  # Estratégia de divisão (None = padrão)
    # === Limiares de Dificuldade por Bloco ===
    "block_small": 2,  # Limite para blocos "pequenos" (busca exaustiva)
    "block_medium": 4,  # Limite para blocos "médios" (beam search reduzido)
    "block_large": 8,  # Limite para blocos "grandes" (beam search completo)
    # === Parâmetros de Busca ===
    "exhaustive_limit": 10000,  # Limite para busca exaustiva (|Σ|^m)
    "beam_width": 32,  # Largura do beam search
    "k_candidates": 5,  # Número de candidatos por bloco
    # === Refinamento Global ===
    "local_search_iters": 3,  # Iterações de busca local (deprecated)
    "local_iters": 3,  # Iterações de busca local (atual)
    # === Controle de Execução ===
    "max_time": 300,  # Tempo máximo em segundos
    "seed": None,  # Semente para reprodutibilidade
    # === Parâmetros Experimentais ===
    "diversity_threshold": 1,  # Limiar de diversidade (não usado atualmente)
    "fallback_enabled": True,  # Habilita fallback para blocos grandes
}

# Mapeamento de parâmetros para compatibilidade
# Alguns parâmetros têm nomes alternativos para compatibilidade com versões anteriores
_PARAMETER_ALIASES = {
    "local_search_iters": "local_iters",  # Alias para compatibilidade
}
