# H³-CSP Configuration
H3_CSP_DEFAULTS = {
    'auto_blocks': True,
    'min_block_size': 2,
    'max_blocks': None,  # pode ser L/2
    'exhaustive_limit': 10000,
    'beam_width': 32,
    'k_candidates': 5,
    'local_search_iters': 3,
    'max_time': 300,
    'diversity_threshold': 1,
    'fallback_enabled': True,
    'seed': 42,
    # Adicionados para evitar KeyError e garantir tipos corretos
    'block_medium': 4,
    'block_large': 8,
    'block_small': 2,
    'block_size': 2,
    'block_strategy': None,
    'local_iters': 3,
}
