"""
Default configurations for the H³-CSP algorithm.

This module defines all default parameters for the H³-CSP algorithm,
including block division parameters, search techniques,
refinement and execution control.

Constants:
    H3_CSP_DEFAULTS (dict): Dictionary with all default parameters
                           for the H³-CSP algorithm.
"""

# H³-CSP Configuration
H3_CSP_DEFAULTS = {
    # === Block Division Parameters ===
    "auto_blocks": True,  # Use automatic division by √L
    "min_block_size": 2,  # Minimum block size
    "max_blocks": None,  # Maximum blocks (None = automatic)
    "block_size": 2,  # Base block size (if not automatic)
    "block_strategy": None,  # Division strategy (None = default)
    # === Difficulty Thresholds per Block ===
    "block_small": 2,  # Limit for "small" blocks (exhaustive search)
    "block_medium": 4,  # Limit for "medium" blocks (reduced beam search)
    "block_large": 8,  # Limit for "large" blocks (full beam search)
    # === Search Parameters ===
    "exhaustive_limit": 10000,  # Limit for exhaustive search (|Σ|^m)
    "beam_width": 32,  # Beam search width
    "k_candidates": 5,  # Number of candidates per block
    # === Global Refinement ===
    "local_search_iters": 3,  # Local search iterations (deprecated)
    "local_iters": 3,  # Local search iterations (current)
    # === Execution Control ===
    "max_time": 300,  # Maximum time in seconds
    "seed": None,  # Seed for reproducibility
    # === Experimental Parameters ===
    "diversity_threshold": 1,  # Diversity threshold (not currently used)
    "fallback_enabled": True,  # Enable fallback for large blocks
}

# Parameter mapping for compatibility
# Some parameters have alternative names for compatibility with previous versions
_PARAMETER_ALIASES = {
    "local_search_iters": "local_iters",  # Alias for compatibility
}
