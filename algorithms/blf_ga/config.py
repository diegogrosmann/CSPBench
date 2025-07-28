"""
Default configurations for BLF-GA algorithm.

BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) is a hybrid metaheuristic
that combines local block learning with global genetic algorithm.

MAIN PARAMETERS:
- Population: Controls search size and quality
- Blocks: Defines local learning granularity
- Operators: Configures genetic evolution
- Adaptation: Dynamic adjustment mechanisms
- Criteria: Controls when to stop/restart

USAGE RECOMMENDATIONS:
- Small instances (n<10): pop_size=50-100, max_gens=50-100
- Medium instances (n=10-50): pop_size=100-200, max_gens=100-200
- Large instances (n>50): pop_size=200-500, max_gens=200-500

Attributes:
    BLF_GA_DEFAULTS (dict): Default BLF-GA parameters.
"""

# BLF-GA Configuration
BLF_GA_DEFAULTS = {
    # --- POPULATION PARAMETERS ---
    "pop_size": 1.5,  # Population size (fixed int or float >=1 for n multiplier)
    # Examples: 100 (fixed), 1.5 (1.5*n), 2.0 (2*n)
    "min_pop_size": 20,  # Minimum population size when using n proportion
    # Prevents very small populations in smaller instances
    "seed": None,  # Seed for reproducibility (None = random)
    # --- BLOCK PARAMETERS (LEARNING) ---
    "initial_blocks": 0.2,  # Number of initial blocks (fixed int or float 0-1 for L proportion)
    # Examples: 5 (fixed), 0.2 (20% of length)
    "min_block_len": 1,  # Minimum block size (prevents very small blocks)
    "rediv_freq": 10,  # Redivision frequency (every X generations)
    # Allows dynamic adaptation of block structure
    # --- GENETIC OPERATORS ---
    "cross_prob": 0.9,  # Crossover probability (0.6 to 0.95)
    # Higher value = more recombination
    "crossover_type": "one_point",  # Crossover type: one_point, uniform, blend_blocks
    "mut_prob": 0.1,  # Mutation probability (0.01 to 0.2)
    # Balances exploration vs exploitation
    "mutation_type": "multi",  # Mutation type: multi, inversion, transposition
    "mutation_multi_n": 2,  # Number of positions for multi mutation
    "elite_rate": 0.05,  # Elite rate (0.01 to 0.1) - % of best preserved
    "tournament_k": 2,  # Tournament size for selection
    # --- DIVERSITY AND IMMIGRANTS ---
    "immigrant_freq": 10,  # Generate immigrants every X generations
    # Injects diversity by replacing worst with random
    "immigrant_ratio": 0.2,  # Immigrant proportion (% of population)
    "diversity_threshold": 0.4,  # Diversity threshold (0-1)
    # Below this value, activates adaptive mutation
    # --- ADAPTIVE MUTATION ---
    "mutation_adapt_N": 10,  # N generations to detect convergence
    "mutation_adapt_factor": 2.0,  # Temporary mutation increase factor
    "mutation_adapt_duration": 5,  # Duration of mutation increase (generations)
    # --- NICHING ---
    "niching": False,  # Activate niching (preserves local diversity)
    "niching_radius": 3,  # Niche radius (minimum distance between solutions)
    # --- LOCAL REFINEMENT ---
    "refinement_type": "greedy",  # Type: greedy, swap, insertion, 2opt
    "refine_elites": "best",  # Who to refine: all (all elites), best (only best)
    "refine_iter_limit": 100,  # Iteration limit per refinement
    # --- STOPPING AND RESTART CRITERIA ---
    "max_gens": 100,  # Maximum number of generations (30 to 100+)
    "max_time": 1200.0,  # Maximum time in seconds (20 minutes)
    "no_improve_patience": 0.2,  # Generations without improvement for early stopping
    # (fixed int or float 0-1 for proportion of max_gens)
    "restart_patience": 20,  # Generations without improvement for partial restart
    "restart_ratio": 0.3,  # Proportion of population to restart
    "disable_elitism_gens": 5,  # Disable elitism every X generations
    # Prevents premature convergence
}
