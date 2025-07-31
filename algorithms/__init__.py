"""
CSP algorithms package.

This module initializes the algorithms package and implements the
auto-discovery system for algorithms. All implemented algorithms are
automatically registered in the global registry.

IMPLEMENTED ALGORITHMS:

1. **Baseline (Simple Algorithms)**:
   - Greedy Consensus: Consensus by majority voting
   - Baselines for performance comparison

2. **BLF-GA (Blockwise Learning Fusion + Genetic Algorithm)**:
   - Advanced hybrid metaheuristic
   - Combines block learning with genetic evolution
   - Adaptive mechanisms for diversity and convergence
   - Efficient parallelization and advanced configurability

3. **CSC (Consensus String Clustering)**:
   - Divide and conquer strategy
   - Hierarchical clustering of strings
   - Local optimization with global recombination

4. **DP-CSP (Dynamic Programming CSP)**:
   - Dynamic programming with pruning
   - Optimality guarantee within resource limits
   - Efficient state modeling

5. **H³-CSP (Hybrid Hierarchical Hamming Search)**:
   - Three-layer hierarchical approach
   - Adaptive technique selection per block
   - Automatic balancing between quality and efficiency

AUTO-DISCOVERY:
The system automatically discovers and registers all algorithms
implemented in subpackages, enabling dynamic usage through the
global registry.

USAGE EXAMPLE:
```python
from cspbench.domain.algorithms import global_registry

# List available algorithms
print("Available algorithms:")
for name, cls in global_registry.items():
    print(f"  {name}: {cls.__doc__.split('.')[0]}")

# Use a specific algorithm
algorithm_class = global_registry["Baseline"]
algorithm = algorithm_class(strings=["ATCG", "ATCC"], alphabet="ATCG")
result_string, max_distance, metadata = algorithm.run()
```

STRUCTURE:
Each algorithm must be in its own subpackage with:
- __init__.py: Main class exposure
- algorithm.py: Wrapper with @register_algorithm decorator
- implementation.py: Algorithm-specific logic
- config.py: Configurations and default parameters
- README.md: Detailed documentation
"""

import importlib
import pkgutil
from pathlib import Path

# Import registry from domain
from src.domain.algorithms import global_registry, register_algorithm


# Auto-discovery and algorithm import
def _discover_algorithms():
    """Discover and automatically import all algorithms."""
    algorithms_path = Path(__file__).parent

    for importer, modname, ispkg in pkgutil.iter_modules([str(algorithms_path)]):
        if ispkg and not modname.startswith("_"):
            try:
                # Import subpackage to activate automatic registration
                importlib.import_module(f"algorithms.{modname}")
            except ImportError as e:
                # Algorithm may have optional dependencies
                print(f"Warning: Algorithm '{modname}' could not be loaded: {e}")


# Execute auto-discovery on import
_discover_algorithms()

__all__ = ["global_registry", "register_algorithm"]

__all__ = ["global_registry", "register_algorithm"]
