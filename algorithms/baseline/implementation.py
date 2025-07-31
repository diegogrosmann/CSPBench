"""
Implementation of the greedy consensus algorithm (Baseline) for CSP.

The Baseline algorithm represents the simplest and most fundamental approach to the
Closest String Problem (CSP). It uses a greedy strategy that constructs
the center string position by position, always choosing the symbol that locally
minimizes the maximum distance at that moment.

GREEDY CONSENSUS ALGORITHM:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                        DETAILED BASELINE ALGORITHM                              │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. INITIALIZATION                                                              │
│   ├── Validate input (non-empty strings, same length)                          │
│   └── Initialize empty consensus string                                         │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. POSITION-BY-POSITION CONSTRUCTION                                          │
│   ├── For each position i from 0 to L-1:                                       │
│   │   ├── For each symbol c in alphabet:                                       │
│   │   │   ├── Calculate partial string = consensus[0:i] + c                    │
│   │   │   ├── Calculate maximum distance of partial string                     │
│   │   │   └── Store if it's the best option so far                            │
│   │   └── Choose symbol that minimizes maximum distance                        │
│   └── Add best symbol to consensus                                             │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 3. VALIDATION AND RETURN                                                      │
│   ├── Calculate final distance of complete consensus string                     │
│   ├── Log result for auditing                                                  │
│   └── Return constructed consensus string                                       │
└─────────────────────────────────────────────────────────────────────────────────┘

ALGORITHM CHARACTERISTICS:

• **SIMPLICITY**: Direct and comprehensible implementation
• **DETERMINISTIC**: Always produces the same result for the same input
• **LOCAL GREEDY**: Optimizes each position independently
• **RELIABLE BASELINE**: Serves as reference for comparison
• **LOW COMPLEXITY**: O(L * |Σ| * n) where L=length, |Σ|=alphabet, n=strings

LIMITATIONS:

• **LOCAL OPTIMUM**: May not find global optimal solution
• **MYOPIA**: Local decisions may harm global quality
• **NO BACKTRACKING**: Does not reconsider previous decisions
• **LIMITED QUALITY**: Generally inferior to advanced metaheuristics

APPLICATION AND USE:

The Baseline is fundamental as:
- **Comparison reference** for other algorithms
- **Initial solution** for more sophisticated methods
- **Implementation validation** (must work correctly)
- **Problem difficulty analysis** (if Baseline solves well, problem is easy)

OPERATION EXAMPLE:

Input: ["ACGT", "AGCT", "ATCT"]
Position 0: A=3/3, best='A' (partial distance=0)
Position 1: C vs G vs T → C minimizes maximum distance
Position 2: G vs C vs T → C minimizes maximum distance
Position 3: T=3/3, best='T' (minimal partial distance)
Result: "ACCT" with maximum distance = 1

Classes:
    None - Direct functional implementation

Main functions:
    greedy_consensus(strings, alphabet): Main greedy consensus algorithm
    max_distance(center, strings): Auxiliary function for distance calculation

Author: Reference implementation based on classic greedy algorithms
Version: Optimized for clarity and to serve as reliable baseline
"""

import logging

logger = logging.getLogger(__name__)


def greedy_consensus(strings: list[str], alphabet: str) -> str:
    """
    Constructs a consensus string using greedy strategy position by position.

    This is the Baseline algorithm implementation for the Closest String Problem.
    The greedy strategy constructs the solution incrementally, making
    the locally optimal decision at each step without considering global impact.

    DETAILED ALGORITHM:

    1. **INPUT VALIDATION**:
       - Checks if there are strings to process
       - Assumes all have the same length (precondition)

    2. **INCREMENTAL CONSTRUCTION**:
       - For each position i from 0 to L-1:
         a) Tests each alphabet symbol at this position
         b) For each candidate symbol:
            - Constructs partial string up to position i
            - Calculates maximum distance of partial string
         c) Chooses symbol that minimizes maximum distance
         d) Adds chosen symbol to consensus

    3. **LOCAL OPTIMIZATION STRATEGY**:
       - At each position, chooses symbol that results in smallest
         maximum distance considering only the constructed prefix
       - Does not consider future impact (greedy characteristic)

    TEMPORAL COMPLEXITY:
    - O(L × |Σ| × n × L) where:
      - L = string length
      - |Σ| = alphabet size
      - n = number of strings
      - Additional L factor comes from partial distance calculation

    CHARACTERISTICS:
    - **Deterministic**: Always produces the same result
    - **Greedy**: Local optimization without backtracking
    - **Incremental**: Constructs solution position by position
    - **Baseline**: Reference for comparison with other algorithms

    LIMITATIONS:
    - May get stuck in local optima
    - Does not guarantee globally optimal solution
    - Quality depends on construction order

    Args:
        strings: List of input strings (all with same length)
        alphabet: String containing all valid alphabet symbols

    Returns:
        str: Consensus string constructed by greedy strategy

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> alphabet = "ACGT"
        >>> greedy_consensus(strings, alphabet)
        "ACCT"  # Maximum distance = 1

    Note:
        The function logs the final result for auditing
        and implementation validation.
    """
    # INPUT VALIDATION
    if not strings:
        return ""

    # INITIALIZATION
    L = len(strings[0])  # String length (assumed uniform)
    consensus = []  # Consensus string being constructed

    # POSITION-BY-POSITION CONSTRUCTION
    for pos in range(L):
        best_char = None  # Best symbol for current position
        best_max_dist = float("inf")  # Smallest maximum distance found

        # TEST EACH ALPHABET SYMBOL
        for char in alphabet:
            # Build partial string with candidate symbol
            partial_consensus = consensus + [char]

            # PARTIAL MAXIMUM DISTANCE CALCULATION
            # Calculates distance only considering already decided positions
            max_dist = 0
            for s in strings:
                # Partial Hamming distance (up to current position + 1)
                dist = sum(1 for i in range(pos + 1) if partial_consensus[i] != s[i])
                max_dist = max(max_dist, dist)

            # BEST CANDIDATE UPDATE
            # If found symbol that results in smaller maximum distance
            if max_dist < best_max_dist:
                best_max_dist = max_dist
                best_char = char

        # GREEDY DECISION: Add best symbol found
        consensus.append(best_char)

    # FINAL STRING CONSTRUCTION
    result = "".join(consensus)

    # RESULT VALIDATION AND LOGGING
    # Local import to avoid name conflict
    from src.domain.metrics import max_distance as calc_max_distance

    final_distance = calc_max_distance(result, strings)
    logger.info("[CONSENSUS] Consensus: %s, distance: %d", result, final_distance)

    return result


def max_distance(center: str, strings: list[str]) -> int:
    """
    Calculates the maximum Hamming distance between a center string and a set of strings.

    This function implements the objective function of the Closest String Problem (CSP).
    The maximum distance is the value we need to minimize to find
    the best possible center string.

    MATHEMATICAL DEFINITION:
    For a center string c and set of strings S = {s₁, s₂, ..., sₙ}:
    max_distance(c, S) = max{d_H(c, sᵢ) | sᵢ ∈ S}

    where d_H(a,b) is the Hamming distance between strings a and b.

    IMPORTANCE IN CSP:
    - **Objective Function**: Value to be minimized
    - **Quality Criterion**: Lower value = better solution
    - **Fitness**: Used to evaluate candidates in evolutionary algorithms
    - **Stopping Criterion**: Algorithm stops when reaches 0 (optimal)

    PROPERTIES:
    - **Symmetric**: max_distance(c, S) = max_distance(c, S')
    - **Bounded**: 0 ≤ result ≤ L (string length)
    - **Monotonic**: More differences → greater distance
    - **Discrete**: Only integer values

    COMPLEXITY:
    - Temporal: O(n × L) where n=number of strings, L=length
    - Spatial: O(1) - requires no additional storage

    SPECIAL CASES:
    - If center ∈ strings, then result ≤ set diameter
    - If strings are identical and center = string, result = 0
    - Worst case: center differs from all at all positions

    Args:
        center: Candidate center string (same length as strings)
        strings: List of input strings for comparison

    Returns:
        int: Largest Hamming distance found between center and strings

    Example:
        >>> center = "ACCT"
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> max_distance(center, strings)
        1  # max(1, 1, 1) = 1

    Note:
        This function is a wrapper to maintain consistency with baseline
        module interface, delegating calculation to optimized implementation.
    """
    # Local import to use optimized implementation
    from src.domain.metrics import hamming_distance

    # CALCULATE ALL DISTANCES
    # Calculates Hamming distance between center and each string
    distances = [hamming_distance(center, s) for s in strings]

    # RETURN MAXIMUM DISTANCE
    # max() finds the largest value in the distance list
    return max(distances)
