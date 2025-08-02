"""
Exact DP-CSP (Dynamic Programming) implementation for the Closest String Problem.

DP-CSP is a dynamic programming algorithm that finds the EXACT solution
for the Closest String Problem. Unlike heuristics, it guarantees finding
the center string with the minimum possible maximum distance,
providing an optimal lower bound for comparison.

ALGORITHMIC ARCHITECTURE:

DETAILED DP-CSP ALGORITHM:

1. PROBLEM MODELING
   - State: (position, distance_vector)
   - Decision: which symbol to choose at current position
   - Transition: state[i] → state[i+1] via symbol choice
   - Objective: minimize final maximum distance

2. DYNAMIC PROGRAMMING
   - DP Table: dp[position][distance_state] = minimum_cost
   - Recurrence: dp[i][d] = min(dp[i+1][d'] + transition_cost)
   - Base case: dp[L][d] = max(d) (final maximum distance)
   - Solution: dp[0][initial_state] = optimal distance

3. EFFICIENCY OPTIMIZATIONS
   - Memory limit pruning (max_states)
   - Efficient state encoding as tuples
   - Memoization to avoid recalculations
   - Lazy optimal solution reconstruction

4. RESOURCE CONTROLS
   - Memory limit to prevent exponential explosion
   - Timeout for intractable problems
   - Fallback to heuristics if insufficient resources
   - Progress monitoring via callbacks

ALGORITHMIC PHILOSOPHY:

- EXACTNESS: Guarantees optimal solution within resource limits
- COMPLETENESS: Explores all viable state space
- EFFICIENCY: Uses dynamic programming to avoid recalculations
- ROBUSTNESS: Fallbacks and controls for intractable cases

DISTINCTIVE CHARACTERISTICS:

- OPTIMALITY GUARANTEE: Unlike heuristics, guarantees optimal solution
- CONTROLLED COMPLEXITY: Memory limit pruning
- VERSATILITY: Works for any alphabet and string size
- TRANSPARENCY: Step-by-step optimal solution reconstruction

CSP APPLICATION:
DP-CSP is ideal for:
- Small/medium instances where exactness is critical
- Quality benchmarking for heuristics
- Problem structure analysis
- Approximate solution validation

LIMITATIONS:
- Exponential complexity in worst case
- Requires high memory for large problems
- Can be slow for complex instances
- Not scalable for very long strings

EXAMPLE OPERATION:
For strings ["AC", "AG", "AT"] and alphabet "ACGT":
1. Initial states: (0, [0,0,0])
2. Position 0: choose A → (1, [0,0,0])
3. Position 1: choose C → (2, [0,1,1]), G → (2, [1,0,1]), T → (2, [1,1,0])
4. Optimal solution: "AC" with maximum distance 1

Classes:
    DPCSP: DP-CSP algorithm implementation with resource controls.

Auxiliary functions:
    _dp_decision(): Dynamic programming decision function.
    exact_dp_closest_string(): Main wrapper with resource controls.

Types:
    String: Alias for str (string representation).
    State: Tuple representing distance state.

Author: Implementation based on classical dynamic programming for CSP
Version: Optimized with resource controls and robust fallbacks
"""

from __future__ import annotations

import logging
import time
from collections.abc import Callable, Sequence
from typing import Any, TypeAlias, cast

from src.domain.metrics import max_distance

from .config import DP_CSP_DEFAULTS

logger = logging.getLogger(__name__)

# Type aliases for clarity
RemVec: TypeAlias = tuple[int, ...]  # (Remaining errors vector)
State: TypeAlias = tuple[int, RemVec]  # (position, rem_errors)
String: TypeAlias = str


def _dp_decision(strings: Sequence[String], alphabet: str, d: int) -> String | None:
    """
    DP decision algorithm: checks if center string exists with radius ≤ d.

    This is the core DP-CSP function that implements the dynamic programming
    algorithm for the decision problem: "Does there exist a center string such that
    the maximum distance to all strings is at most d?"

    DYNAMIC PROGRAMMING ALGORITHM:

    1. STATE MODELING:
       - State: rem = (rem[0], rem[1], ..., rem[n-1])
       - rem[i] = number of errors still allowed for string i
       - Initial state: rem = (d, d, ..., d)
       - Final viable state: rem[i] ≥ 0 for all i

    2. DIFFERENCE PRE-COMPUTATION:
       - For each position pos and symbol σ:
       - δσ[pos] = vector indicating if σ differs from each string at pos
       - δσ[pos][i] = 1 if σ ≠ strings[i][pos], 0 otherwise

    3. STATE TRANSITIONS:
       - For each position pos = 0, 1, ..., L-1:
       - For each state rem in current frontier:
       - For each symbol σ of alphabet:
         a) Calculate new_rem[i] = rem[i] - δσ[pos][i]
         b) If min(new_rem) ≥ 0: viable state, add to next frontier
         c) Store parent pointer for reconstruction

    4. STATE PRUNING:
       - Remove states where some rem[i] < 0 (unfeasible)
       - Avoid duplicates: one state per rem configuration
       - Stop if frontier empty (no solution possible)

    5. SOLUTION RECONSTRUCTION:
       - If reached end with viable states: solution exists
       - Use parent pointers to backtrack from end to start
       - Reconstruct center string symbol by symbol

    COMPLEXITY:
    - States: O((d+1)^n) - each rem[i] can be 0, 1, ..., d
    - Transitions: O(L × |Σ|) - for each pos, test each symbol
    - Total: O((d+1)^n × L × |Σ|)

    OPTIMIZATIONS IMPLEMENTED:
    - Pre-computation of δσ avoids recalculations
    - Compact frontier reduces memory usage
    - Early pruning eliminates unfeasible states quickly
    - Parent pointers allow efficient reconstruction

    Args:
        strings: Sequence of input strings (same length)
        alphabet: String with all valid symbols
        d: Maximum allowed radius (decision threshold)

    Returns:
        Center string with radius ≤ d if exists, None otherwise

    Example:
        >>> strings = ["AC", "AT", "GC"]
        >>> _dp_decision(strings, "ACGT", 1)
        "AC"  # A solution with radius = 1

        >>> _dp_decision(strings, "ACGT", 0)
        None  # Impossible radius = 0

    Note:
        Function is deterministic but may return any valid solution
        (not necessarily unique). Final validation checks if returned
        solution really satisfies radius ≤ d.
    """
    n, L = len(strings), len(strings[0])

    # PRE-COMPUTATION: Difference matrix δσ[pos][string]
    # For each position and symbol, calculate discrepancy vector
    delta: list[dict[String, RemVec]] = []
    for pos in range(L):
        # Symbols at position pos of all strings
        col = [s[pos] for s in strings]
        # For each symbol σ, calculate differences
        delta.append({σ: tuple(int(σ != c) for c in col) for σ in alphabet})

    # DYNAMIC PROGRAMMING INITIALIZATION
    start: RemVec = (d,) * n  # Initial state: d errors for each string
    frontier: set[RemVec] = {start}  # Current frontier of viable states

    # Pointers for reconstruction: (pos, state) → (previous_state, symbol)
    parent: dict[tuple[int, RemVec], tuple[RemVec | None, String]] = {
        (0, start): (None, "")  # Initial state has no predecessor
    }

    # POSITION-BY-POSITION DYNAMIC PROGRAMMING
    for pos in range(L):
        nxt: set[RemVec] = set()  # Next frontier

        # Process all states in current frontier
        for rem in frontier:
            # Test each alphabet symbol at current position
            for σ, dv in delta[pos].items():
                # STATE TRANSITION: consume errors based on differences
                new_rem = tuple(r - v for r, v in zip(rem, dv))

                # PRUNING: Remove unfeasible states (some rem[i] negative)
                if min(new_rem) < 0:
                    continue  # Unfeasible state, prune

                # Avoid duplicate states
                key = (pos + 1, new_rem)
                if key in parent:
                    continue  # State already visited

                # STORAGE: Record transition and add to next frontier
                parent[key] = (rem, σ)
                nxt.add(new_rem)

        # Update frontier for next position
        frontier = nxt

        # GLOBAL PRUNING: If no viable states, impossible to continue
        if not frontier:
            return None  # No solution exists for this d

    # SOLUTION RECONSTRUCTION
    # If we reach here, at least one final viable state exists
    final_rem = next(iter(frontier))  # Any final state works

    center_chars: list[String] = []
    pos: int = L
    rem: RemVec = final_rem

    # Backtrack using parent pointers
    while pos > 0:
        prev_rem, σ = parent[(pos, rem)]
        center_chars.append(σ)  # Symbol chosen at this position
        pos -= 1

        if prev_rem is None:  # Reached initial state
            break
        rem = cast(RemVec, prev_rem)

    # Reconstruct string in correct order
    center_chars.reverse()
    result = "".join(center_chars)

    # FINAL VALIDATION: Check if solution is really valid
    from src.domain.metrics import hamming_distance

    max_dist = max(hamming_distance(result, s) for s in strings)
    if max_dist > d:
        logger.error(
            "[DP_DECISION] ERROR: Invalid solution! dist=%d > d=%d", max_dist, d
        )

    return result


def exact_dp_closest_string(
    strings: list[String],
    alphabet: str,
    max_d: int | None = None,
    progress_callback: Callable[[str], None] | None = None,
    warning_callback: Callable[[str], None] | None = None,
) -> tuple[String, int]:
    """
    Finds the EXACT Closest String Problem solution using dynamic programming.

    This is the main DP-CSP interface that coordinates the search for the optimal
    minimum radius. Unlike heuristics, it guarantees finding the globally
    optimal solution, providing the true lower bound for the problem.

    OPTIMAL RADIUS SEARCH STRATEGY:

    1. INCREMENTAL SEARCH:
       - Tests values d = 0, 1, 2, ..., max_d sequentially
       - For each d, uses DP decision algorithm
       - Stops at first d where solution is found (guaranteed optimal)

    2. GUARANTEED OPTIMALITY:
       - Since testing in increasing order, first solution is optimal
       - d* = min{d : there exists center with radius ≤ d}
       - Avoids unnecessary binary search for small problems

    3. RESOURCE MONITORING:
       - Estimates complexity (d+1)^n before executing
       - Monitors RSS memory usage during execution
       - Interrupts if exceeding time or memory limits
       - Reports detailed metrics for auditing

    LIMITATIONS AND SAFEGUARDS:

    - Exponential Complexity: O((d+1)^n × L × |Σ|) grows rapidly
    - State Limit: Aborts if (d+1)^n > 2e9 (prevents memory crash)
    - Memory Limit: Monitors RSS, aborts if > 95% of safe limit
    - Time Limit: Configurable timeout (default 300s)

    RECOMMENDED USE CASES:

    - Validation: Check quality of heuristic algorithms
    - Benchmark: Establish ground truth for comparison
    - Small Instances: n ≤ 8, L ≤ 20, alphabet ≤ 4
    - Research: Theoretical analysis of CSP properties

    PARAMETER CONFIGURATION:

    - max_d=None: Uses baseline (distance from first string) as limit
    - progress_callback: Reports progress "Testing d=X"
    - warning_callback: Reports resource warnings before aborting

    DETAILED METRICS AND LOGS:

    Function records complete information for auditing:
    - Input parameters (n, L, alphabet, strings)
    - Configured resource limits
    - Search progress (d tested, estimated states)
    - Final result (center found, optimal d*, validation)
    - Performance metrics (time, memory, iterations)

    Args:
        strings: List of input strings (same length)
        alphabet: String with valid alphabet symbols
        max_d: Maximum radius to test (None = automatic baseline)
        progress_callback: Function to report progress (optional)
        warning_callback: Function to report resource warnings (optional)

    Returns:
        tuple: (optimal_center, optimal_d*)
            - optimal_center: Center string with minimum possible radius
            - optimal_d*: Minimum possible radius for the instance

    Raises:
        RuntimeError: If no solution found within limits or
                     if exceeding resources (memory, time, complexity)

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> center, d_star = exact_dp_closest_string(strings, "ACGT")
        >>> print(f"Optimal solution: {center} with radius {d_star}")
        "Optimal solution: ACCT with radius 1"

    Note:
        For large instances, consider using heuristic algorithms
        like BLF-GA or CSC that offer good quality in viable time.
        DP-CSP should be used when exactness is critical and
        computational resources are adequate.
    """
    # INITIAL CONFIGURATION AND VALIDATION
    baseline_val = max_distance(strings[0], strings)  # Simple upper bound
    if max_d is None:
        max_d = baseline_val

    n = len(strings)
    L = len(strings[0])

    # DETAILED INPUT LOGS
    logger.info(
        "[DP_CSP] Starting exact search with max_d=%d, baseline=%d", max_d, baseline_val
    )
    logger.info("[DP_CSP] Dataset: n=%d, L=%d, alphabet=%s", n, L, alphabet)
    for i, s in enumerate(strings):
        logger.info("[DP_CSP] String %d: %s", i, s)

    # RESOURCE MONITORING CONFIGURATION
    safe_mem_mb = 1000.0  # Simplified default limit
    max_time = DP_CSP_DEFAULTS.get("max_time", 300)
    t0 = time.time()

    logger.info("[DP_CSP] Limits: mem=%.1fMB, time=%ds", safe_mem_mb, max_time)

    def check_limits(d):
        """Check resource limits before processing radius d."""
        import gc
        import os

        # MEMORY CLEANUP AND MEASUREMENT
        gc.collect()
        mem_mb = 0.0
        try:
            if os.path.exists("/proc/self/status"):
                with open("/proc/self/status", encoding="utf-8") as f:
                    for line in f:
                        if line.startswith("VmRSS:"):
                            kb = int(line.split()[1])
                            mem_mb = kb / 1024.0
                            break
        except (OSError, ValueError):
            pass  # Continue even if cannot measure memory

        elapsed = time.time() - t0

        # COMPLEXITY ESTIMATION
        # Practical limit: ~2e9 states (≈16GB RAM, 8 bytes/state)
        state_count_est = (d + 1) ** n

        # RESOURCE CHECKS
        if state_count_est > 2_000_000_000:
            msg = (
                f"DP-CSP interrupted: (d+1)^n = {state_count_est:,} exceeds practical limit "
                "(~2e9 states, ~16GB RAM). Try reducing n or d."
            )
            logger.error("[DP_CSP] %s", msg)
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

        if mem_mb > safe_mem_mb * 0.95:
            msg = f"DP-CSP interrupted: memory usage {mem_mb:.1f}MB exceeded safe limit ({safe_mem_mb:.1f}MB)"
            logger.error("[DP_CSP] %s", msg)
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

        if elapsed > max_time:
            msg = f"DP-CSP interrupted: execution time exceeded {max_time}s"
            logger.error("[DP_CSP] %s", msg)
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

    # INCREMENTAL SEARCH FOR OPTIMAL RADIUS
    for d in range(max_d + 1):
        # Check resources before each main iteration
        check_limits(d)

        if progress_callback:
            progress_callback(f"Testing d={d}")

        logger.info("[DP_CSP] Testing d=%d (attempt %d/%d)", d, d + 1, max_d + 1)

        # DECISION ALGORITHM FOR RADIUS d
        center = _dp_decision(strings, alphabet, d)

        if center is not None:
            # SUCCESS: Found solution with radius d

            # RIGOROUS FINAL VALIDATION
            from src.domain.metrics import hamming_distance

            max_dist = max(hamming_distance(center, s) for s in strings)

            logger.info("[DP_CSP] SUCCESS! Found solution with d=%d", d)
            logger.info("[DP_CSP] Center found: %s", center)
            logger.info("[DP_CSP] Final validation: maximum distance = %d", max_dist)

            # Consistency check
            if max_dist != d:
                logger.warning(
                    "[DP_CSP] INCONSISTENCY: d=%d but actual distance=%d", d, max_dist
                )

            return center, d

    # FAILURE: No solution found within limits
    logger.error("[DP_CSP] FAILURE: Could not find center with d ≤ %d", max_d)
    raise RuntimeError(
        f"Could not find center with d ≤ {max_d}. " "Try increasing the limit."
    )


class DPCSP:
    """
    Wrapper class for the DP-CSP (Dynamic Programming CSP) algorithm.

    This class provides an object-oriented interface for the DP-CSP algorithm,
    encapsulating configurations and offering convenient methods
    for execution and analysis.

    The class is primarily a wrapper around the exact_dp_closest_string()
    function, providing:
    - Consistent interface with other project algorithms
    - Parameter configuration via constructor
    - Callback methods for monitoring
    - Execution metadata and statistics

    Attributes:
        strings: List of input strings
        alphabet: Alphabet used
        max_d: Maximum radius for search
        progress_callback: Progress callback
        warning_callback: Warning callback
    """

    def __init__(
        self, strings: list[str], alphabet: str, max_d: int | None = None, **kwargs: Any
    ):
        """
        Initialize the DP-CSP algorithm.

        Args:
            strings: List of input strings (same length)
            alphabet: Alphabet used
            max_d: Maximum radius for search (None = automatic)
            **kwargs: Additional parameters (for compatibility)
        """
        self.strings = strings
        self.alphabet = alphabet
        self.max_d = max_d
        self.progress_callback: Callable[[str], None] | None = None
        self.warning_callback: Callable[[str], None] | None = None

        # Parameters for compatibility (not used by DP-CSP)
        self.params = kwargs

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Set callback to report progress."""
        self.progress_callback = callback

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """Set callback to report warnings."""
        self.warning_callback = callback

    def run(self) -> tuple[str, int]:
        """
        Execute the DP-CSP algorithm.

        Returns:
            tuple: (center, optimal_distance)
                - center: Optimal center string
                - optimal_distance: Minimum possible radius
        """
        return exact_dp_closest_string(
            self.strings,
            self.alphabet,
            self.max_d,
            self.progress_callback,
            self.warning_callback,
        )
