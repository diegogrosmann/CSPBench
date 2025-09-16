"""Core DP decision function for DP-CSP refactored algorithm.

Mantemos apenas a função `_dp_decision` (problema de decisão) utilizada pelo
`DPCSPAlgorithm`. O controle de busca por d, progressos e tratamento
de complexidade agora vivem em `algorithm.py`.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import TypeAlias, cast

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
    # Cálculo inline de distância de Hamming para validação final
    def _hamming(a: str, b: str) -> int:
        return sum(c1 != c2 for c1, c2 in zip(a, b))

    max_dist = max(_hamming(result, s) for s in strings)
    if max_dist > d:
        logger.error(
            "[DP_DECISION] ERROR: Invalid solution! dist=%d > d=%d", max_dist, d
        )

    return result


# Removed legacy functions: exact_dp_closest_string / DPCSP wrapper.
