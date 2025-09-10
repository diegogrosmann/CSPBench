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

This module now contains both the Baseline algorithm class and its greedy
consensus functional implementation (previously in implementation.py) to
simplify structure as requested. If external code/tests previously imported
`greedy_consensus` from `algorithms.baseline.implementation`, update imports
to use `algorithms.baseline.algorithm` (or simply `algorithms.baseline` which
re-exports it via `__init__`).

Classes:
    BaselineAlg: Baseline algorithm implementation.
Functions:
    greedy_consensus: Greedy construction of consensus string.
"""

from src.domain.algorithms import AlgorithmResult, CSPAlgorithm, register_algorithm

from .config import BASELINE_DEFAULTS


@register_algorithm
class BaselineAlg(CSPAlgorithm):
    """
    Greedy consensus algorithm (Baseline) for the Closest String Problem.

    The Baseline algorithm uses a greedy strategy that constructs the center
    string by position by position, always choosing the symbol that locally
    minimizes the maximum distance at that moment.

    This algorithm serves as a fundamental reference for comparison with
    more sophisticated algorithms.

    Args:
        strings: List of input strings
        alphabet: Used alphabet
        distance_calculator: Distance calculator instance
        store: WorkStatePersistence for progress/warnings reporting
        global_seed: Global seed that overrides local seeds
        internal_jobs: Number of internal parallel jobs (default: 1)
        **params: Algorithm-specific parameters

    Methods:
        run(): Execute the algorithm and return AlgorithmResult.
    """

    name = "Baseline"
    default_params: dict = BASELINE_DEFAULTS
    supports_internal_parallel = False

    def run(self) -> AlgorithmResult:
        """
        Execute the Baseline greedy consensus algorithm.

        Returns:
            AlgorithmResult: Dictionary containing center string, max distance,
                           parameters, execution metadata, success status and error info
        """
        import time

        start_time = time.time()

        try:
            # Report algorithm start
            if self._monitor:
                self._monitor.on_progress(0.0, "Starting Baseline algorithm")

            # Verificação inicial de cancelamento
            if self._monitor and self._monitor.is_cancelled():
                return self._build_cancelled_result(start_time)

            # Validate input
            if not self.strings:
                raise ValueError("String list cannot be empty")

            if not self.alphabet:
                raise ValueError("Alphabet cannot be empty")

            self.tie_break = self.params.get("tie_break", "lex")
            # Execute greedy consensus algorithm
            center_string = self.greedy_consensus()

            # Verificar cancelamento após construção
            if self._monitor and self._monitor.is_cancelled():
                return self._build_cancelled_result(start_time)

            # Calculate final metrics
            max_distance = self.max_distance(center_string)
            avg_distance = self.average_distance(center_string)
            total_distance = self.total_distance(center_string)

            end_time = time.time()
            execution_time = end_time - start_time

            # Report completion
            if self._monitor:
                self._monitor.on_progress(1.0, "Baseline algorithm finished")

            # Build successful result
            return AlgorithmResult(
                success=True,
                center_string=center_string,
                max_distance=max_distance,
                parameters=self.get_actual_params(),
                error=None,
                metadata={
                    "algorithm_name": self.name,
                    "execution_time": execution_time,
                    "avg_distance": avg_distance,
                    "total_distance": total_distance,
                    "num_strings": len(self.strings),
                    "string_length": len(self.strings[0]) if self.strings else 0,
                    "alphabet_size": len(self.alphabet),
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                },
            )

        except Exception as e:
            # Calculate execution time even on error
            end_time = time.time()
            execution_time = end_time - start_time

            error_message = f"Error executing Baseline algorithm: {str(e)}"

            # Report error
            if self._monitor:
                self._monitor.on_progress(0.0, f"Error: {error_message}")
                self._monitor.on_warning(error_message)

            # Build error result
            return AlgorithmResult(
                success=False,
                center_string="",  # Empty string on error
                max_distance=-1,  # Invalid distance to indicate error
                parameters=self.get_actual_params(),
                error=error_message,
                metadata={
                    "algorithm_name": self.name,
                    "execution_time": execution_time,
                    "avg_distance": -1,
                    "total_distance": -1,
                    "num_strings": len(self.strings) if self.strings else 0,
                    "string_length": (
                        len(self.strings[0]) if self.strings and self.strings[0] else 0
                    ),
                    "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                    "error_type": type(e).__name__,
                },
            )

    def greedy_consensus(self) -> str:
        """Construct a consensus string using a greedy position-by-position strategy.

        This is the Baseline algorithm implementation for the Closest String Problem.
        The greedy strategy constructs the solution incrementally, making the locally
        optimal decision at each step without considering global impact.

        Args:
            strings: List of input strings (all with same length)
            alphabet: String containing all valid alphabet symbols
            algorithm: Optional CSPAlgorithm instance for distance calculation / progress
            tie_break: Tie-breaking strategy ('lex', 'random', 'first')

        Returns:
            Consensus string constructed by greedy strategy (empty if no strings)
        """
        if not self.strings:
            return ""

        L = len(self.strings[0])  # Assumed uniform length (validated elsewhere)
        consensus: list[str] = []

        if self._monitor:
            self._monitor.on_progress(
                0.05, f"Building consensus string position by position (L={L})"
            )

        for pos in range(L):
            # Verificar cancelamento durante construção posição-por-posição
            if self._monitor and self._monitor.is_cancelled():
                # Retornar consenso parcial se cancelado
                return "".join(consensus)

            # Frequency counting approach (majority vote with tie-break)
            counts: dict[str, int] = dict.fromkeys(self.alphabet, 0)
            for s in self.strings:
                counts[s[pos]] = counts.get(s[pos], 0) + 1

            max_count = max(counts[c] for c in self.alphabet)
            candidates = [
                c for c in self.alphabet if counts[c] == max_count and max_count > 0
            ]

            # If all counts are zero (pathological – shouldn't happen), fallback to lex order entire alphabet
            if not candidates:  # pragma: no cover
                candidates = list(self.alphabet)

            if len(candidates) == 1:
                best_char = candidates[0]
            else:
                if self.tie_break == "lex":
                    best_char = min(candidates)
                elif self.tie_break == "random":
                    best_char = self.rng.choice(candidates)
                elif self.tie_break == "first":
                    best_char = candidates[0]
                else:  # fallback
                    best_char = min(candidates)

            consensus.append(best_char)

            progress = 0.05 + (pos + 1) / L * 0.90
            if self._monitor:
                self._monitor.on_progress(
                    progress,
                    f"Position {pos + 1}/{L}: '{best_char}' (freq={max_count}, candidates={candidates})",
                )

        result = "".join(consensus)

        final_distance = self.max_distance(result)
        if self._monitor:
            self._monitor.on_progress(
                0.95,
                f"Consensus built: '{result}' (final max distance: {final_distance})",
            )
        return result

    def _build_cancelled_result(self, start_time: float) -> AlgorithmResult:
        """Constrói resultado para execução cancelada."""
        end_time = time.time()
        execution_time = end_time - start_time
        return AlgorithmResult(
            success=False,
            center_string="",
            max_distance=-1,
            parameters=self.get_actual_params(),
            error="Algorithm execution was cancelled",
            metadata={
                "algorithm_name": self.name,
                "execution_time": execution_time,
                "avg_distance": -1,
                "total_distance": -1,
                "num_strings": len(self.strings) if self.strings else 0,
                "string_length": (
                    len(self.strings[0]) if self.strings and self.strings[0] else 0
                ),
                "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                "seed": self.seed,
                "internal_jobs": self.internal_jobs,
                "termination_reason": "cancelled",
            },
        )
