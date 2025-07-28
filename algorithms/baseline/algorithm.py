"""
Greedy consensus algorithm (Baseline) for the Closest String Problem.

Classes:
    BaselineAlg: Baseline algorithm implementation.
"""

from src.domain.algorithms import CSPAlgorithm, register_algorithm

from .implementation import greedy_consensus, max_distance


@register_algorithm
class BaselineAlg(CSPAlgorithm):
    """
    Greedy consensus algorithm (Baseline) for the Closest String Problem.

    Args:
        strings (list[str]): List of input strings.
        alphabet (str): Used alphabet.

    Methods:
        run(): Execute the algorithm and return (center, max_distance, metadata).
    """

    name = "Baseline"
    default_params: dict = {}
    is_deterministic = True
    supports_internal_parallel = False  # Baseline doesn't support internal parallelism

    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Initialize the Baseline algorithm.

        Args:
            strings: List of dataset strings
            alphabet: Used alphabet
            **params: Algorithm-specific parameters
        """
        super().__init__(strings, alphabet, **params)

    def run(self) -> tuple[str, int, dict]:
        """
        Execute greedy consensus and calculate maximum Hamming distance.

        Returns:
            tuple[str, int, dict]: (center string, max distance, metadata)
        """
        import time

        start_time = time.time()

        self._report_progress("Starting greedy consensus...")

        # Save initial state in history
        self._save_history_entry(
            0,
            phase="initialization",
            total_strings=len(self.strings),
            string_length=len(self.strings[0]) if self.strings else 0,
        )

        center = greedy_consensus(self.strings, self.alphabet)

        self._report_progress("Computing maximum distance...")
        dist = max_distance(center, self.strings)

        end_time = time.time()
        execution_time = end_time - start_time

        # Save final state to history
        self._save_history_entry(
            1,
            phase="completion",
            center_found=center,
            max_distance=dist,
            execution_time=execution_time,
        )

        metadata = {
            "iterations": 1,  # Deterministic algorithm executed once
            "center_found": center,
            "total_strings": len(self.strings),
            "execution_time": execution_time,
            "history": self.get_history() if self.save_history else [],
        }

        return center, dist, metadata
