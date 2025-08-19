"""
CSC: Consensus String Clustering for CSP.

Classes:
    CSCAlgorithm: CSC algorithm implementation.
"""

from src.domain.algorithms import CSPAlgorithm, register_algorithm
from src.domain.metrics import max_distance

from .config import CSC_DEFAULTS
from .implementation import heuristic_closest_string


@register_algorithm
class CSCAlgorithm(CSPAlgorithm):
    """
    CSC: Consensus String Clustering for the Closest String Problem.

    Args:
        strings (list[str]): List of input strings.
        alphabet (str): Alphabet used.
        **params: Algorithm parameters.

    Methods:
        run(): Executes CSC and returns (center, maximum distance, metadata).
    """

    name = "CSC"
    default_params = CSC_DEFAULTS
    supports_internal_parallel = False  # CSC doesn't support internal parallelism
    is_deterministic = True  # CSC is deterministic

    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Initialize CSC algorithm.

        Args:
            strings: List of dataset strings
            alphabet: Alphabet used
            **params: Algorithm-specific parameters

        Raises:
            ValueError: If strings have different lengths (CSC requirement)
        """
        # Validate that all strings have the same length (CSC requirement)
        if strings:
            string_lengths = [len(s) for s in strings]
            if len(set(string_lengths)) > 1:
                min_len, max_len = min(string_lengths), max(string_lengths)
                raise ValueError(
                    f"CSC requires all strings to have the same length. "
                    f"Found lengths ranging from {min_len} to {max_len}. "
                    f"Consider preprocessing the dataset to normalize string lengths."
                )

        super().__init__(strings, alphabet, **params)

    def run(self) -> tuple[str, int, dict]:
        """
        Execute CSC algorithm and return center string, maximum distance and metadata.

        Returns:
            tuple[str, int, dict]: (center_string, maximum_distance, metadata)
        """
        # Save initial state to history if enabled
        if self.save_history:
            self._save_history_entry(
                0,
                phase="initialization",
                parameters=self.params,
                message="Starting CSC algorithm",
            )

        self._report_progress(
            f"Starting CSC (d={self.params.get('d')}, n_blocks={self.params.get('n_blocks')})"
        )

        center = heuristic_closest_string(
            self.strings, d=self.params.get("d"), n_blocks=self.params.get("n_blocks")
        )

        if center:
            dist = max_distance(center, self.strings)
            metadata = {
                "iteracoes": 1,
                "parametros_usados": {
                    "d": self.params.get("d"),
                    "n_blocks": self.params.get("n_blocks"),
                },
                "center_found": center,
                "success": True,
            }

            # Save final state to history if enabled
            if self.save_history:
                self._save_history_entry(
                    1,
                    phase="completion",
                    best_fitness=dist,
                    best_solution=center,
                    message="CSC algorithm completed successfully",
                )

                # Add history to metadata
                metadata["history"] = self.get_history()

            return center, dist, metadata
        else:
            # Fallback para primeira string se falhar
            self._report_warning("CSC falhou, usando fallback para primeira string")
            fallback_center = self.strings[0]
            fallback_dist = max_distance(fallback_center, self.strings)

            metadata = {
                "iteracoes": 1,
                "parametros_usados": {
                    "d": self.params.get("d"),
                    "n_blocks": self.params.get("n_blocks"),
                },
                "center_found": fallback_center,
                "success": False,
                "fallback_used": True,
            }

            # Save final state to history if enabled (with fallback)
            if self.save_history:
                self._save_history_entry(
                    1,
                    phase="completion",
                    best_fitness=fallback_dist,
                    best_solution=fallback_center,
                    message="CSC algorithm failed, using fallback",
                )

                # Add history to metadata
                metadata["history"] = self.get_history()

            return fallback_center, fallback_dist, metadata
