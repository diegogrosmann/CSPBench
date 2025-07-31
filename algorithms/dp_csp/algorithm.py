"""
DP-CSP: Exact Dynamic Programming for the Closest String Problem.

Classes:
    DPCSPAlgorithm: Exact dynamic programming algorithm implementation.
"""

from src.domain.algorithms import CSPAlgorithm, register_algorithm
from src.domain.metrics import max_distance

from .config import DP_CSP_DEFAULTS
from .implementation import exact_dp_closest_string


@register_algorithm
class DPCSPAlgorithm(CSPAlgorithm):
    """
    DP-CSP: Exact dynamic programming solution for the Closest String Problem.

    Args:
        strings (list[str]): List of input strings.
        alphabet (str): Alphabet used.
        **params: Algorithm parameters.

    Methods:
        run(): Executes DP-CSP and returns (center, maximum distance, metadata).
    """

    name = "DP-CSP"
    default_params = DP_CSP_DEFAULTS
    supports_internal_parallel = False  # DP-CSP doesn't support internal parallelism
    is_deterministic = True

    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Initialize DP-CSP algorithm.

        Args:
            strings: List of dataset strings
            alphabet: Alphabet used
            **params: Algorithm-specific parameters
        """
        super().__init__(strings, alphabet, **params)

    def run(self) -> tuple[str, int, dict]:
        """
        Execute DP-CSP algorithm and return center string, maximum distance and metadata.

        Returns:
            tuple[str, int, dict]: (center_string, maximum_distance, metadata)
        """
        # Save initial state to history if enabled
        if self.save_history:
            self._save_history_entry(
                0,
                phase="initialization",
                parameters=self.params,
                message="Starting DP-CSP algorithm",
            )

        max_d = self.params.get("max_d")
        if max_d is None:
            # Use baseline as upper bound
            max_d = max_distance(self.strings[0], self.strings)

        self._report_progress(f"Starting DP-CSP with max_d={max_d}")

        try:
            center, dist = exact_dp_closest_string(
                self.strings,
                self.alphabet,
                max_d,
                progress_callback=self._report_progress,
            )

            metadata = {
                "iteracoes": 1,
                "max_d_usado": max_d,
                "solucao_exata": True,
                "center_found": center,
            }

            # Save final state to history if enabled
            if self.save_history:
                self._save_history_entry(
                    1,
                    phase="completion",
                    best_fitness=dist,
                    best_solution=center,
                    message="DP-CSP algorithm completed successfully",
                )

                # Add history to metadata
                metadata["history"] = self.get_history()

            return center, dist, metadata
        except RuntimeError as e:
            # DP-CSP failed due to memory/time limitations
            self._report_warning(f"DP-CSP failed: {e}")

            # Save failure state to history if enabled
            if self.save_history:
                self._save_history_entry(
                    1, phase="error", message=f"DP-CSP algorithm failed: {e}"
                )

            # Re-raise the exception to indicate real failure
            raise RuntimeError(f"DP-CSP failed: {e}") from e
