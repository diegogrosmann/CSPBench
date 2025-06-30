from typing import Callable
from algorithms.base import Algorithm, register_algorithm
from .config import DP_CSP_DEFAULTS
from .implementation import exact_dp_closest_string
from utils.distance import max_hamming

@register_algorithm
class DPCSPAlgorithm(Algorithm):
    """DP-CSP: Dynamic Programming exact solution for CSP"""
    name = "DP-CSP"
    default_params = DP_CSP_DEFAULTS
    is_deterministic = True

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}
        self.progress_callback = None

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Armazena o callback para ser usado no run()."""
        self.progress_callback = callback

    def run(self) -> tuple[str, int]:
        """Execute DP-CSP algorithm"""
        max_d = self.params.get('max_d')
        if max_d is None:
            # Use baseline as upper bound
            max_d = max_hamming(self.strings[0], self.strings)
        
        try:
            center, dist = exact_dp_closest_string(
                self.strings, 
                self.alphabet, 
                max_d,
                progress_callback=self.progress_callback
            )
            return center, dist
        except RuntimeError:
            # Algorithm failed, return fallback
            return self.strings[0], max_hamming(self.strings[0], self.strings)
