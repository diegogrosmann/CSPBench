from algorithms.base import Algorithm, register_algorithm
from .config import CSC_DEFAULTS
from .implementation import heuristic_closest_string
from utils.distance import max_hamming

@register_algorithm
class CSCAlgorithm(Algorithm):
    """CSC: Consensus String Clustering"""
    name = "CSC"
    default_params = CSC_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}

    def run(self) -> tuple[str, int]:
        """Execute CSC algorithm"""
        center = heuristic_closest_string(
            self.strings, 
            d=self.params.get('d'),
            n_blocks=self.params.get('n_blocks')
        )
        if center:
            dist = max_hamming(center, self.strings)
            return center, dist
        else:
            # Fallback to first string if algorithm fails
            return self.strings[0], max_hamming(self.strings[0], self.strings)
