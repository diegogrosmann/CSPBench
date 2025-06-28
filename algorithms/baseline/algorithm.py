from algorithms.base import Algorithm, register_algorithm
from .implementation import greedy_consensus, max_distance

@register_algorithm
class BaselineAlg(Algorithm):
    """Algoritmo consenso ganancioso como baseline."""
    name = "Baseline"
    default_params: dict = {}
    is_deterministic = True

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet

    def run(self) -> tuple[str, int]:
        center = greedy_consensus(self.strings, self.alphabet)
        dist = max_distance(center, self.strings)
        return center, dist
