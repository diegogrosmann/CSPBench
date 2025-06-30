from typing import Callable
from algorithms.base import Algorithm, register_algorithm
from .config import H3_CSP_DEFAULTS
from .implementation import H3CSP

@register_algorithm
class H3CSPAlgorithm(Algorithm):
    """H³-CSP: Hybrid Hierarchical Hamming Search"""
    name = "H³-CSP"
    default_params = H3_CSP_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}
        self.h3_csp_instance = H3CSP(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Passa o callback para a instância do H3CSP."""
        self.h3_csp_instance.set_progress_callback(callback)

    def run(self) -> tuple[str, int]:
        """Execute H³-CSP algorithm"""
        return self.h3_csp_instance.run()
