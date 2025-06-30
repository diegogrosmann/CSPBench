from typing import Callable
from algorithms.base import Algorithm, register_algorithm
from .config import BLF_GA_DEFAULTS
from .implementation import BLFGA

@register_algorithm
class BLFGAAlgorithm(Algorithm):
    """BLF-GA: Blockwise Learning Fusion + Genetic Algorithm"""
    name = "BLF-GA"
    default_params = BLF_GA_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        # Merge default params with user params
        self.params = {**self.default_params, **params}
        self.blf_ga_instance = BLFGA(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Passa o callback para a instância do BLFGA."""
        self.blf_ga_instance.set_progress_callback(callback)

    def run(self) -> tuple[str, int]:
        """Execute BLF-GA algorithm"""
        return self.blf_ga_instance.run()
