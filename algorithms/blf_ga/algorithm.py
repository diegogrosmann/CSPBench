"""
BLF-GA: Blockwise Learning Fusion + Genetic Algorithm para CSP.

Classes:
    BLFGAAlgorithm: Wrapper para integração do BLF-GA ao framework CSP.
"""

from typing import Callable
from algorithms.base import Algorithm, register_algorithm
from .config import BLF_GA_DEFAULTS
from .implementation import BLFGA

@register_algorithm
class BLFGAAlgorithm(Algorithm):
    """
    BLF-GA: Blockwise Learning Fusion + Genetic Algorithm para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        set_progress_callback(callback): Define callback de progresso.
        run(): Executa o BLF-GA e retorna (centro, distância máxima).
    """
    name = "BLF-GA"
    default_params = BLF_GA_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        # Merge default params with user params
        self.params = {**self.default_params, **params}
        self.blf_ga_instance = BLFGA(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Passa o callback para a instância do BLFGA.

        Args:
            callback (Callable[[str], None]): Função de callback de progresso.
        """
        self.blf_ga_instance.set_progress_callback(callback)

    def run(self) -> tuple[str, int]:
        """
        Executa o BLF-GA e retorna a string central e a distância máxima.

        Returns:
            tuple[str, int]: (string_central, distancia_maxima)
        """
        return self.blf_ga_instance.run()
