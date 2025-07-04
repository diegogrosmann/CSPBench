"""
H³-CSP: Hybrid Hierarchical Hamming Search para CSP.

Classes:
    H3CSPAlgorithm: Wrapper para integração do H³-CSP ao framework CSP.
"""

from collections.abc import Callable

from algorithms.base import Algorithm, register_algorithm

from .config import H3_CSP_DEFAULTS
from .implementation import H3CSP


@register_algorithm
class H3CSPAlgorithm(Algorithm):
    """
    H³-CSP: Hybrid Hierarchical Hamming Search para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        set_progress_callback(callback): Define callback de progresso.
        run(): Executa o H³-CSP e retorna (centro, distância máxima).
    """

    name = "H³-CSP"
    default_params = H3_CSP_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}
        self.h3_csp_instance = H3CSP(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Passa o callback para a instância do H3CSP.

        Args:
            callback (Callable[[str], None]): Função de callback de progresso.
        """
        self.h3_csp_instance.set_progress_callback(callback)

    def run(self) -> tuple[str, int]:
        """
        Executa o H³-CSP e retorna a string central e a distância máxima.

        Returns:
            tuple[str, int]: (string_central, distancia_maxima)
        """
        return self.h3_csp_instance.run()
