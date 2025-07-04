"""
H³-CSP: Hybrid Hierarchical Hamming Search para CSP.

Classes:
    H3CSPAlgorithm: Wrapper para integração do H³-CSP ao framework CSP.
"""

from collections.abc import Callable

from algorithms.base import CSPAlgorithm, register_algorithm

from .config import H3_CSP_DEFAULTS
from .implementation import H3CSP


@register_algorithm
class H3CSPAlgorithm(CSPAlgorithm):
    """
    H³-CSP: Hybrid Hierarchical Hamming Search para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        run(): Executa o H³-CSP e retorna (centro, distância máxima, metadata).
    """

    name = "H³-CSP"
    default_params = H3_CSP_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)
        self.h3_csp_instance = H3CSP(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Passa o callback para a instância do H3CSP.

        Args:
            callback (Callable[[str], None]): Função de callback de progresso.
        """
        super().set_progress_callback(callback)
        self.h3_csp_instance.set_progress_callback(callback)

    def run(self) -> tuple[str, int, dict]:
        """
        Executa o H³-CSP e retorna a string central, distância máxima e metadata.

        Returns:
            tuple[str, int, dict]: (string_central, distancia_maxima, metadata)
        """
        center, dist = self.h3_csp_instance.run()

        metadata = {
            "iteracoes": getattr(self.h3_csp_instance, "iterations", 1),
            "algoritmo": "H³-CSP",
            "parametros_usados": self.params,
            "centro_encontrado": center,
        }

        return center, dist, metadata
