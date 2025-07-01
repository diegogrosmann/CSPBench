"""
DP-CSP: Dynamic Programming exata para o Closest String Problem.

Classes:
    DPCSPAlgorithm: Implementação do algoritmo exato por programação dinâmica.
"""
from typing import Callable
from algorithms.base import Algorithm, register_algorithm
from .config import DP_CSP_DEFAULTS
from .implementation import exact_dp_closest_string
from utils.distance import max_hamming

@register_algorithm
class DPCSPAlgorithm(Algorithm):
    """
    DP-CSP: Solução exata por programação dinâmica para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        set_progress_callback(callback): Define callback de progresso.
        run(): Executa o DP-CSP e retorna (centro, distância máxima).
    """
    name = "DP-CSP"
    default_params = DP_CSP_DEFAULTS
    is_deterministic = True

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}
        self.progress_callback = None

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Armazena o callback para ser usado no run().

        Args:
            callback (Callable[[str], None]): Função de callback de progresso.
        """
        self.progress_callback = callback

    def run(self) -> tuple[str, int]:
        """
        Executa o algoritmo DP-CSP e retorna a string central e a distância máxima.

        Returns:
            tuple[str, int]: (string_central, distancia_maxima)
        """
        max_d = self.params.get('max_d')
        if max_d is None:
            # Usa baseline como upper bound
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
            # Fallback para primeira string se falhar
            return self.strings[0], max_hamming(self.strings[0], self.strings)
