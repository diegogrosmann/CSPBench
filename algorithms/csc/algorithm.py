"""
CSC: Consensus String Clustering para CSP.

Classes:
    CSCAlgorithm: Implementação do algoritmo CSC.
"""
from algorithms.base import Algorithm, register_algorithm
from .config import CSC_DEFAULTS
from .implementation import heuristic_closest_string
from utils.distance import max_hamming

@register_algorithm
class CSCAlgorithm(Algorithm):
    """
    CSC: Consensus String Clustering para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        run(): Executa o CSC e retorna (centro, distância máxima).
    """
    name = "CSC"
    default_params = CSC_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}

    def run(self) -> tuple[str, int]:
        """
        Executa o algoritmo CSC e retorna a string central e a distância máxima.

        Returns:
            tuple[str, int]: (string_central, distancia_maxima)
        """
        center = heuristic_closest_string(
            self.strings, 
            d=self.params.get('d'),
            n_blocks=self.params.get('n_blocks')
        )
        if center:
            dist = max_hamming(center, self.strings)
            return center, dist
        else:
            # Fallback para primeira string se falhar
            return self.strings[0], max_hamming(self.strings[0], self.strings)
