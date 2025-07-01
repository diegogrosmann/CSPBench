"""
Algoritmo de consenso ganancioso (Baseline) para o Closest String Problem.

Classes:
    BaselineAlg: Implementação do algoritmo baseline.
"""
from algorithms.base import Algorithm, register_algorithm
from .implementation import greedy_consensus, max_distance

@register_algorithm
class BaselineAlg(Algorithm):
    """
    Algoritmo de consenso ganancioso (Baseline) para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.

    Métodos:
        run(): Executa o algoritmo e retorna (centro, distância máxima).
    """
    name = "Baseline"
    default_params: dict = {}
    is_deterministic = True

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet

    def run(self) -> tuple[str, int]:
        """
        Executa o consenso guloso e calcula a maior distância de Hamming.

        Returns:
            tuple[str, int]: (string centro, distância máxima)
        """
        center = greedy_consensus(self.strings, self.alphabet)
        dist = max_distance(center, self.strings)
        return center, dist
