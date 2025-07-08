"""
Algoritmo de consenso ganancioso (Baseline) para o Closest String Problem.

Classes:
    BaselineAlg: Implementação do algoritmo baseline.
"""

from algorithms.base import CSPAlgorithm, register_algorithm

from .implementation import greedy_consensus, max_distance


@register_algorithm
class BaselineAlg(CSPAlgorithm):
    """
    Algoritmo de consenso ganancioso (Baseline) para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.

    Métodos:
        run(): Executa o algoritmo e retorna (centro, distância máxima, metadata).
    """

    name = "Baseline"
    default_params: dict = {}
    is_deterministic = True
    supports_internal_parallel = False  # Baseline não suporta paralelismo interno

    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Inicializa o algoritmo Baseline.

        Args:
            strings: Lista de strings do dataset
            alphabet: Alfabeto utilizado
            **params: Parâmetros específicos do algoritmo
        """
        super().__init__(strings, alphabet, **params)

    def run(self) -> tuple[str, int, dict]:
        """
        Executa o consenso guloso e calcula a maior distância de Hamming.

        Returns:
            tuple[str, int, dict]: (string centro, distância máxima, metadata)
        """
        self._report_progress("Iniciando consenso ganancioso...")
        center = greedy_consensus(self.strings, self.alphabet)

        self._report_progress("Calculando distância máxima...")
        dist = max_distance(center, self.strings)

        metadata = {
            "iteracoes": 1,  # Algoritmo determinístico executado uma vez
            "centro_encontrado": center,
            "total_strings": len(self.strings),
        }

        return center, dist, metadata
