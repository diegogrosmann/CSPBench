"""
Algoritmo de consenso ganancioso (Baseline) para o Closest String Problem.

Classes:
    BaselineAlg: Implementação do algoritmo baseline.
"""

from src.domain.algorithms import CSPAlgorithm, register_algorithm

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
        import time

        start_time = time.time()

        self._report_progress("Iniciando consenso ganancioso...")

        # Salvar estado inicial no histórico
        self._save_history_entry(
            0,
            phase="initialization",
            total_strings=len(self.strings),
            string_length=len(self.strings[0]) if self.strings else 0,
        )

        center = greedy_consensus(self.strings, self.alphabet)

        self._report_progress("Calculando distância máxima...")
        dist = max_distance(center, self.strings)

        end_time = time.time()
        execution_time = end_time - start_time

        # Salvar estado final no histórico
        self._save_history_entry(
            1,
            phase="completion",
            center_found=center,
            max_distance=dist,
            execution_time=execution_time,
        )

        metadata = {
            "iteracoes": 1,  # Algoritmo determinístico executado uma vez
            "centro_encontrado": center,
            "total_strings": len(self.strings),
            "execution_time": execution_time,
            "history": self.get_history() if self.save_history else [],
        }

        return center, dist, metadata
