"""
CSC: Consensus String Clustering para CSP.

Classes:
    CSCAlgorithm: Implementação do algoritmo CSC.
"""

from algorithms.base import CSPAlgorithm, register_algorithm
from src.utils.distance import max_distance

from .config import CSC_DEFAULTS
from .implementation import heuristic_closest_string


@register_algorithm
class CSCAlgorithm(CSPAlgorithm):
    """
    CSC: Consensus String Clustering para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        run(): Executa o CSC e retorna (centro, distância máxima, metadata).
    """

    name = "CSC"
    default_params = CSC_DEFAULTS
    supports_internal_parallel = False  # CSC não suporta paralelismo interno

    def run(self) -> tuple[str, int, dict]:
        """
        Executa o algoritmo CSC e retorna a string central, distância máxima e metadata.

        Returns:
            tuple[str, int, dict]: (string_central, distancia_maxima, metadata)
        """
        self._report_progress(f"Iniciando CSC (d={self.params.get('d')}, n_blocks={self.params.get('n_blocks')})")

        center = heuristic_closest_string(self.strings, d=self.params.get("d"), n_blocks=self.params.get("n_blocks"))

        if center:
            dist = max_distance(center, self.strings)
            metadata = {
                "iteracoes": 1,
                "parametros_usados": {
                    "d": self.params.get("d"),
                    "n_blocks": self.params.get("n_blocks"),
                },
                "centro_encontrado": center,
                "sucesso": True,
            }
            return center, dist, metadata
        else:
            # Fallback para primeira string se falhar
            self._report_warning("CSC falhou, usando fallback para primeira string")
            fallback_center = self.strings[0]
            fallback_dist = max_distance(fallback_center, self.strings)

            metadata = {
                "iteracoes": 1,
                "parametros_usados": {
                    "d": self.params.get("d"),
                    "n_blocks": self.params.get("n_blocks"),
                },
                "centro_encontrado": fallback_center,
                "sucesso": False,
                "fallback_usado": True,
            }

            return fallback_center, fallback_dist, metadata
