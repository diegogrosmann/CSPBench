"""
DP-CSP: Dynamic Programming exata para o Closest String Problem.

Classes:
    DPCSPAlgorithm: Implementação do algoritmo exato por programação dinâmica.
"""


from algorithms.base import CSPAlgorithm, register_algorithm
from csp_blfga.utils.distance import max_distance

from .config import DP_CSP_DEFAULTS
from .implementation import exact_dp_closest_string


@register_algorithm
class DPCSPAlgorithm(CSPAlgorithm):
    """
    DP-CSP: Solução exata por programação dinâmica para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        run(): Executa o DP-CSP e retorna (centro, distância máxima, metadata).
    """

    name = "DP-CSP"
    default_params = DP_CSP_DEFAULTS
    is_deterministic = True

    def __init__(self, strings: list[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)

    def run(self) -> tuple[str, int, dict]:
        """
        Executa o algoritmo DP-CSP e retorna a string central, distância máxima e metadata.

        Returns:
            tuple[str, int, dict]: (string_central, distancia_maxima, metadata)
        """
        max_d = self.params.get("max_d")
        if max_d is None:
            # Usa baseline como upper bound
            max_d = max_distance(self.strings[0], self.strings)

        self._report_progress(f"Iniciando DP-CSP com max_d={max_d}")

        try:
            center, dist = exact_dp_closest_string(
                self.strings,
                self.alphabet,
                max_d,
                progress_callback=self._report_progress,
            )

            metadata = {
                "iteracoes": 1,
                "max_d_usado": max_d,
                "solucao_exata": True,
                "centro_encontrado": center,
            }

            return center, dist, metadata
        except RuntimeError:
            # Fallback para primeira string se falhar
            self._report_warning("DP-CSP falhou, usando fallback para primeira string")
            fallback_center = self.strings[0]
            fallback_dist = max_distance(fallback_center, self.strings)

            metadata = {
                "iteracoes": 1,
                "max_d_usado": max_d,
                "solucao_exata": False,
                "fallback_usado": True,
                "centro_encontrado": fallback_center,
            }

            return fallback_center, fallback_dist, metadata
