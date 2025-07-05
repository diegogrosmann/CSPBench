"""
BLF-GA: Blockwise Learning Fusion + Genetic Algorithm para CSP.

Classes:
    BLFGAAlgorithm: Wrapper para integração do BLF-GA ao framework CSP.
"""

from collections.abc import Callable

from algorithms.base import CSPAlgorithm, register_algorithm

from .config import BLF_GA_DEFAULTS
from .implementation import BLFGA


@register_algorithm
class BLFGAAlgorithm(CSPAlgorithm):
    """
    BLF-GA: Blockwise Learning Fusion + Genetic Algorithm para o Closest String Problem.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        run(): Executa o BLF-GA e retorna (centro, distância máxima, metadata).
    """

    name = "BLF-GA"
    default_params = BLF_GA_DEFAULTS
    supports_internal_parallel = True  # BLF-GA pode usar paralelismo interno

    def __init__(self, strings: list[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)
        self.blf_ga_instance = BLFGA(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Passa o callback para a instância do BLFGA.

        Args:
            callback (Callable[[str], None]): Função de callback de progresso.
        """
        super().set_progress_callback(callback)
        self.blf_ga_instance.set_progress_callback(callback)

    def run(self) -> tuple[str, int, dict]:
        """
        Executa o BLF-GA e retorna a string central, a distância máxima e metadata.
        """
        best, best_val, history = self.blf_ga_instance.run()

        metadata = {
            "iteracoes": len(history) if history else 0,
            "melhor_distancia": best_val,
            "historico_completo": len(history) if history else 0,
        }

        return best, best_val, metadata

    def run_with_history(self) -> tuple[str, int, list]:
        """
        Executa o BLF-GA e retorna a string central, a distância máxima e o histórico de distâncias.
        """
        return self.blf_ga_instance.run()
