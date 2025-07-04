"""
Base e registry global para algoritmos CSP.

Classes:
    Algorithm: Interface abstrata para algoritmos CSP.

Funções:
    register_algorithm(cls): Decorador para registrar algoritmos no registry global.

Atributos:
    global_registry (dict): Registry global de algoritmos disponíveis.
"""

from abc import ABC, abstractmethod
from collections.abc import Callable

# Registry for all algorithms
global_registry: dict[str, type] = {}


def register_algorithm(cls: type) -> type:
    """
    Decorador para registrar uma classe de algoritmo no registry global.

    Args:
        cls (type): Classe do algoritmo a ser registrada.
    Returns:
        type: A própria classe, para uso como decorador.
    """
    name = getattr(cls, "name", cls.__name__)
    global_registry[name] = cls
    return cls


class Algorithm(ABC):
    """
    Classe base abstrata para todos os algoritmos CSP.

    Cada algoritmo deve definir:
        - name (str): Nome do algoritmo.
        - default_params (dict): Parâmetros padrão.
        - is_deterministic (bool, opcional): Se é determinístico.

    Métodos obrigatórios:
        - __init__(self, strings, alphabet, **params)
        - run(self) -> tuple[str, int]
    """

    name: str
    default_params: dict
    is_deterministic: bool = False

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Inicializa o algoritmo com as strings e alfabeto.

        Args:
            strings (list[str]): Lista de strings do dataset.
            alphabet (str): Alfabeto utilizado.
            **params: Parâmetros específicos do algoritmo.
        """
        pass

    def set_progress_callback(
        self, callback: Callable[[str], None]
    ) -> None:  # noqa: B027
        """
        (Opcional) Define um callback para relatar o progresso do algoritmo.

        Args:
            callback (Callable[[str], None]): Função de callback.
        """
        # Implementação opcional - subclasses podem sobrescrever
        pass

    @abstractmethod
    def run(self) -> tuple[str, int]:
        """
        Executa o algoritmo, retornando a string central e a distância máxima.

        Returns:
            tuple[str, int]: (string_central, distancia_maxima)
        """
        pass
