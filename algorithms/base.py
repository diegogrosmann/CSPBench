"""
Base e registry global para algoritmos CSP.

Classes:
    CSPAlgorithm: Interface abstrata moderna para algoritmos CSP.
    Algorithm: Interface abstrata legacy (mantida para compatibilidade).

Funções:
    register_algorithm(cls): Decorador para registrar algoritmos no registry global.

Atributos:
    global_registry (dict): Registry global de algoritmos disponíveis.
"""

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Any

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


class CSPAlgorithm(ABC):
    """
    Classe base abstrata moderna para todos os algoritmos CSP.

    Esta é a nova interface padrão que suporta:
    - Execução com timeout
    - Callbacks de progresso e warning
    - Compatibilidade com ProcessPoolExecutor
    - Metadata estruturada

    Atributos obrigatórios:
        name (str): Nome do algoritmo
        default_params (dict): Parâmetros padrão
        is_deterministic (bool): Se é determinístico (padrão: False)
    """

    name: str
    default_params: dict
    is_deterministic: bool = False

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Inicializa o algoritmo com as strings e alfabeto.

        Args:
            strings: Lista de strings do dataset
            alphabet: Alfabeto utilizado
            **params: Parâmetros específicos do algoritmo
        """
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}
        self.progress_callback: Callable[[str], None] | None = None
        self.warning_callback: Callable[[str], None] | None = None

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Define callback para relatar progresso do algoritmo."""
        self.progress_callback = callback

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """Define callback para relatar warnings do algoritmo."""
        self.warning_callback = callback

    def _report_progress(self, message: str) -> None:
        """Relata progresso se callback estiver definido."""
        if self.progress_callback:
            self.progress_callback(message)

    def _report_warning(self, message: str) -> None:
        """Relata warning se callback estiver definido."""
        if self.warning_callback:
            self.warning_callback(message)

    @abstractmethod
    def run(self) -> tuple[str, int, dict[str, Any]]:
        """
        Executa o algoritmo.

        Returns:
            tuple contendo:
            - str: String central encontrada
            - int: Distância máxima (valor objetivo)
            - dict: Metadata da execução (iterações, etc.)
        """
        pass

    def get_metadata(self) -> dict[str, Any]:
        """
        Retorna metadados do algoritmo.

        Returns:
            Dicionário com informações sobre o algoritmo
        """
        return {
            "name": self.name,
            "params": self.params,
            "is_deterministic": self.is_deterministic,
            "input_size": len(self.strings),
            "string_length": len(self.strings[0]) if self.strings else 0,
            "alphabet_size": len(self.alphabet),
        }
