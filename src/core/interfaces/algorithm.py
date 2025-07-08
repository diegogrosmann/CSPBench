"""
Interface padronizada para algoritmos CSP.

Define o protocolo IAlgorithm que todos os algoritmos devem implementar
para garantir compatibilidade com o sistema de execução.
"""

from typing import Any, Callable, Optional, Protocol, TypedDict


class Result(TypedDict):
    """
    Resultado padronizado de execução de algoritmo.

    Attributes:
        center: String centro encontrada
        distance: Distância do centro
        metadata: Metadados adicionais (iterações, tempo, etc.)
    """

    center: str
    distance: float
    metadata: dict[str, Any]


class IAlgorithm(Protocol):
    """
    Interface padronizada para algoritmos CSP.

    Define os métodos obrigatórios que todos os algoritmos devem implementar
    para funcionar com o sistema de execução padronizado.
    """

    # Atributo opcional de classe
    is_deterministic: bool = False

    def run(self) -> Result:
        """
        Executa o algoritmo e retorna o resultado.

        Returns:
            Result: Resultado da execução com center, distance e metadata
        """
        ...

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define callback para atualizações de progresso.

        Args:
            callback: Função que recebe mensagens de progresso
        """
        ...

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define callback para avisos durante a execução.

        Args:
            callback: Função que recebe mensagens de aviso
        """
        ...


def is_algorithm_compatible(obj: Any) -> bool:
    """
    Verifica se um objeto implementa a interface IAlgorithm.

    Args:
        obj: Objeto a ser verificado

    Returns:
        bool: True se o objeto implementa a interface
    """
    return (
        hasattr(obj, "run")
        and callable(obj.run)
        and hasattr(obj, "set_progress_callback")
        and callable(obj.set_progress_callback)
        and hasattr(obj, "set_warning_callback")
        and callable(obj.set_warning_callback)
    )


def get_algorithm_deterministic(obj: Any) -> bool:
    """
    Obtém o atributo is_deterministic de um algoritmo.

    Args:
        obj: Objeto algoritmo

    Returns:
        bool: True se o algoritmo é determinístico
    """
    return getattr(obj, "is_deterministic", False)
