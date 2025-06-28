from abc import ABC, abstractmethod
from typing import Callable

"""
Base e registry global para algoritmos CSP.

Classes:
    Algorithm: Interface abstrata para algoritmos CSP.

Funções:
    register_algorithm(cls): Decorador para registrar algoritmos no registry global.

Atributos:
    global_registry (dict): Registry global de algoritmos disponíveis.
"""

# Registry for all algorithms
global_registry: dict[str, type] = {}

def register_algorithm(cls: type) -> type:
    """
    Decorator to register an algorithm class in the global registry.
    Usage:
        @register_algorithm
        class MyAlg(Algorithm): ...
    """
    name = getattr(cls, 'name', cls.__name__)
    global_registry[name] = cls
    return cls

class Algorithm(ABC):
    """
    Base class for all CSP algorithms.
    Each algorithm must set:
      - name: str
      - default_params: dict
      - is_deterministic: bool (optional, default False)

    Must implement:
      - __init__(self, strings: list[str], alphabet: str, **params)
      - run(self) -> tuple[str, int]
    """
    name: str
    default_params: dict
    is_deterministic: bool = False

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str, **params):
        pass

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Opcional: Define um callback para relatar o progresso."""
        pass

    @abstractmethod
    def run(self) -> tuple[str, int]:
        """Execute the algorithm, returning (center_string, max_distance)"""
        pass
