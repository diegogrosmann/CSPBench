"""
Domain: CSP Algorithms

This module contains the core logic for Closest String Problem (CSP) algorithms,
including abstract interfaces, algorithm registry, and genetic operators.
Free from external dependencies following hexagonal architecture.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, TYPE_CHECKING, TypedDict

from .distance import DistanceCalculator
from .monitoring import AlgorithmMonitor


# =============================================================================
# ALGORITHM RESULT TYPES
# =============================================================================
class AlgorithmResult(TypedDict):
    """
    Complete algorithm result type.

    This represents the full structured result returned by CSP algorithms.
    """

    success: bool  # Whether the algorithm completed successfully
    center_string: str  # The center string solution
    max_distance: int  # Maximum distance from center to input strings
    parameters: Dict[str, Any]  # Algorithm parameters used
    error: str | None  # Error message if any
    metadata: Dict[str, Any]  # Detailed execution metadata


# =============================================================================
# ALGORITHM REGISTRY
# =============================================================================

global_registry: dict[str, type] = {}


def register_algorithm(cls: type) -> type:
    """
    Decorator to register an algorithm class in the global registry.

    Args:
        cls: Algorithm class to be registered

    Returns:
        type: The class itself, allowing use as decorator
    """
    algorithm_name = getattr(cls, "name", cls.__name__)
    global_registry[algorithm_name] = cls
    return cls


# =============================================================================
# ALGORITHM INTERFACES
# =============================================================================


class CSPAlgorithm(ABC):
    """Abstract base interface for all CSP algorithms."""

    # Required class attributes
    name: str
    default_params: dict = {}  # Default empty dict to ensure it exists
    supports_internal_parallel: bool = False

    def __init__(
        self,
        strings: list[str],
        alphabet: str,
        distance_calculator: DistanceCalculator,
        monitor: Optional[AlgorithmMonitor] = None,
        seed: Optional[int] = None,
        internal_jobs: int = 1,
        store: Any | None = None,  # adaptador simples usado nos testes legado
        **params,
    ):
        """
        Initialize algorithm with strings, alphabet, and dependencies.

        Args:
            strings: List of dataset strings
            alphabet: Used alphabet
            distance_calculator: Distance calculator instance
            monitor: AlgorithmMonitor para progresso/avisos/erros (opcional)
            global_seed: Global seed that overrides local seeds
            internal_jobs: Number of internal parallel jobs (default: 1)
            **params: Algorithm-specific parameters
        """
        self.strings = strings
        self.alphabet = alphabet
        # Se 'store' for passado (contrato legado dos testes de algoritmos), criamos
        # um adaptador leve para o protocolo AlgorithmMonitor. Isso permite que
        # os algoritmos usem self._monitor.on_progress / on_warning enquanto os
        # testes capturam via store.report_algorithm_progress / store.warning.

        self.params = {**self.default_params, **params}

        if monitor is None and store is not None:

            class _StoreMonitor:
                def __init__(self, _store):
                    self._store = _store

                # Compatível com AlgorithmMonitor
                def on_progress(
                    self, progress: float, message: str, /, **data: Any
                ) -> None:  # pragma: no cover - simples
                    rap = getattr(self._store, "report_algorithm_progress", None)
                    if callable(rap):
                        try:
                            rap(progress, message, **data)
                        except Exception:
                            pass

                def on_warning(
                    self, message: str, /, **data: Any
                ) -> None:  # pragma: no cover - simples
                    warn = getattr(self._store, "warning", None)
                    if callable(warn):
                        try:
                            warn(message)
                        except Exception:
                            pass

                def on_error(
                    self, message: str, exc: Exception | None = None, /, **data: Any
                ) -> None:  # pragma: no cover - simples
                    # Reaproveita warning caso exista
                    self.on_warning(f"ERROR: {message}")

                def is_cancelled(self) -> bool:  # pragma: no cover - simples
                    chk = getattr(self._store, "is_cancelled", None)
                    if callable(chk):
                        try:
                            return bool(chk())
                        except Exception:
                            return False
                    return False

            monitor = _StoreMonitor(store)

        self._monitor = monitor

        # Resource control
        self.internal_jobs = internal_jobs

        # Use provided distance calculator
        self._distance_calc = distance_calculator

        # Configure random generator with seed
        self.seed = seed
        self._setup_random_generator()

    def __getattr__(self, name: str):
        """
        Delegate attribute access to the distance calculator if the attribute
        is not found in the CSPAlgorithm instance.

        Args:
            name: Name of the attribute to access

        Returns:
            The attribute from the distance calculator instance
        """
        if hasattr(self._distance_calc, name):
            return getattr(self._distance_calc, name)
        elif self._monitor and hasattr(self._monitor, name):
            return getattr(self._monitor, name)
        else:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute '{name}'"
            )

    def _setup_random_generator(self) -> None:
        """Configure random generator with provided seed."""
        import random

        if self.seed is not None:
            self.rng = random.Random(self.seed)
            # Configure numpy if available
            try:
                import numpy as np

                np.random.seed(self.seed)
            except ImportError:
                pass
        else:
            self.rng = random.Random()

    def _get_timestamp(self) -> float:
        """Return current timestamp for history."""
        import time

        return time.time()

    # ------------------------------------------------------------------
    # Helper de aviso usado pelos algoritmos concretos (_report_warning)
    # ------------------------------------------------------------------
    def _report_warning(self, message: str, /, **data: Any) -> None:
        """Encapsula emissão de aviso sem depender rigidamente do monitor.

        Os algoritmos chamam self._report_warning(...) para registrar avisos
        sem precisar verificar existência de _monitor. Aqui garantimos que
        nenhuma exceção se propague dessa camada.
        """
        try:  # pragma: no branch - caminho simples
            if self._monitor and hasattr(self._monitor, "on_warning"):
                self._monitor.on_warning(message, **data)
        except Exception:  # pragma: no cover - segurança
            pass

    @abstractmethod
    def run(self) -> AlgorithmResult:
        """
        Execute algorithm and return structured result.

        Returns:
            AlgorithmResult: Dictionary containing:
                - center_string: The center string solution
                - max_distance: Maximum distance from center to input strings
                - parameters: Algorithm parameters used
                - metadata: Detailed execution metadata
        """

    def get_actual_params(self) -> dict[str, Any]:
        """Return actual parameters used (merged defaults and received)."""
        return self.params.copy()
