"""
Implementações das interfaces padronizadas.
"""

from .algorithm_adapter import AlgorithmAdapter
from .parallel_executor import ParallelExecutor
from .simple_console import SimpleConsole

__all__ = [
    "ParallelExecutor",
    "AlgorithmAdapter",
    "SimpleConsole",
]
