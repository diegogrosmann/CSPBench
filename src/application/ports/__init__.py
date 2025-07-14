"""
Módulo de Portas da Aplicação

Define interfaces (portas) que devem ser implementadas pela infraestrutura.
"""

from .repositories import AlgorithmRegistry, DatasetRepository, ExecutorPort, ExportPort

__all__ = [
    "DatasetRepository",
    "AlgorithmRegistry",
    "ExportPort",
    "ExecutorPort",
]
