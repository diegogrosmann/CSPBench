"""
Módulo de Persistência da Infraestrutura

Implementações para armazenamento e recuperação de dados.
"""

from .algorithm_registry import DomainAlgorithmRegistry
from .dataset_repository import FileDatasetRepository

__all__ = [
    "FileDatasetRepository",
    "DomainAlgorithmRegistry",
]
