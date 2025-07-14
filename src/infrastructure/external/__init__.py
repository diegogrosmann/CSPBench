"""
Módulo External da Infraestrutura

Adaptadores para sistemas externos como NCBI, bases de dados, etc.
"""

from .dataset_entrez import fetch_dataset

__all__ = [
    "fetch_dataset",
]
