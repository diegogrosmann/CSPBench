"""
Infrastructure External Module

Adapters for external systems like NCBI, databases, etc.
"""

from .dataset_entrez import fetch_dataset

__all__ = [
    "fetch_dataset",
]
