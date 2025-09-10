"""
Infrastructure External Module

Adapters for external systems like NCBI, databases, etc.
"""

from .dataset_entrez import ENTREZ_DEFAULTS, EntrezDatasetDownloader

__all__ = [
    "EntrezDatasetDownloader",
    "ENTREZ_DEFAULTS",
]
