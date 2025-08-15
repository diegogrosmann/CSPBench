"""
Infrastructure External Module

Adapters for external systems like NCBI, databases, etc.
"""

from .dataset_entrez import EntrezDatasetDownloader, ENTREZ_DEFAULTS

__all__ = [
    "EntrezDatasetDownloader",
    "ENTREZ_DEFAULTS",
]
