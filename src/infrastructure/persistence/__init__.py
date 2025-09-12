"""
Infrastructure Persistence Module

This module provides implementations for data storage and retrieval in CSPBench.
It includes repositories for datasets and algorithms, supporting file-based and
database-backed persistence strategies.

Available components:
    - FileDatasetRepository: File-based dataset storage using FASTA format
    - DomainAlgorithmRegistry: Algorithm registration and management
"""

from .algorithm_registry import DomainAlgorithmRegistry
from .dataset_repository import FileDatasetRepository

__all__ = [
    "FileDatasetRepository",
    "DomainAlgorithmRegistry",
]
