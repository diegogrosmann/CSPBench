"""
Infrastructure External Module.

Adapters and clients for external systems integration including NCBI databases,
web services, and other third-party data sources. This module provides
infrastructure-level components for accessing external resources.

Components:
    - EntrezDatasetDownloader: NCBI Entrez API client for sequence downloads
    - ENTREZ_DEFAULTS: Default configuration for NCBI connections

Features:
    - External API integration with proper error handling
    - Configurable connection parameters and authentication
    - Data transformation from external formats to domain objects
    - Caching and rate limiting for external service calls

The external module follows the hexagonal architecture pattern, providing
adapters that translate between external service APIs and internal domain
models, ensuring loose coupling and testability.
"""

from .dataset_entrez import ENTREZ_DEFAULTS, EntrezDatasetDownloader

__all__ = [
    "EntrezDatasetDownloader",
    "ENTREZ_DEFAULTS",
]
