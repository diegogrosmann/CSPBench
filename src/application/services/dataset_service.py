"""
Dataset Service - Simplified Dataset Loading.

Provides a unified interface for loading and generating datasets from
various sources including synthetic generation, file loading, and
external database retrieval.

This service acts as a facade that coordinates between different dataset
providers and ensures consistent dataset loading across the application.

Features:
- Synthetic dataset generation with configurable parameters
- File-based dataset loading with DATASET_DIRECTORY support
- External database integration (Entrez)
- Consistent return format (Dataset, parameters)

Public API:
    load_dataset(cfg) -> tuple[Dataset, dict]
        - synthetic: generates random dataset with optional noise
        - file: reads FASTA by logical name/path (uses DATASET_DIRECTORY if defined)
        - entrez: support via external module (dataset_entrez)
"""

from __future__ import annotations

from src.application.services.dataset_generator import SyntheticDatasetGenerator
from src.domain import Dataset
from src.domain.config import (
    DatasetAny,
    EntrezDatasetConfig,
    FileDatasetConfig,
    SyntheticDatasetConfig,
)
from src.infrastructure.external.dataset_entrez import EntrezDatasetDownloader
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository


def load_dataset(cfg: DatasetAny) -> tuple[Dataset, dict]:
    """
    Load/generate a Dataset from a DatasetConfig.
    
    This is the main entry point for dataset loading. It dispatches
    to the appropriate loader based on the configuration type.
    
    Args:
        cfg: Dataset configuration specifying type and parameters
        
    Returns:
        tuple: (Dataset object, parameters dict)
            - Dataset: Loaded or generated dataset
            - dict: Parameters used for loading/generation
            
    Raises:
        TypeError: If dataset configuration type is not supported
        
    Supported Types:
        - SyntheticDatasetConfig: Generates synthetic data
        - FileDatasetConfig: Loads from file system
        - EntrezDatasetConfig: Downloads from NCBI Entrez
    """
    if isinstance(cfg, SyntheticDatasetConfig):
        return SyntheticDatasetGenerator.generate_from_config(cfg)

    if isinstance(cfg, FileDatasetConfig):
        dataset, params = FileDatasetRepository.load(cfg.filename)
        if hasattr(cfg, "name") and cfg.name:
            dataset.name = cfg.name
        if hasattr(cfg, "id") and cfg.id:
            dataset.id = cfg.id
        return dataset, params

    if isinstance(cfg, EntrezDatasetConfig):
        dataset, params = EntrezDatasetDownloader.download(cfg)
        if hasattr(cfg, "name") and cfg.name:
            dataset.name = cfg.name
        if hasattr(cfg, "id") and cfg.id:
            dataset.id = cfg.id
        return dataset, params

    raise TypeError(f"Unsupported dataset type: {type(cfg)!r}")
