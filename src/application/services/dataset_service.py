"""
Dataset Service simplificado

Apenas uma função pública: load_dataset(cfg) -> Dataset
- synthetic: gera dataset aleatório com ruído opcional
- file: lê FASTA por nome lógico/caminho (usa DATASET_DIRECTORY se definido)
- entrez: suporte via módulo externo (dataset_entrez)
"""

from __future__ import annotations

from src.application.services.dataset_generator import SyntheticDatasetGenerator
from src.infrastructure.external.dataset_entrez import EntrezDatasetDownloader
from src.domain import Dataset
from src.domain.config import (
    DatasetAny,
    EntrezDatasetConfig,
    FileDatasetConfig,
    SyntheticDatasetConfig,
)

from src.infrastructure.persistence.dataset_repository import FileDatasetRepository
from src.domain.errors import DatasetNotFoundError


def load_dataset(cfg: DatasetAny) -> tuple[Dataset, dict]:
    """
    Carrega/gera um Dataset a partir de um DatasetConfig.
    """
    if isinstance(cfg, SyntheticDatasetConfig):
        return SyntheticDatasetGenerator.generate_from_config(cfg)

    if isinstance(cfg, FileDatasetConfig):
        dataset, params = FileDatasetRepository.load(cfg.filename)
        if hasattr(cfg, "name") and cfg.name:
            dataset.name = cfg.name
        return dataset, params

    if isinstance(cfg, EntrezDatasetConfig):
        dataset, params = EntrezDatasetDownloader.download(cfg)
        if hasattr(cfg, "name") and cfg.name:
            dataset.name = cfg.name
        return dataset, params

    raise TypeError(f"Tipo de dataset não suportado: {type(cfg)!r}")
