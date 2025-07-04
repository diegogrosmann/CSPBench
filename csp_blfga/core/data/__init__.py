"""
Módulo de dados para o sistema CSP-BLFGA.

Este módulo contém factories e utilitários para criação e gerenciamento
de datasets de diferentes tipos.
"""

from .dataset_factory import DatasetError, DatasetFactory, DatasetType

__all__ = ["DatasetFactory", "DatasetType", "DatasetError"]
