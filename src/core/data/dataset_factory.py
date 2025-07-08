"""
Factory para criação de datasets de diferentes tipos.

Classes:
    DatasetFactory: Factory para criação de datasets
    DatasetType: Enum com tipos de datasets suportados

Funcionalidades:
    - Criação de datasets a partir de configurações
    - Suporte a múltiplos tipos (file, entrez, synthetic)
    - Validação de parâmetros
    - Cache de datasets
    - Logging detalhado
"""

import logging
from enum import Enum
from pathlib import Path
from typing import Any

from src.datasets.dataset_entrez import fetch_dataset_silent
from src.datasets.dataset_file import load_dataset_with_params
from src.datasets.dataset_synthetic import generate_dataset_with_params

logger = logging.getLogger(__name__)


class DatasetType(Enum):
    """Enum para tipos de datasets suportados."""

    FILE = "file"
    ENTREZ = "entrez"
    SYNTHETIC = "synthetic"


class DatasetError(Exception):
    """Exceção personalizada para erros de dataset."""


class DatasetFactory:
    """
    Factory para criação de datasets de diferentes tipos.

    Centraliza a criação de datasets baseado em configurações,
    fornecendo uma interface unificada para todos os tipos.
    """

    def __init__(self, cache_enabled: bool = True):
        """
        Inicializa a factory de datasets.

        Args:
            cache_enabled: Se deve usar cache para datasets
        """
        self.cache_enabled = cache_enabled
        self._cache: dict[str, tuple[list[str], dict[str, Any]]] = {}

    def create_dataset(
        self, dataset_config: dict[str, Any], base_index: int | None = None
    ) -> tuple[list[str], dict[str, Any]]:
        """
        Cria um dataset baseado na configuração.

        Args:
            dataset_config: Configuração do dataset
            base_index: Índice da base (para múltiplas bases)

        Returns:
            Tupla (sequências, metadados)

        Raises:
            DatasetError: Se houve erro na criação do dataset
        """
        if not isinstance(dataset_config, dict):
            raise DatasetError("Configuração de dataset deve ser um dicionário")

        if "tipo" not in dataset_config:
            raise DatasetError("Configuração de dataset deve ter campo 'tipo'")

        dataset_type = dataset_config["tipo"]

        # Verificar cache
        cache_key = self._get_cache_key(dataset_config, base_index)
        if self.cache_enabled and cache_key in self._cache:
            logger.debug("Usando dataset do cache: %s", cache_key)
            return self._cache[cache_key]

        logger.info("Criando dataset tipo '%s' (base %s)", dataset_type, base_index)

        try:
            if dataset_type == DatasetType.FILE.value:
                sequences, metadata = self._create_file_dataset(dataset_config)
            elif dataset_type == DatasetType.ENTREZ.value:
                sequences, metadata = self._create_entrez_dataset(
                    dataset_config, base_index
                )
            elif dataset_type == DatasetType.SYNTHETIC.value:
                sequences, metadata = self._create_synthetic_dataset(
                    dataset_config, base_index
                )
            else:
                raise DatasetError(f"Tipo de dataset não suportado: {dataset_type}")

            # Adicionar metadata adicional
            metadata["dataset_type"] = dataset_type
            metadata["base_index"] = base_index
            metadata["cache_key"] = cache_key

            # Salvar no cache
            if self.cache_enabled:
                self._cache[cache_key] = (sequences, metadata)

            logger.info(
                "Dataset criado: %d sequências, L=%s",
                len(sequences),
                metadata.get("L", "N/A"),
            )
            return sequences, metadata

        except Exception as e:
            logger.error("Erro ao criar dataset: %s", e)
            raise DatasetError(f"Erro ao criar dataset: {e}") from e

    def _create_file_dataset(
        self, config: dict[str, Any]
    ) -> tuple[list[str], dict[str, Any]]:
        """
        Cria dataset a partir de arquivo.

        Args:
            config: Configuração do dataset

        Returns:
            Tupla (sequências, metadados)
        """
        if "filepath" not in config:
            raise DatasetError("Dataset tipo 'file' deve ter campo 'filepath'")

        filepath = config["filepath"]
        if not Path(filepath).exists():
            raise DatasetError(f"Arquivo não encontrado: {filepath}")

        params = {"filepath": filepath}
        sequences, metadata = load_dataset_with_params(params)

        return sequences, metadata

    def _create_entrez_dataset(
        self, config: dict[str, Any], base_index: int | None = None
    ) -> tuple[list[str], dict[str, Any]]:
        """
        Cria dataset a partir do NCBI Entrez.

        Args:
            config: Configuração do dataset
            base_index: Índice da base (para seed)

        Returns:
            Tupla (sequências, metadados)
        """
        required_fields = ["email", "db", "term"]
        missing_fields = [f for f in required_fields if f not in config]
        if missing_fields:
            raise DatasetError(
                f"Dataset tipo 'entrez' deve ter campos: {missing_fields}"
            )

        params = config.copy()

        # Adicionar seed baseado no índice da base
        if base_index is not None and "seed" not in params:
            params["seed"] = 42 + base_index

        sequences, metadata = fetch_dataset_silent(params)

        return sequences, metadata

    def _create_synthetic_dataset(
        self, config: dict[str, Any], base_index: int | None = None
    ) -> tuple[list[str], dict[str, Any]]:
        """
        Cria dataset sintético.

        Args:
            config: Configuração do dataset
            base_index: Índice da base (para seed)

        Returns:
            Tupla (sequências, metadados)
        """
        required_fields = ["n", "L", "alphabet"]
        missing_fields = [f for f in required_fields if f not in config]
        if missing_fields:
            raise DatasetError(
                f"Dataset tipo 'synthetic' deve ter campos: {missing_fields}"
            )

        params = config.copy()

        # Adicionar seed baseado no índice da base
        if base_index is not None:
            base_seed = params.get("seed", 42)
            params["seed"] = base_seed + base_index

        sequences, metadata = generate_dataset_with_params(params)

        return sequences, metadata

    def _get_cache_key(
        self, config: dict[str, Any], base_index: int | None = None
    ) -> str:
        """
        Gera chave de cache para o dataset.

        Args:
            config: Configuração do dataset
            base_index: Índice da base

        Returns:
            Chave de cache
        """
        # Incluir campos relevantes na chave
        key_parts = [
            config.get("tipo", "unknown"),
            str(base_index) if base_index is not None else "0",
        ]

        # Adicionar campos específicos por tipo
        if config.get("tipo") == DatasetType.FILE.value:
            key_parts.append(config.get("filepath", ""))
        elif config.get("tipo") == DatasetType.ENTREZ.value:
            key_parts.extend(
                [
                    config.get("db", ""),
                    config.get("term", ""),
                    str(config.get("n", 0)),
                    str(config.get("seed", 0)),
                ]
            )
        elif config.get("tipo") == DatasetType.SYNTHETIC.value:
            key_parts.extend(
                [
                    str(config.get("n", 0)),
                    str(config.get("L", 0)),
                    config.get("alphabet", ""),
                    str(config.get("seed", 0)),
                    str(config.get("noise", 0)),
                ]
            )

        return "_".join(key_parts)

    def clear_cache(self) -> None:
        """Limpa o cache de datasets."""
        self._cache.clear()
        logger.info("Cache de datasets limpo")

    def get_cache_info(self) -> dict[str, Any]:
        """
        Retorna informações sobre o cache.

        Returns:
            Dicionário com informações do cache
        """
        return {
            "enabled": self.cache_enabled,
            "size": len(self._cache),
            "keys": list(self._cache.keys()),
        }

    def validate_config(self, config: dict[str, Any]) -> None:
        """
        Valida configuração de dataset.

        Args:
            config: Configuração a ser validada

        Raises:
            DatasetError: Se a configuração estiver inválida
        """
        if not isinstance(config, dict):
            raise DatasetError("Configuração de dataset deve ser um dicionário")

        if "tipo" not in config:
            raise DatasetError("Configuração de dataset deve ter campo 'tipo'")

        dataset_type = config["tipo"]
        if dataset_type not in [t.value for t in DatasetType]:
            raise DatasetError(f"Tipo de dataset não suportado: {dataset_type}")

        # Validações específicas por tipo
        if dataset_type == DatasetType.FILE.value:
            if "filepath" not in config:
                raise DatasetError("Dataset tipo 'file' deve ter campo 'filepath'")

        elif dataset_type == DatasetType.ENTREZ.value:
            required_fields = ["email", "db", "term"]
            missing_fields = [f for f in required_fields if f not in config]
            if missing_fields:
                raise DatasetError(
                    f"Dataset tipo 'entrez' deve ter campos: {missing_fields}"
                )

        elif dataset_type == DatasetType.SYNTHETIC.value:
            required_fields = ["n", "L", "alphabet"]
            missing_fields = [f for f in required_fields if f not in config]
            if missing_fields:
                raise DatasetError(
                    f"Dataset tipo 'synthetic' deve ter campos: {missing_fields}"
                )

    def create_multiple_datasets(
        self, dataset_config: dict[str, Any], num_bases: int
    ) -> list[tuple[list[str], dict[str, Any]]]:
        """
        Cria múltiplos datasets (bases diferentes).

        Args:
            dataset_config: Configuração base do dataset
            num_bases: Número de bases a criar

        Returns:
            Lista de tuplas (sequências, metadados)
        """
        datasets = []
        for i in range(num_bases):
            try:
                dataset = self.create_dataset(dataset_config, base_index=i)
                datasets.append(dataset)
            except Exception as e:
                logger.error("Erro ao criar base %s: %s", i, e)
                raise DatasetError(f"Erro ao criar base {i}: {e}") from e

        return datasets

    def get_dataset_info(self, dataset_config: dict[str, Any]) -> dict[str, Any]:
        """
        Retorna informações sobre um dataset sem carregá-lo.

        Args:
            dataset_config: Configuração do dataset

        Returns:
            Dicionário com informações do dataset
        """
        info = {
            "tipo": dataset_config.get("tipo", "unknown"),
            "valid": False,
            "error": None,
        }

        try:
            self.validate_config(dataset_config)
            info["valid"] = True

            # Informações específicas por tipo
            if dataset_config["tipo"] == DatasetType.FILE.value:
                filepath = dataset_config.get("filepath", "")
                path = Path(filepath)
                info.update(
                    {
                        "filepath": filepath,
                        "exists": path.exists(),
                        "size": path.stat().st_size if path.exists() else 0,
                    }
                )

            elif dataset_config["tipo"] == DatasetType.ENTREZ.value:
                info.update(
                    {
                        "db": dataset_config.get("db", ""),
                        "term": dataset_config.get("term", ""),
                        "n": dataset_config.get("n", 0),
                    }
                )

            elif dataset_config["tipo"] == DatasetType.SYNTHETIC.value:
                info.update(
                    {
                        "n": dataset_config.get("n", 0),
                        "L": dataset_config.get("L", 0),
                        "alphabet": dataset_config.get("alphabet", ""),
                        "noise": dataset_config.get("noise", 0),
                    }
                )

        except DatasetError as e:
            info["error"] = str(e)

        return info
