"""
Carregador de configurações para sistema de execução em lote.

Classes:
    ConfigLoader: Carrega e valida configurações de batch
    ConfigError: Exceção para erros de configuração

Funcionalidades:
    - Carregamento de arquivos YAML
    - Validação de esquema de configuração
    - Resolução de paths relativos
    - Merge de configurações
    - Validação de parâmetros
"""

import logging
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)


class ConfigError(Exception):
    """Exceção personalizada para erros de configuração."""


class ConfigLoader:
    """
    Carregador e validador de configurações de batch.

    Responsável por carregar arquivos YAML, validar esquemas,
    resolver paths e fazer merge de configurações.
    """

    def __init__(self, base_path: str | Path | None = None):
        """
        Inicializa o carregador de configurações.

        Args:
            base_path: Caminho base para resolução de paths relativos
        """
        self.base_path = Path(base_path) if base_path else Path.cwd()
        self.required_fields = {
            "nome",
            "algoritmos",
            "dataset",
            "runs_per_algorithm_per_base",
            "num_bases",
        }
        self.optional_fields = {
            "timeout",
            "max_workers",
            "seed",
            "export_format",
            "output_dir",
            "verbose",
        }

    def load_config(self, config_path: str | Path) -> dict[str, Any]:
        """
        Carrega configuração de um arquivo YAML.

        Args:
            config_path: Caminho para o arquivo de configuração

        Returns:
            Dicionário com a configuração carregada

        Raises:
            ConfigError: Se houve erro no carregamento ou validação
        """
        config_path = Path(config_path)

        if not config_path.is_absolute():
            config_path = self.base_path / config_path

        if not config_path.exists():
            raise ConfigError(f"Arquivo de configuração não encontrado: {config_path}")

        logger.info("Carregando configuração: %s", config_path)

        try:
            with open(config_path, encoding="utf-8") as f:
                config = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ConfigError(f"Erro ao parsear YAML: {e}") from e
        except Exception as e:
            raise ConfigError(f"Erro ao ler arquivo: {e}") from e

        # Resolver paths relativos
        config = self._resolve_paths(config, config_path.parent)

        # Validar configuração
        self._validate_config(config)

        logger.info("Configuração carregada com sucesso: %s", config.get("nome", "sem nome"))
        return config

    def load_multiple_configs(self, config_paths: list[str | Path]) -> list[dict[str, Any]]:
        """
        Carrega múltiplas configurações.

        Args:
            config_paths: Lista de caminhos para arquivos de configuração

        Returns:
            Lista de dicionários com configurações carregadas

        Raises:
            ConfigError: Se houve erro no carregamento
        """
        configs = []
        for path in config_paths:
            try:
                config = self.load_config(path)
                configs.append(config)
            except ConfigError as e:
                logger.error("Erro ao carregar %s: %s", path, e)
                raise

        return configs

    def merge_configs(self, base_config: dict[str, Any], override_config: dict[str, Any]) -> dict[str, Any]:
        """
        Faz merge de duas configurações.

        Args:
            base_config: Configuração base
            override_config: Configuração para sobrescrever

        Returns:
            Configuração resultante do merge
        """
        merged = base_config.copy()
        merged.update(override_config)

        # Validar configuração resultante
        self._validate_config(merged)

        return merged

    def _resolve_paths(self, config: dict[str, Any], config_dir: Path) -> dict[str, Any]:
        """
        Resolve paths relativos na configuração.

        Args:
            config: Configuração a ser processada
            config_dir: Diretório do arquivo de configuração

        Returns:
            Configuração com paths resolvidos
        """
        resolved = config.copy()

        # Resolver paths em campos específicos
        path_fields = ["output_dir"]

        for field in path_fields:
            if field in resolved and resolved[field]:
                path = Path(resolved[field])
                if not path.is_absolute():
                    resolved[field] = str(config_dir / path)

        # Resolver paths em dataset
        if "dataset" in resolved and isinstance(resolved["dataset"], dict):
            dataset = resolved["dataset"].copy()
            if "filepath" in dataset and dataset["filepath"]:
                path = Path(dataset["filepath"])
                if not path.is_absolute():
                    dataset["filepath"] = str(config_dir / path)
            resolved["dataset"] = dataset

        return resolved

    def _validate_config(self, config: dict[str, Any]) -> None:
        """
        Valida se a configuração está correta.

        Args:
            config: Configuração a ser validada

        Raises:
            ConfigError: Se a configuração estiver inválida
        """
        if not isinstance(config, dict):
            raise ConfigError("Configuração deve ser um dicionário")

        # Verificar campos obrigatórios com suporte a retro-compatibilidade
        missing_fields = self.required_fields - config.keys()

        # Verificar se existe runs_per_algorithm_per_base ou execucoes_por_algoritmo_por_base
        if "runs_per_algorithm_per_base" in missing_fields:
            if "execucoes_por_algoritmo_por_base" in config:
                missing_fields.remove("runs_per_algorithm_per_base")
            elif "execucoes_por_algoritmo" in config:
                missing_fields.remove("runs_per_algorithm_per_base")

        if missing_fields:
            raise ConfigError(f"Campos obrigatórios ausentes: {missing_fields}")

        # Verificar campos desconhecidos
        all_fields = self.required_fields | self.optional_fields
        # Adicionar campos antigos permitidos para compatibilidade
        all_fields.add("execucoes_por_algoritmo_por_base")
        all_fields.add("execucoes_por_algoritmo")
        unknown_fields = config.keys() - all_fields
        if unknown_fields:
            logger.warning("Campos desconhecidos ignorados: %s", unknown_fields)

        # Validações específicas
        self._validate_algorithms(config["algoritmos"])
        self._validate_dataset(config["dataset"])
        self._validate_numeric_fields(config)

    def _validate_algorithms(self, algorithms: Any) -> None:
        """Valida lista de algoritmos."""
        if not isinstance(algorithms, list):
            raise ConfigError("Campo 'algoritmos' deve ser uma lista")

        if not algorithms:
            raise ConfigError("Lista de algoritmos não pode estar vazia")

        for alg in algorithms:
            if not isinstance(alg, str):
                raise ConfigError("Nomes de algoritmos devem ser strings")

    def _validate_dataset(self, dataset: Any) -> None:
        """Valida configuração de dataset."""
        if not isinstance(dataset, dict):
            raise ConfigError("Campo 'dataset' deve ser um dicionário")

        if "tipo" not in dataset:
            raise ConfigError("Dataset deve ter campo 'tipo'")

        valid_types = ["file", "entrez", "synthetic"]
        if dataset["tipo"] not in valid_types:
            raise ConfigError(f"Tipo de dataset deve ser um de: {valid_types}")

        # Validações específicas por tipo
        if dataset["tipo"] == "file":
            if "filepath" not in dataset:
                raise ConfigError("Dataset tipo 'file' deve ter campo 'filepath'")

        elif dataset["tipo"] == "entrez":
            required_entrez = ["email", "db", "term"]
            missing = [f for f in required_entrez if f not in dataset]
            if missing:
                raise ConfigError(f"Dataset tipo 'entrez' deve ter campos: {missing}")

        elif dataset["tipo"] == "synthetic":
            required_synthetic = ["n", "L", "alphabet"]
            missing = [f for f in required_synthetic if f not in dataset]
            if missing:
                raise ConfigError(f"Dataset tipo 'synthetic' deve ter campos: {missing}")

    def _validate_numeric_fields(self, config: dict[str, Any]) -> None:
        """Valida campos numéricos."""
        numeric_fields = {
            "runs_per_algorithm_per_base": (int, lambda x: x > 0),
            "execucoes_por_algoritmo_por_base": (int, lambda x: x > 0),  # Retro-compatibilidade
            "num_bases": (int, lambda x: x > 0),
            "timeout": (int, lambda x: x > 0),
            "max_workers": (int, lambda x: x > 0),
            "seed": (int, lambda x: x >= 0),
        }

        for field, (expected_type, validator) in numeric_fields.items():
            if field in config:
                value = config[field]
                if not isinstance(value, expected_type):
                    raise ConfigError(f"Campo '{field}' deve ser do tipo {expected_type.__name__}")
                if not validator(value):
                    raise ConfigError(f"Campo '{field}' tem valor inválido: {value}")

    def create_default_config(self) -> dict[str, Any]:
        """
        Cria uma configuração padrão.

        Returns:
            Configuração padrão
        """
        return {
            "nome": "config_padrao",
            "algoritmos": ["Baseline", "BLF-GA"],
            "dataset": {
                "tipo": "file",
                "filepath": "saved_datasets/dataset_custom.fasta",
            },
            "runs_per_algorithm_per_base": 3,
            "num_bases": 1,
            "timeout": 300,
            "max_workers": 4,
            "seed": 42,
            "export_format": "csv",
            "output_dir": "results",
            "verbose": False,
        }

    def save_config(self, config: dict[str, Any], output_path: str | Path) -> None:
        """
        Salva configuração em arquivo YAML.

        Args:
            config: Configuração a ser salva
            output_path: Caminho do arquivo de saída
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Validar antes de salvar
        self._validate_config(config)

        try:
            with open(output_path, "w", encoding="utf-8") as f:
                yaml.safe_dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
            logger.info("Configuração salva em: %s", output_path)
        except Exception as e:
            raise ConfigError(f"Erro ao salvar configuração: {e}") from e
