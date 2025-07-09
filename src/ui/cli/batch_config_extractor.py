"""
Extrator centralizado para configurações de batch unificado.

Este módulo processa arquivos de configuração YAML no formato padronizado
e extrai informações para execução, otimização e análise de sensibilidade.
"""

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import yaml

from src.datasets.dataset_entrez import fetch_dataset
from src.datasets.dataset_file import load_dataset
from src.datasets.dataset_synthetic import generate_dataset_from_params

logger = logging.getLogger(__name__)


class BatchConfigExtractor:
    """Extrator centralizado para configurações de batch."""

    def __init__(self, config_path: str):
        """
        Inicializa o extrator com o caminho do arquivo de configuração.

        Args:
            config_path: Caminho para o arquivo de configuração YAML
        """
        self.config_path = config_path
        self.config = self._load_config()
        self._validate_config()

    def _load_config(self) -> Dict[str, Any]:
        """Carrega configuração do arquivo YAML."""
        if not os.path.exists(self.config_path):
            raise FileNotFoundError(
                f"Arquivo de configuração não encontrado: {self.config_path}"
            )

        try:
            with open(self.config_path, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)
            return config
        except yaml.YAMLError as e:
            raise ValueError(f"Erro ao carregar YAML: {e}")

    def _validate_config(self) -> None:
        """Valida estrutura básica da configuração."""
        required_sections = ["batch_info", "datasets", "task", "algorithms"]

        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Seção obrigatória '{section}' não encontrada")

        # Validar tipo de tarefa
        task_type = self.config["task"].get("type")
        if task_type not in ["execution", "sensitivity", "optimization"]:
            raise ValueError(f"Tipo de tarefa inválido: {task_type}")

        # Validar se a seção específica da tarefa existe
        if task_type not in self.config["task"]:
            raise ValueError(
                f"Seção '{task_type}' não encontrada na configuração da tarefa"
            )

    def get_batch_info(self) -> Dict[str, Any]:
        """Retorna informações do batch."""
        return self.config["batch_info"]

    def get_datasets(self) -> Dict[str, Dict[str, Any]]:
        """Retorna configurações dos datasets."""
        return {ds["id"]: ds for ds in self.config["datasets"]}

    def get_task_type(self) -> str:
        """Retorna o tipo de tarefa."""
        return self.config["task"]["type"]

    def get_algorithms(self) -> List[str]:
        """Retorna lista de algoritmos."""
        return self.config["algorithms"]

    def get_algorithm_params(self) -> Dict[str, Any]:
        """Retorna parâmetros fixos dos algoritmos."""
        return self.config.get("algorithm_params", {})

    def get_output_config(self) -> Dict[str, Any]:
        """Retorna configuração de saída."""
        return self.config.get("output", {})

    def get_advanced_config(self) -> Dict[str, Any]:
        """Retorna configurações avançadas."""
        return self.config.get("advanced", {})

    def generate_dataset(
        self, dataset_id: str, silent: bool = True
    ) -> Tuple[List[str], Dict[str, Any]]:
        """
        Gera dataset baseado no ID fornecido.

        Args:
            dataset_id: ID do dataset na configuração
            silent: Se True, não exibe mensagens

        Returns:
            Tupla (sequências, parâmetros)
        """
        datasets = self.get_datasets()

        if dataset_id not in datasets:
            raise ValueError(f"Dataset '{dataset_id}' não encontrado na configuração")

        dataset_config = datasets[dataset_id]
        dataset_type = dataset_config.get("tipo", "synthetic")
        params = dataset_config.get("parametros", {})

        if dataset_type == "synthetic":
            seqs, generated_params = generate_dataset_from_params(
                n=params.get("n", 20),
                L=params.get("L", 100),
                alphabet=params.get("alphabet", "ACGT"),
                noise=params.get("noise", 0.1),
                fully_random=params.get("fully_random", False),
                seed=params.get("seed", None),
            )

            batch_params = {"dataset_source": "1", "dataset_id": dataset_id}
            batch_params.update(generated_params)

            return seqs, batch_params

        elif dataset_type == "file":
            # Configurar arquivo se especificado
            if "filename" in params:
                os.environ["DATASET_FILE"] = params["filename"]

            seqs, file_params = load_dataset(silent=silent)

            batch_params = {"dataset_source": "2", "dataset_id": dataset_id}
            batch_params.update(file_params)

            return seqs, batch_params

        elif dataset_type == "entrez":
            # Configurar parâmetros Entrez
            if "query" in params:
                os.environ["ENTREZ_QUERY"] = params["query"]
            if "db" in params:
                os.environ["ENTREZ_DB"] = params["db"]
            if "retmax" in params:
                os.environ["ENTREZ_RETMAX"] = str(params["retmax"])

            seqs, entrez_params = fetch_dataset()

            batch_params = {"dataset_source": "3", "dataset_id": dataset_id}
            batch_params.update(entrez_params)

            return seqs, batch_params

        else:
            raise ValueError(f"Tipo de dataset não suportado: {dataset_type}")

    def get_execution_configs(self) -> List[Dict[str, Any]]:
        """Retorna configurações de execução."""
        if self.get_task_type() != "execution":
            raise ValueError("Configuração não é do tipo 'execution'")

        return self.config["task"]["execution"]["executions"]

    def get_sensitivity_configs(self) -> List[Dict[str, Any]]:
        """Retorna configurações de análise de sensibilidade."""
        if self.get_task_type() != "sensitivity":
            raise ValueError("Configuração não é do tipo 'sensitivity'")

        return self.config["task"]["sensitivity"]["analyses"]

    def get_optimization_configs(self) -> List[Dict[str, Any]]:
        """Retorna configurações de otimização."""
        if self.get_task_type() != "optimization":
            raise ValueError("Configuração não é do tipo 'optimization'")

        return self.config["task"]["optimization"]["studies"]

    def get_global_timeout(self) -> int:
        """Retorna timeout global."""
        return self.config["batch_info"].get("timeout_global", 1800)

    def should_use_curses(self) -> bool:
        """Retorna se deve usar interface curses."""
        return self.get_advanced_config().get("use_curses", True)

    def get_parallel_config(self) -> Dict[str, Any]:
        """Retorna configuração de paralelismo."""
        return self.get_advanced_config().get("parallel", {"n_jobs": -1})

    def get_log_level(self) -> str:
        """Retorna nível de log."""
        return self.get_advanced_config().get("log_level", "INFO")
