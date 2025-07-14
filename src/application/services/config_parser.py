"""
Módulo de Parsing de Configuração CSPBench

Fornece classes e funções modulares para parsing e validação de configurações
de batch, otimização e análise de sensibilidade.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml

from src.domain.errors import (
    BatchConfigurationError,
    OptimizationConfigurationError,
    SensitivityConfigurationError,
)


@dataclass
class BatchMetadata:
    """Metadados do batch."""

    nome: str
    descricao: str
    autor: str
    versao: str
    data_criacao: str
    tags: List[str]
    timeout_global: int


@dataclass
class OptimizationConfig:
    """Configuração de uma otimização específica."""

    nome: str
    study_name: str
    direction: str
    n_trials: int
    timeout_per_trial: int
    target_datasets: List[str]
    target_algorithm: str
    parameters: Dict[str, Any]
    optuna_config: Optional[Dict[str, Any]] = None


@dataclass
class SensitivityConfig:
    """Configuração de uma análise de sensibilidade específica."""

    nome: str
    analysis_method: str
    target_datasets: List[str]
    target_algorithm: str
    n_samples: int
    repetitions_per_sample: int
    parameters: Dict[str, Any]
    output_metrics: List[str]
    method_config: Optional[Dict[str, Any]] = None


class ConfigurationParser:
    """Parser modular para configurações CSPBench."""

    @staticmethod
    def load_file(file_path: Union[str, Path]) -> Dict[str, Any]:
        """Carrega arquivo de configuração YAML ou JSON."""
        file_path = Path(file_path)

        if not file_path.exists():
            raise BatchConfigurationError(f"Arquivo não encontrado: {file_path}")

        with open(file_path, "r", encoding="utf-8") as f:
            if file_path.suffix.lower() in [".yaml", ".yml"]:
                config = yaml.safe_load(f)
            elif file_path.suffix.lower() == ".json":
                config = json.load(f)
            else:
                raise BatchConfigurationError(
                    f"Formato não suportado: {file_path.suffix}"
                )

        if not isinstance(config, dict):
            raise BatchConfigurationError("Configuração deve ser um objeto/dicionário")

        return config

    @staticmethod
    def parse_metadata(config: Dict[str, Any]) -> BatchMetadata:
        """Extrai metadados do batch."""
        # Suporte a ambos batch_info (legado) e metadados (novo)
        metadata_section = config.get("metadados", config.get("batch_info", {}))

        if not metadata_section:
            raise BatchConfigurationError("Seção metadados/batch_info não encontrada")

        return BatchMetadata(
            nome=metadata_section.get("nome", "sem nome"),
            descricao=metadata_section.get("descricao", ""),
            autor=metadata_section.get("autor", ""),
            versao=metadata_section.get("versao", "1.0"),
            data_criacao=metadata_section.get("data_criacao", ""),
            tags=metadata_section.get("tags", []),
            timeout_global=metadata_section.get("timeout_global", 3600),
        )

    @staticmethod
    def parse_optimization_configs(config: Dict[str, Any]) -> List[OptimizationConfig]:
        """Extrai configurações de otimização."""
        optimization_section = config.get("optimization", {})

        if not optimization_section:
            raise OptimizationConfigurationError("Seção optimization não encontrada")

        # Suporte a estrutura nova (optimizations) e legada (campos diretos)
        optimizations_list = optimization_section.get("optimizations")

        if optimizations_list is not None:
            # Nova estrutura: lista de otimizações
            return [
                ConfigurationParser._parse_single_optimization(opt_config)
                for opt_config in optimizations_list
            ]
        else:
            # Estrutura legada: configuração única
            # Converte para lista de uma otimização
            legacy_config = optimization_section.copy()

            # Garantir campos obrigatórios para estrutura legada
            if "target_dataset" in legacy_config:
                legacy_config["target_datasets"] = [legacy_config.pop("target_dataset")]

            if "nome" not in legacy_config:
                legacy_config["nome"] = (
                    f"Otimização {legacy_config.get('target_algorithm', 'Unknown')}"
                )

            return [ConfigurationParser._parse_single_optimization(legacy_config)]

    @staticmethod
    def _parse_single_optimization(opt_config: Dict[str, Any]) -> OptimizationConfig:
        """Parse uma única configuração de otimização."""
        required_fields = ["target_algorithm", "parameters"]
        for field in required_fields:
            if field not in opt_config:
                raise OptimizationConfigurationError(
                    f"Campo obrigatório ausente: {field}"
                )

        return OptimizationConfig(
            nome=opt_config.get("nome", "Otimização sem nome"),
            study_name=opt_config.get("study_name", "default_study"),
            direction=opt_config.get("direction", "minimize"),
            n_trials=opt_config.get("n_trials", 100),
            timeout_per_trial=opt_config.get("timeout_per_trial", 300),
            target_datasets=opt_config.get("target_datasets", []),
            target_algorithm=opt_config["target_algorithm"],
            parameters=opt_config["parameters"],
            optuna_config=opt_config.get("optuna_config"),
        )

    @staticmethod
    def parse_sensitivity_configs(config: Dict[str, Any]) -> List[SensitivityConfig]:
        """Extrai configurações de análise de sensibilidade."""
        sensitivity_section = config.get("sensitivity", {})

        if not sensitivity_section:
            raise SensitivityConfigurationError("Seção sensitivity não encontrada")

        # Extrair configurações globais
        global_config = sensitivity_section.get("global_salib_config", {})
        global_n_samples = global_config.get("n_samples", 1000)
        global_repetitions = global_config.get("repetitions_per_sample", 3)

        # Suporte a estrutura nova (analyses) e legada (campos diretos)
        analyses_list = sensitivity_section.get("analyses")

        if analyses_list is not None:
            # Nova estrutura: lista de análises
            return [
                ConfigurationParser._parse_single_sensitivity(
                    sens_config, global_n_samples, global_repetitions
                )
                for sens_config in analyses_list
            ]
        else:
            # Estrutura legada: configuração única
            legacy_config = sensitivity_section.copy()

            # Garantir campos obrigatórios para estrutura legada
            if "target_dataset" in legacy_config:
                legacy_config["target_datasets"] = [legacy_config.pop("target_dataset")]

            if "nome" not in legacy_config:
                legacy_config["nome"] = (
                    f"Análise {legacy_config.get('target_algorithm', 'Unknown')}"
                )

            return [
                ConfigurationParser._parse_single_sensitivity(
                    legacy_config, global_n_samples, global_repetitions
                )
            ]

    @staticmethod
    def _parse_single_sensitivity(
        sens_config: Dict[str, Any],
        global_n_samples: int = 1000,
        global_repetitions: int = 3,
    ) -> SensitivityConfig:
        """Parse uma única configuração de análise de sensibilidade."""
        required_fields = ["target_algorithm", "parameters"]
        for field in required_fields:
            if field not in sens_config:
                raise SensitivityConfigurationError(
                    f"Campo obrigatório ausente: {field}"
                )

        # Coletar configuração específica do método
        analysis_method = sens_config.get("analysis_method", "morris")
        method_config = None

        if analysis_method == "morris" and "morris" in sens_config:
            method_config = sens_config["morris"]
        elif analysis_method == "sobol" and "sobol" in sens_config:
            method_config = sens_config["sobol"]
        elif analysis_method == "fast" and "fast" in sens_config:
            method_config = sens_config["fast"]

        return SensitivityConfig(
            nome=sens_config.get("nome", "Análise sem nome"),
            analysis_method=analysis_method,
            target_datasets=sens_config.get("target_datasets", []),
            target_algorithm=sens_config["target_algorithm"],
            n_samples=sens_config.get("n_samples", global_n_samples),
            repetitions_per_sample=sens_config.get(
                "repetitions_per_sample", global_repetitions
            ),
            parameters=sens_config["parameters"],
            output_metrics=sens_config.get(
                "output_metrics", ["distance", "execution_time"]
            ),
            method_config=method_config,
        )

    @staticmethod
    def parse_resources_config(config: Dict[str, Any]) -> Dict[str, Any]:
        """Parse configurações de recursos do batch."""
        resources_section = config.get("resources", {})
        parallel_config = resources_section.get("parallel", {})

        return {
            "enabled": parallel_config.get("enabled", True),
            "max_workers": parallel_config.get("max_workers", 4),
            "internal_jobs": parallel_config.get("internal_jobs", 4),
            "backend": parallel_config.get("backend", "multiprocessing"),
            "timeouts": resources_section.get("timeouts", {}),
        }


class ConfigurationValidator:
    """Validador modular para configurações CSPBench."""

    @staticmethod
    def validate_batch_structure(config: Dict[str, Any]) -> str:
        """Valida e retorna o tipo de estrutura do batch."""
        # Detectar tipo de estrutura
        task_type = config.get("task", {}).get("type")

        if task_type in ["execution", "optimization", "sensitivity"]:
            # Nova estrutura unificada
            required_sections = ["metadados", "datasets", "algorithms", "task"]
            missing = [
                section for section in required_sections if section not in config
            ]

            if missing:
                # Verificar se é batch_info em vez de metadados (retrocompatibilidade)
                if "batch_info" in config and "metadados" in missing:
                    missing.remove("metadados")

            if missing:
                raise BatchConfigurationError(
                    f"Seções obrigatórias ausentes: {missing}"
                )

            return task_type

        elif "experiments" in config:
            # Estrutura legada de processamento
            return "execution"

        else:
            raise BatchConfigurationError(
                "Estrutura de configuração não reconhecida. "
                "Verifique se há seção 'task' com 'type' válido ou seção 'experiments'"
            )

    @staticmethod
    def validate_optimization_config(opt_config: OptimizationConfig) -> None:
        """Valida configuração de otimização."""
        if not opt_config.target_datasets:
            raise OptimizationConfigurationError(
                "Lista target_datasets não pode estar vazia"
            )

        if not opt_config.target_algorithm:
            raise OptimizationConfigurationError("Campo target_algorithm é obrigatório")

        if not opt_config.parameters:
            raise OptimizationConfigurationError(
                "Seção parameters não pode estar vazia"
            )

        # Validar parâmetros
        for param_name, param_config in opt_config.parameters.items():
            if not isinstance(param_config, dict):
                raise OptimizationConfigurationError(
                    f"Configuração do parâmetro {param_name} deve ser um objeto"
                )

            param_type = param_config.get("type")
            if param_type not in ["int", "uniform", "categorical"]:
                raise OptimizationConfigurationError(
                    f"Tipo de parâmetro inválido para {param_name}: {param_type}"
                )

    @staticmethod
    def validate_sensitivity_config(sens_config: SensitivityConfig) -> None:
        """Valida configuração de análise de sensibilidade."""
        if not sens_config.target_datasets:
            raise SensitivityConfigurationError(
                "Lista target_datasets não pode estar vazia"
            )

        if not sens_config.target_algorithm:
            raise SensitivityConfigurationError("Campo target_algorithm é obrigatório")

        if not sens_config.parameters:
            raise SensitivityConfigurationError("Seção parameters não pode estar vazia")

        # Validar método de análise
        valid_methods = ["morris", "sobol", "fast", "delta"]
        if sens_config.analysis_method not in valid_methods:
            raise SensitivityConfigurationError(
                f"Método de análise inválido: {sens_config.analysis_method}. "
                f"Métodos válidos: {valid_methods}"
            )

        # Validar parâmetros
        for param_name, param_config in sens_config.parameters.items():
            if not isinstance(param_config, dict):
                raise SensitivityConfigurationError(
                    f"Configuração do parâmetro {param_name} deve ser um objeto"
                )

            param_type = param_config.get("type")
            if param_type not in ["integer", "float", "categorical"]:
                raise SensitivityConfigurationError(
                    f"Tipo de parâmetro inválido para {param_name}: {param_type}"
                )
