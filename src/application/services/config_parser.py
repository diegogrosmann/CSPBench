"""
Compat configuration parser for tests.

Provides minimal dataclasses and parsing/validation helpers expected by tests
in tests/unit/application/test_config_parser.py. This is independent from the
new unified config used elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from src.domain.errors import (
    BatchConfigurationError,
    OptimizationConfigurationError,
    SensitivityConfigurationError,
)


# -------------------- Dataclasses used in tests --------------------
@dataclass
class BatchMetadata:
    nome: str
    descricao: str
    autor: str
    versao: str
    data_criacao: str | None = None
    tags: List[str] = field(default_factory=list)
    timeout_global: int | None = None


@dataclass
class InfrastructureConfig:
    history: Dict[str, Any] = field(default_factory=dict)
    result: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ExportConfig:
    enabled: bool = True
    destination: str = "outputs"
    formats: Dict[str, Any] = field(default_factory=lambda: {"json": True})


@dataclass
class PlotsConfig:
    enabled: bool = True
    plot_convergence: bool = True
    style: str = "seaborn-v0_8"


@dataclass
class MonitoringConfig:
    enabled: bool = True
    interface: str = "simple"
    update_interval: int = 5


@dataclass
class LoggingConfig:
    level: str = "INFO"
    output: Dict[str, Any] = field(default_factory=lambda: {"console": True})


@dataclass
class SystemConfig:
    reproducibility: Dict[str, Any] = field(
        default_factory=lambda: {"global_seed": None, "strict_mode": False}
    )
    checkpointing: Dict[str, Any] = field(
        default_factory=lambda: {"enabled": False, "interval": 0}
    )


@dataclass
class OptimizationConfig:
    nome: str
    study_name: str
    direction: str
    n_trials: int
    timeout_per_trial: int
    target_datasets: List[str]
    target_algorithm: str
    parameters: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SensitivityConfig:
    nome: str
    analysis_method: str
    target_datasets: List[str]
    target_algorithm: str
    n_samples: int
    repetitions_per_sample: int
    parameters: Dict[str, Any] = field(default_factory=dict)
    output_metrics: List[str] = field(default_factory=list)


# Lightweight container used by ExperimentService (compat)
@dataclass
class BatchConfig:
    metadata: Any
    task: Any
    datasets: List[Any] = field(default_factory=list)
    algorithms: List[Any] = field(default_factory=list)
    experiment: Optional[Dict[str, Any]] = None
    optimization: Optional[Dict[str, Any]] = None
    sensitivity: Optional[Dict[str, Any]] = None
    export: ExportConfig = field(default_factory=ExportConfig)
    infrastructure: InfrastructureConfig = field(default_factory=InfrastructureConfig)
    plots: PlotsConfig = field(default_factory=PlotsConfig)
    monitoring: MonitoringConfig = field(default_factory=MonitoringConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    system: SystemConfig = field(default_factory=SystemConfig)
    resources: Optional[Any] = None


class ConfigurationParser:
    """Parser util expected by unit tests."""

    @staticmethod
    def load_file(path: str | Path) -> Dict[str, Any]:
        p = Path(path)
        if not p.exists():
            raise BatchConfigurationError("File not found")
        if p.suffix.lower() in {".yaml", ".yml"}:
            with p.open(encoding="utf-8") as fh:
                data = yaml.safe_load(fh) or {}
        elif p.suffix.lower() == ".json":
            with p.open(encoding="utf-8") as fh:
                data = json.load(fh)
        else:
            raise BatchConfigurationError("Unsupported format")
        if not isinstance(data, dict):
            raise BatchConfigurationError("Configuration root must be a mapping")
        return data

    # Aliases for ExperimentService
    load_config = load_file

    @staticmethod
    def parse_metadata(config: Dict[str, Any]) -> BatchMetadata:
        meta = config.get("metadados") or config.get("batch_info")
        if not meta:
            raise BatchConfigurationError("metadados/batch_info section not found")
        return BatchMetadata(
            nome=meta.get("nome", ""),
            descricao=meta.get("descricao", ""),
            autor=meta.get("autor", ""),
            versao=meta.get("versao", ""),
            data_criacao=meta.get("data_criacao"),
            tags=list(meta.get("tags", []) or []),
            timeout_global=meta.get("timeout_global"),
        )

    @staticmethod
    def parse_infrastructure_config(config: Dict[str, Any]) -> InfrastructureConfig:
        infra = config.get("infrastructure", {})
        return InfrastructureConfig(
            history=dict(infra.get("history", {})), result=dict(infra.get("result", {}))
        )

    @staticmethod
    def parse_export_config(config: Dict[str, Any]) -> ExportConfig:
        exp = config.get("export", {})
        return ExportConfig(
            enabled=bool(exp.get("enabled", True)),
            destination=str(exp.get("destination", "outputs")),
            formats=dict(exp.get("formats", {"json": True})),
        )

    @staticmethod
    def parse_plots_config(config: Dict[str, Any]) -> PlotsConfig:
        plots = config.get("plots", {})
        return PlotsConfig(
            enabled=bool(plots.get("enabled", True)),
            plot_convergence=bool(plots.get("plot_convergence", True)),
            style=str(plots.get("style", "seaborn-v0_8")),
        )

    @staticmethod
    def parse_monitoring_config(config: Dict[str, Any]) -> MonitoringConfig:
        mon = config.get("monitoring", {})
        return MonitoringConfig(
            enabled=bool(mon.get("enabled", True)),
            interface=str(mon.get("interface", "simple")),
            update_interval=int(mon.get("update_interval", 5)),
        )

    @staticmethod
    def parse_logging_config(config: Dict[str, Any]) -> LoggingConfig:
        log = config.get("logging", {})
        return LoggingConfig(
            level=str(log.get("level", "INFO")),
            output=dict(log.get("output", {"console": True})),
        )

    @staticmethod
    def parse_system_config(config: Dict[str, Any]) -> SystemConfig:
        sysc = config.get("system", {})
        return SystemConfig(
            reproducibility=dict(
                sysc.get("reproducibility", {"global_seed": None, "strict_mode": False})
            ),
            checkpointing=dict(
                sysc.get("checkpointing", {"enabled": False, "interval": 0})
            ),
        )

    @staticmethod
    def parse_resources_config(config: Dict[str, Any]) -> Dict[str, Any]:
        res = config.get("resources", {})
        par = res.get("parallel", {})
        return {
            "enabled": bool(par.get("enabled", True)),
            "max_workers": int(par.get("max_workers", 1)),
            "internal_jobs": int(par.get("internal_jobs", 1)),
        }

    @staticmethod
    def parse_optimization_configs(config: Dict[str, Any]) -> List[OptimizationConfig]:
        opt = config.get("optimization", {})
        lst = opt.get("optimizations", []) or []
        out: List[OptimizationConfig] = []
        for item in lst:
            out.append(
                OptimizationConfig(
                    nome=item.get("nome", item.get("name", "")),
                    study_name=item.get("study_name", "study"),
                    direction=item.get("direction", "minimize"),
                    n_trials=int(item.get("n_trials", item.get("trials", 0))),
                    timeout_per_trial=int(item.get("timeout_per_trial", 0)),
                    target_datasets=list(item.get("target_datasets", [])),
                    target_algorithm=item.get("target_algorithm", ""),
                    parameters=dict(item.get("parameters", {})),
                )
            )
        return out

    @staticmethod
    def parse_sensitivity_configs(config: Dict[str, Any]) -> List[SensitivityConfig]:
        sen = config.get("sensitivity", {})
        glob = sen.get("global_salib_config", {})
        lst = sen.get("analyses", []) or []
        out: List[SensitivityConfig] = []
        for item in lst:
            out.append(
                SensitivityConfig(
                    nome=item.get("nome", item.get("name", "")),
                    analysis_method=item.get("analysis_method", "morris"),
                    target_datasets=list(item.get("target_datasets", [])),
                    target_algorithm=item.get("target_algorithm", ""),
                    n_samples=int(glob.get("n_samples", 0)),
                    repetitions_per_sample=int(glob.get("repetitions_per_sample", 1)),
                    parameters=dict(item.get("parameters", {})),
                    output_metrics=list(item.get("output_metrics", [])),
                )
            )
        return out


# Simple facade used by ExperimentService (compat)
class ConfigParser:
    parse_config = staticmethod(ConfigurationParser.load_file)
    load_config = staticmethod(ConfigurationParser.load_file)


class ConfigurationValidator:
    @staticmethod
    def validate_batch_structure(config: Dict[str, Any]) -> str:
        # New structure with explicit task.type
        task = config.get("task", {})
        if task and "type" in task:
            t = task["type"]
            if t in {"experiment", "optimization", "sensitivity"}:
                # Minimal required sections for experiment
                if t == "experiment":
                    if "datasets" not in config or "algorithms" not in config:
                        raise BatchConfigurationError("Required sections missing")
                return t

        # Legacy: presence of 'experiments' list implies experiment
        if "experiments" in config:
            return "experiment"

        raise BatchConfigurationError("Configuration structure not recognized")

    @staticmethod
    def validate_optimization_config(config: OptimizationConfig) -> None:
        if not config.target_datasets:
            raise OptimizationConfigurationError("target_datasets list cannot be empty")
        if not isinstance(config.parameters, dict):
            raise OptimizationConfigurationError("parameters must be a mapping")

    @staticmethod
    def validate_sensitivity_config(config: SensitivityConfig) -> None:
        valid_methods = {"morris", "sobol", "fast", "delta"}
        if config.analysis_method not in valid_methods:
            raise SensitivityConfigurationError("Invalid analysis method")
        if not config.target_datasets:
            raise SensitivityConfigurationError("target_datasets list cannot be empty")
