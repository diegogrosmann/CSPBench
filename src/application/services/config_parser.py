"""
Compatibility Configuration Parser for Tests.

Provides minimal dataclasses and parsing/validation helpers expected by tests
in tests/unit/application/test_config_parser.py. This is independent from the
new unified config used elsewhere in the application.

This module maintains backward compatibility with existing test infrastructure
while providing a clean separation from the main configuration system.

Features:
    - Legacy configuration parsing
    - Batch metadata extraction
    - Infrastructure configuration parsing
    - Optimization and sensitivity analysis configuration
    - Validation helpers for different configuration types
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
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
    """
    Batch metadata configuration for legacy compatibility.
    
    Attributes:
        nome (str): Batch name.
        descricao (str): Batch description.
        autor (str): Batch author.
        versao (str): Batch version.
        data_criacao (Optional[str]): Creation date.
        tags (List[str]): List of tags.
        timeout_global (Optional[int]): Global timeout in seconds.
    """
    
    nome: str
    descricao: str
    autor: str
    versao: str
    data_criacao: str | None = None
    tags: List[str] = field(default_factory=list)
    timeout_global: int | None = None


@dataclass
class InfrastructureConfig:
    """
    Infrastructure configuration for legacy compatibility.
    
    Attributes:
        history (Dict[str, Any]): History configuration.
        result (Dict[str, Any]): Result configuration.
    """
    
    history: Dict[str, Any] = field(default_factory=dict)
    result: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ExportConfig:
    """
    Export configuration for legacy compatibility.
    
    Attributes:
        enabled (bool): Whether export is enabled.
        destination (str): Export destination path.
        formats (Dict[str, Any]): Export format configuration.
    """
    
    enabled: bool = True
    destination: str = "outputs"
    formats: Dict[str, Any] = field(default_factory=lambda: {"json": True})


@dataclass
class PlotsConfig:
    """
    Plots configuration for legacy compatibility.
    
    Attributes:
        enabled (bool): Whether plotting is enabled.
        plot_convergence (bool): Whether to plot convergence.
        style (str): Plot style to use.
    """
    
    enabled: bool = True
    plot_convergence: bool = True
    style: str = "seaborn-v0_8"


@dataclass
class MonitoringConfig:
    """
    Monitoring configuration for legacy compatibility.
    
    Attributes:
        enabled (bool): Whether monitoring is enabled.
        interface (str): Monitoring interface type.
        update_interval (int): Update interval in seconds.
    """
    
    enabled: bool = True
    interface: str = "simple"
    update_interval: int = 5


@dataclass
class LoggingConfig:
    """
    Logging configuration for legacy compatibility.
    
    Attributes:
        level (str): Logging level.
        output (Dict[str, Any]): Output configuration.
    """
    
    level: str = "INFO"
    output: Dict[str, Any] = field(default_factory=lambda: {"console": True})


@dataclass
class SystemConfig:
    """
    System configuration for legacy compatibility.
    
    Attributes:
        reproducibility (Dict[str, Any]): Reproducibility settings.
        checkpointing (Dict[str, Any]): Checkpointing settings.
    """
    
    reproducibility: Dict[str, Any] = field(
        default_factory=lambda: {"global_seed": None, "strict_mode": False}
    )
    checkpointing: Dict[str, Any] = field(
        default_factory=lambda: {"enabled": False, "interval": 0}
    )


@dataclass
class OptimizationConfig:
    """
    Optimization configuration for legacy compatibility.
    
    Attributes:
        nome (str): Optimization name.
        study_name (str): Study name.
        direction (str): Optimization direction ("minimize" or "maximize").
        n_trials (int): Number of trials.
        timeout_per_trial (int): Timeout per trial in seconds.
        target_datasets (List[str]): Target dataset names.
        target_algorithm (str): Target algorithm name.
        parameters (Dict[str, Any]): Optimization parameters.
    """
    
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
    """
    Sensitivity analysis configuration for legacy compatibility.
    
    Attributes:
        nome (str): Analysis name.
        analysis_method (str): Analysis method ("morris", "sobol", etc.).
        target_datasets (List[str]): Target dataset names.
        target_algorithm (str): Target algorithm name.
        n_samples (int): Number of samples.
        repetitions_per_sample (int): Repetitions per sample.
        parameters (Dict[str, Any]): Analysis parameters.
        output_metrics (List[str]): Output metrics to analyze.
    """
    
    nome: str
    analysis_method: str
    target_datasets: List[str]
    target_algorithm: str
    n_samples: int
    repetitions_per_sample: int
    parameters: Dict[str, Any] = field(default_factory=dict)
    output_metrics: List[str] = field(default_factory=list)


@dataclass
class BatchConfig:
    """
    Lightweight batch configuration container for ExperimentService compatibility.
    
    Attributes:
        metadata: Batch metadata.
        task: Task configuration.
        datasets (List[Any]): Dataset configurations.
        algorithms (List[Any]): Algorithm configurations.
        experiment (Optional[Dict[str, Any]]): Experiment configuration.
        optimization (Optional[Dict[str, Any]]): Optimization configuration.
        sensitivity (Optional[Dict[str, Any]]): Sensitivity configuration.
        export (ExportConfig): Export configuration.
        infrastructure (InfrastructureConfig): Infrastructure configuration.
        plots (PlotsConfig): Plots configuration.
        monitoring (MonitoringConfig): Monitoring configuration.
        logging (LoggingConfig): Logging configuration.
        system (SystemConfig): System configuration.
        resources (Optional[Any]): Resource configuration.
    """
    
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
    """
    Parser utility expected by unit tests.
    
    Provides static methods for parsing various configuration file formats
    and extracting different configuration sections for legacy compatibility.
    """

    @staticmethod
    def load_file(path: str | Path) -> Dict[str, Any]:
        """
        Load configuration from file.
        
        Supports YAML and JSON formats with automatic format detection
        based on file extension.
        
        Args:
            path (Union[str, Path]): Path to configuration file.
            
        Returns:
            Dict[str, Any]: Parsed configuration data.
            
        Raises:
            BatchConfigurationError: If file not found, unsupported format,
                                   or invalid structure.
        """
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
        """
        Extract batch metadata from configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            BatchMetadata: Extracted metadata.
            
        Raises:
            BatchConfigurationError: If metadata section not found.
        """
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
        """
        Extract infrastructure configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            InfrastructureConfig: Extracted infrastructure configuration.
        """
        infra = config.get("infrastructure", {})
        return InfrastructureConfig(
            history=dict(infra.get("history", {})), result=dict(infra.get("result", {}))
        )

    @staticmethod
    def parse_export_config(config: Dict[str, Any]) -> ExportConfig:
        """
        Extract export configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            ExportConfig: Extracted export configuration.
        """
        exp = config.get("export", {})
        return ExportConfig(
            enabled=bool(exp.get("enabled", True)),
            destination=str(exp.get("destination", "outputs")),
            formats=dict(exp.get("formats", {"json": True})),
        )

    @staticmethod
    def parse_plots_config(config: Dict[str, Any]) -> PlotsConfig:
        """
        Extract plots configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            PlotsConfig: Extracted plots configuration.
        """
        plots = config.get("plots", {})
        return PlotsConfig(
            enabled=bool(plots.get("enabled", True)),
            plot_convergence=bool(plots.get("plot_convergence", True)),
            style=str(plots.get("style", "seaborn-v0_8")),
        )

    @staticmethod
    def parse_monitoring_config(config: Dict[str, Any]) -> MonitoringConfig:
        """
        Extract monitoring configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            MonitoringConfig: Extracted monitoring configuration.
        """
        mon = config.get("monitoring", {})
        return MonitoringConfig(
            enabled=bool(mon.get("enabled", True)),
            interface=str(mon.get("interface", "simple")),
            update_interval=int(mon.get("update_interval", 5)),
        )

    @staticmethod
    def parse_logging_config(config: Dict[str, Any]) -> LoggingConfig:
        """
        Extract logging configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            LoggingConfig: Extracted logging configuration.
        """
        log = config.get("logging", {})
        return LoggingConfig(
            level=str(log.get("level", "INFO")),
            output=dict(log.get("output", {"console": True})),
        )

    @staticmethod
    def parse_system_config(config: Dict[str, Any]) -> SystemConfig:
        """
        Extract system configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            SystemConfig: Extracted system configuration.
        """
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
        """
        Extract resources configuration.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            Dict[str, Any]: Resources configuration containing parallel settings.
        """
        res = config.get("resources", {})
        par = res.get("parallel", {})
        return {
            "enabled": bool(par.get("enabled", True)),
            "max_workers": int(par.get("max_workers", 1)),
            "internal_jobs": int(par.get("internal_jobs", 1)),
        }

    @staticmethod
    def parse_optimization_configs(config: Dict[str, Any]) -> List[OptimizationConfig]:
        """
        Extract optimization configurations.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            List[OptimizationConfig]: List of optimization configurations.
        """
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
        """
        Extract sensitivity analysis configurations.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            List[SensitivityConfig]: List of sensitivity configurations.
        """
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


# Simple facade used by ExperimentService (compatibility)
class ConfigParser:
    """
    Simple configuration parser facade for ExperimentService compatibility.
    
    Provides aliases to ConfigurationParser methods for backward compatibility.
    """
    
    parse_config = staticmethod(ConfigurationParser.load_file)
    load_config = staticmethod(ConfigurationParser.load_file)


class ConfigurationValidator:
    """
    Configuration validation utilities for legacy compatibility.
    
    Provides validation methods for different configuration types and structures.
    """
    
    @staticmethod
    def validate_batch_structure(config: Dict[str, Any]) -> str:
        """
        Validate batch configuration structure and determine type.
        
        Args:
            config (Dict[str, Any]): Configuration dictionary to validate.
            
        Returns:
            str: Configuration type ('experiment', 'optimization', 'sensitivity').
            
        Raises:
            BatchConfigurationError: If structure is invalid or unrecognized.
        """
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
        """
        Validate optimization configuration.
        
        Args:
            config (OptimizationConfig): Optimization configuration to validate.
            
        Raises:
            OptimizationConfigurationError: If configuration is invalid.
        """
        if not config.target_datasets:
            raise OptimizationConfigurationError("target_datasets list cannot be empty")
        if not isinstance(config.parameters, dict):
            raise OptimizationConfigurationError("parameters must be a mapping")

    @staticmethod
    def validate_sensitivity_config(config: SensitivityConfig) -> None:
        """
        Validate sensitivity analysis configuration.
        
        Args:
            config (SensitivityConfig): Sensitivity configuration to validate.
            
        Raises:
            SensitivityConfigurationError: If configuration is invalid.
        """
        valid_methods = {"morris", "sobol", "fast", "delta"}
        if config.analysis_method not in valid_methods:
            raise SensitivityConfigurationError("Invalid analysis method")
        if not config.target_datasets:
            raise SensitivityConfigurationError("target_datasets list cannot be empty")
