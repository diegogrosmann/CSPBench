"""
CSPBench Configuration Parsing Module

Provides modular classes and functions for parsing and validating
batch, optimization, and sensitivity analysis configurations.
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
    """Batch metadata with standardized field names."""

    name: str  # Standardized English field names
    description: str
    author: str
    version: str
    creation_date: str
    tags: List[str]
    timeout_global: int  # Keep for backward compatibility


@dataclass
class OptimizationConfig:
    """Configuration for a specific optimization."""

    name: str  # Standardized English field names
    study_name: str
    direction: str
    trials: int  # Renamed from n_trials
    timeout_per_trial: int
    repetitions: int  # Added repetitions support
    target_datasets: List[str]
    target_algorithm: str
    parameters: Dict[str, Any]
    optuna_config: Optional[Dict[str, Any]] = None


@dataclass
class SensitivityConfig:
    """Configuration for a specific sensitivity analysis."""

    name: str  # Standardized English field names
    method: str  # Renamed from analysis_method
    target_datasets: List[str]
    target_algorithm: str
    samples: int  # Renamed from n_samples
    repetitions: int  # Renamed from repetitions_per_sample
    parameters: Dict[str, Any]
    output_metrics: List[str]
    method_config: Optional[Dict[str, Any]] = None


@dataclass
class InfrastructureConfig:
    """Infrastructure configuration."""

    history: Optional[Dict[str, Any]] = None
    result: Optional[Dict[str, Any]] = None


@dataclass
class ExportConfig:
    """Export configuration."""

    enabled: bool = True
    destination: str = "outputs/{session_folder}"
    formats: Optional[Dict[str, Any]] = None
    csv_config: Optional[Dict[str, Any]] = None
    json_config: Optional[Dict[str, Any]] = None
    directory_structure: Optional[Dict[str, Any]] = None
    include: Optional[List[str]] = None


@dataclass
class PlotsConfig:
    """Visualization/plots configuration."""

    enabled: bool = True
    plot_convergence: bool = True
    plot_comparison: bool = True
    plot_boxplots: bool = True
    plot_scatter: bool = True
    plot_heatmap: bool = True
    plot_runtime: bool = True
    plot_success_rate: bool = True
    plot_optimization_history: bool = True
    plot_parameter_importance: bool = True
    plot_parallel_coordinate: bool = True
    plot_sensitivity_indices: bool = True
    plot_morris_trajectories: bool = True
    plot_interaction_effects: bool = True
    style: str = "seaborn-v0_8"
    figure_size: List[int] = None
    dpi: int = 300
    color_palette: str = "Set2"
    formats: List[str] = None
    font_size: int = 12
    title_size: int = 16
    label_size: int = 14
    legend_size: int = 10
    tight_layout: bool = True


@dataclass
class MonitoringConfig:
    """Monitoring configuration."""

    enabled: bool = True
    interface: str = "simple"
    update_interval: int = 3


@dataclass
class LoggingConfig:
    """Logging configuration."""

    level: str = "INFO"
    output: Optional[Dict[str, Any]] = None
    file_config: Optional[Dict[str, Any]] = None
    format: Optional[Dict[str, Any]] = None
    components: Optional[Dict[str, Any]] = None
    filters: Optional[Dict[str, Any]] = None
    structured: Optional[Dict[str, Any]] = None


@dataclass
class SystemConfig:
    """System configuration."""

    reproducibility: Optional[Dict[str, Any]] = None
    checkpointing: Optional[Dict[str, Any]] = None
    progress: Optional[Dict[str, Any]] = None
    early_stopping: Optional[Dict[str, Any]] = None
    error_handling: Optional[Dict[str, Any]] = None
    cleanup: Optional[Dict[str, Any]] = None


class ConfigurationParser:
    """Modular parser for CSPBench configurations."""

    @staticmethod
    def load_file(file_path: Union[str, Path]) -> Dict[str, Any]:
        """Load YAML or JSON configuration file."""
        file_path = Path(file_path)

        if not file_path.exists():
            raise BatchConfigurationError(f"File not found: {file_path}")

        with open(file_path, encoding="utf-8") as f:
            if file_path.suffix.lower() in [".yaml", ".yml"]:
                config = yaml.safe_load(f)
            elif file_path.suffix.lower() == ".json":
                config = json.load(f)
            else:
                raise BatchConfigurationError(f"Unsupported format: {file_path.suffix}")

        if not isinstance(config, dict):
            raise BatchConfigurationError("Configuration must be an object/dictionary")

        return config

    @staticmethod
    def parse_metadata(config: Dict[str, Any]) -> BatchMetadata:
        """Extract batch metadata with support for new standardized format."""
        # Support new metadata format (English) and legacy formats (Portuguese)
        metadata_section = config.get(
            "metadata", config.get("metadados", config.get("batch_info", {}))
        )

        if not metadata_section:
            raise BatchConfigurationError("metadata section not found")

        return BatchMetadata(
            name=metadata_section.get(
                "name", metadata_section.get("nome", "untitled batch")
            ),
            description=metadata_section.get(
                "description", metadata_section.get("descricao", "")
            ),
            author=metadata_section.get("author", metadata_section.get("autor", "")),
            version=metadata_section.get(
                "version", metadata_section.get("versao", "1.0")
            ),
            creation_date=metadata_section.get(
                "creation_date", metadata_section.get("data_criacao", "")
            ),
            tags=metadata_section.get("tags", []),
            timeout_global=metadata_section.get("timeout_global", 3600),
        )

    @staticmethod
    def parse_optimization_configs(config: Dict[str, Any]) -> List[OptimizationConfig]:
        """Extract optimization configurations with support for new standardized format."""
        optimization_section = config.get("optimization", {})

        if not optimization_section:
            raise OptimizationConfigurationError("optimization section not found")

        # Support new structure (optimizations) and legacy (direct fields)
        optimizations_list = optimization_section.get("optimizations")

        if optimizations_list is not None:
            # New structure: list of optimizations
            return [
                ConfigurationParser._parse_single_optimization(opt_config)
                for opt_config in optimizations_list
            ]
        else:
            # Legacy structure: single configuration
            # Convert to list of one optimization
            legacy_config = optimization_section.copy()

            # Ensure required fields for legacy structure
            if "target_dataset" in legacy_config:
                legacy_config["datasets"] = [legacy_config.pop("target_dataset")]
            if "target_datasets" in legacy_config:
                legacy_config["datasets"] = legacy_config.pop("target_datasets")

            if "nome" not in legacy_config and "name" not in legacy_config:
                legacy_config["name"] = (
                    f"Optimization {legacy_config.get('target_algorithm', legacy_config.get('algorithm', 'Unknown'))}"
                )

            return [ConfigurationParser._parse_single_optimization(legacy_config)]

    @staticmethod
    def _parse_single_optimization(opt_config: Dict[str, Any]) -> OptimizationConfig:
        """Parse a single optimization configuration with new standardized format."""
        # Verificar se já é um OptimizationConfig
        if (
            hasattr(opt_config, "__class__")
            and "OptimizationConfig" in opt_config.__class__.__name__
        ):
            return opt_config

        required_fields = ["algorithm", "parameters"]  # Updated field names
        for field in required_fields:
            # Support both new and legacy field names
            if field not in opt_config and f"target_{field}" not in opt_config:
                raise OptimizationConfigurationError(
                    f"Required field missing: {field} or target_{field}"
                )

        return OptimizationConfig(
            name=opt_config.get("name", opt_config.get("nome", "Unnamed optimization")),
            study_name=opt_config.get("study_name", "default_study"),
            direction=opt_config.get("direction", "minimize"),
            trials=opt_config.get("trials", opt_config.get("n_trials", 100)),
            timeout_per_trial=opt_config.get("timeout_per_trial", 300),
            repetitions=opt_config.get("repetitions", 1),
            target_datasets=opt_config.get(
                "datasets", opt_config.get("target_datasets", [])
            ),
            target_algorithm=opt_config.get(
                "algorithm", opt_config.get("target_algorithm", "")
            ),
            parameters=opt_config["parameters"],
            optuna_config=opt_config.get("optuna_config"),
        )

    @staticmethod
    def parse_sensitivity_configs(config: Dict[str, Any]) -> List[SensitivityConfig]:
        """Extract sensitivity analysis configurations with support for new standardized format."""
        sensitivity_section = config.get("sensitivity", {})

        if not sensitivity_section:
            raise SensitivityConfigurationError("sensitivity section not found")

        # Extract global configurations with new format support
        global_config = sensitivity_section.get(
            "salib_defaults", sensitivity_section.get("global_salib_config", {})
        )
        global_samples = global_config.get(
            "samples", global_config.get("n_samples", 1000)
        )
        global_repetitions = global_config.get(
            "repetitions", global_config.get("repetitions_per_sample", 3)
        )

        # Support new structure (analyses) and legacy (direct fields)
        analyses_list = sensitivity_section.get("analyses")

        if analyses_list is not None:
            # New structure: list of analyses
            return [
                ConfigurationParser._parse_single_sensitivity(
                    sens_config, global_samples, global_repetitions
                )
                for sens_config in analyses_list
            ]
        else:
            # Legacy structure: single configuration
            legacy_config = sensitivity_section.copy()

            # Ensure required fields for legacy structure
            if "target_dataset" in legacy_config:
                legacy_config["datasets"] = [legacy_config.pop("target_dataset")]
            if "target_datasets" in legacy_config:
                legacy_config["datasets"] = legacy_config.pop("target_datasets")

            if "nome" not in legacy_config and "name" not in legacy_config:
                legacy_config["name"] = (
                    f"Analysis {legacy_config.get('target_algorithm', legacy_config.get('algorithm', 'Unknown'))}"
                )

            return [
                ConfigurationParser._parse_single_sensitivity(
                    legacy_config, global_samples, global_repetitions
                )
            ]

    @staticmethod
    def _parse_single_sensitivity(
        sens_config: Dict[str, Any],
        global_samples: int = 1000,
        global_repetitions: int = 3,
    ) -> SensitivityConfig:
        """Parse a single sensitivity analysis configuration with new standardized format."""
        # Verificar se já é um SensitivityConfig
        if (
            hasattr(sens_config, "__class__")
            and "SensitivityConfig" in sens_config.__class__.__name__
        ):
            return sens_config

        required_fields = ["algorithm", "parameters"]  # Updated field names
        for field in required_fields:
            # Support both new and legacy field names
            if field not in sens_config and f"target_{field}" not in sens_config:
                raise SensitivityConfigurationError(
                    f"Required field missing: {field} or target_{field}"
                )

        # Collect method-specific configuration with new format support
        analysis_method = sens_config.get(
            "method", sens_config.get("analysis_method", "morris")
        )
        method_config = None

        if analysis_method == "morris" and "morris" in sens_config:
            method_config = sens_config["morris"]
        elif analysis_method == "sobol" and "sobol" in sens_config:
            method_config = sens_config["sobol"]
        elif analysis_method == "fast" and "fast" in sens_config:
            method_config = sens_config["fast"]

        return SensitivityConfig(
            name=sens_config.get("name", sens_config.get("nome", "Unnamed analysis")),
            method=analysis_method,
            target_datasets=sens_config.get(
                "datasets", sens_config.get("target_datasets", [])
            ),
            target_algorithm=sens_config.get(
                "algorithm", sens_config.get("target_algorithm", "")
            ),
            samples=sens_config.get(
                "samples", sens_config.get("n_samples", global_samples)
            ),
            repetitions=sens_config.get(
                "repetitions",
                sens_config.get("repetitions_per_sample", global_repetitions),
            ),
            parameters=sens_config["parameters"],
            output_metrics=sens_config.get(
                "output_metrics", ["distance", "execution_time"]
            ),
            method_config=method_config,
        )

    @staticmethod
    def parse_resources_config(config: Dict[str, Any]) -> Dict[str, Any]:
        """Parse batch resource configurations with support for new standardized format."""
        resources_section = config.get("resources", {})
        parallel_config = resources_section.get("parallel", {})
        timeouts_config = resources_section.get("timeouts", {})

        # Support new timeout structure
        timeout_per_algorithm = timeouts_config.get("timeout_per_algorithm", 3600)
        timeout_total_batch = timeouts_config.get("timeout_total_batch", 86400)

        return {
            "enabled": parallel_config.get("enabled", True),
            "max_workers": parallel_config.get("max_workers", 4),
            "internal_jobs": parallel_config.get("internal_jobs", 4),
            "backend": parallel_config.get("backend", "multiprocessing"),
            "timeouts": {
                "timeout_per_algorithm": timeout_per_algorithm,
                "timeout_total_batch": timeout_total_batch,
                # Legacy support
                "global_timeout": timeouts_config.get(
                    "global_timeout", timeout_per_algorithm
                ),
                "per_algorithm_run": timeout_per_algorithm,
                "total_batch": timeout_total_batch,
            },
        }

    @staticmethod
    def parse_infrastructure_config(config: Dict[str, Any]) -> InfrastructureConfig:
        """Parse infrastructure configuration."""
        infrastructure_section = config.get("infrastructure", {})

        return InfrastructureConfig(
            history=infrastructure_section.get("history"),
            result=infrastructure_section.get("result"),
        )

    @staticmethod
    def parse_export_config(config: Dict[str, Any]) -> ExportConfig:
        """Parse export configuration."""
        export_section = config.get("export", {})

        return ExportConfig(
            enabled=export_section.get("enabled", True),
            destination=export_section.get("destination", "outputs/{session_folder}"),
            formats=export_section.get("formats"),
            csv_config=export_section.get("csv_config"),
            json_config=export_section.get("json_config"),
            directory_structure=export_section.get("directory_structure"),
            include=export_section.get("include"),
        )

    @staticmethod
    def parse_plots_config(config: Dict[str, Any]) -> PlotsConfig:
        """Parse plots configuration."""
        plots_section = config.get("plots", {})

        return PlotsConfig(
            enabled=plots_section.get("enabled", True),
            plot_convergence=plots_section.get("plot_convergence", True),
            plot_comparison=plots_section.get("plot_comparison", True),
            plot_boxplots=plots_section.get("plot_boxplots", True),
            plot_scatter=plots_section.get("plot_scatter", True),
            plot_heatmap=plots_section.get("plot_heatmap", True),
            plot_runtime=plots_section.get("plot_runtime", True),
            plot_success_rate=plots_section.get("plot_success_rate", True),
            plot_optimization_history=plots_section.get(
                "plot_optimization_history", True
            ),
            plot_parameter_importance=plots_section.get(
                "plot_parameter_importance", True
            ),
            plot_parallel_coordinate=plots_section.get(
                "plot_parallel_coordinate", True
            ),
            plot_sensitivity_indices=plots_section.get(
                "plot_sensitivity_indices", True
            ),
            plot_morris_trajectories=plots_section.get(
                "plot_morris_trajectories", True
            ),
            plot_interaction_effects=plots_section.get(
                "plot_interaction_effects", True
            ),
            style=plots_section.get("style", "seaborn-v0_8"),
            figure_size=plots_section.get("figure_size", [12, 8]),
            dpi=plots_section.get("dpi", 300),
            color_palette=plots_section.get("color_palette", "Set2"),
            formats=plots_section.get("formats", ["png"]),
            font_size=plots_section.get("font_size", 12),
            title_size=plots_section.get("title_size", 16),
            label_size=plots_section.get("label_size", 14),
            legend_size=plots_section.get("legend_size", 10),
            tight_layout=plots_section.get("tight_layout", True),
        )

    @staticmethod
    def parse_monitoring_config(config: Dict[str, Any]) -> MonitoringConfig:
        """Parse monitoring configuration."""
        monitoring_section = config.get("monitoring", {})

        return MonitoringConfig(
            enabled=monitoring_section.get("enabled", True),
            interface=monitoring_section.get("interface", "simple"),
            update_interval=monitoring_section.get("update_interval", 3),
        )

    @staticmethod
    def parse_logging_config(config: Dict[str, Any]) -> LoggingConfig:
        """Parse logging configuration."""
        logging_section = config.get("logging", {})

        return LoggingConfig(
            level=logging_section.get("level", "INFO"),
            output=logging_section.get("output"),
            file_config=logging_section.get("file_config"),
            format=logging_section.get("format"),
            components=logging_section.get("components"),
            filters=logging_section.get("filters"),
            structured=logging_section.get("structured"),
        )

    @staticmethod
    def parse_system_config(config: Dict[str, Any]) -> SystemConfig:
        """Parse system configuration."""
        system_section = config.get("system", {})

        return SystemConfig(
            reproducibility=system_section.get("reproducibility"),
            checkpointing=system_section.get("checkpointing"),
            progress=system_section.get("progress"),
            early_stopping=system_section.get("early_stopping"),
            error_handling=system_section.get("error_handling"),
            cleanup=system_section.get("cleanup"),
        )

    @staticmethod
    def parse_infrastructure_config(config: Dict[str, Any]) -> InfrastructureConfig:
        """Parse infrastructure configuration."""
        infrastructure_section = config.get("infrastructure", {})

        return InfrastructureConfig(
            history=infrastructure_section.get("history"),
            result=infrastructure_section.get("result"),
        )

    @staticmethod
    def parse_export_config(config: Dict[str, Any]) -> ExportConfig:
        """Parse export configuration."""
        export_section = config.get("export", {})

        return ExportConfig(
            enabled=export_section.get("enabled", True),
            destination=export_section.get("destination", "outputs/{session_folder}"),
            formats=export_section.get("formats"),
            csv_config=export_section.get("csv_config"),
            json_config=export_section.get("json_config"),
            directory_structure=export_section.get("directory_structure"),
            include=export_section.get("include"),
        )

    @staticmethod
    def parse_plots_config(config: Dict[str, Any]) -> PlotsConfig:
        """Parse plots/visualization configuration."""
        plots_section = config.get("plots", {})

        return PlotsConfig(
            enabled=plots_section.get("enabled", True),
            plot_convergence=plots_section.get("plot_convergence", True),
            plot_comparison=plots_section.get("plot_comparison", True),
            plot_boxplots=plots_section.get("plot_boxplots", True),
            plot_scatter=plots_section.get("plot_scatter", True),
            plot_heatmap=plots_section.get("plot_heatmap", True),
            plot_runtime=plots_section.get("plot_runtime", True),
            plot_success_rate=plots_section.get("plot_success_rate", True),
            plot_optimization_history=plots_section.get(
                "plot_optimization_history", True
            ),
            plot_parameter_importance=plots_section.get(
                "plot_parameter_importance", True
            ),
            plot_parallel_coordinate=plots_section.get(
                "plot_parallel_coordinate", True
            ),
            plot_sensitivity_indices=plots_section.get(
                "plot_sensitivity_indices", True
            ),
            plot_morris_trajectories=plots_section.get(
                "plot_morris_trajectories", True
            ),
            plot_interaction_effects=plots_section.get(
                "plot_interaction_effects", True
            ),
            style=plots_section.get("style", "seaborn-v0_8"),
            figure_size=plots_section.get("figure_size", [12, 8]),
            dpi=plots_section.get("dpi", 300),
            color_palette=plots_section.get("color_palette", "Set2"),
            formats=plots_section.get("formats", ["png"]),
            font_size=plots_section.get("font_size", 12),
            title_size=plots_section.get("title_size", 16),
            label_size=plots_section.get("label_size", 14),
            legend_size=plots_section.get("legend_size", 10),
            tight_layout=plots_section.get("tight_layout", True),
        )

    @staticmethod
    def parse_monitoring_config(config: Dict[str, Any]) -> MonitoringConfig:
        """Parse monitoring configuration."""
        monitoring_section = config.get("monitoring", {})

        return MonitoringConfig(
            enabled=monitoring_section.get("enabled", True),
            interface=monitoring_section.get("interface", "simple"),
            update_interval=monitoring_section.get("update_interval", 3),
        )

    @staticmethod
    def parse_logging_config(config: Dict[str, Any]) -> LoggingConfig:
        """Parse logging configuration."""
        logging_section = config.get("logging", {})

        return LoggingConfig(
            level=logging_section.get("level", "INFO"),
            output=logging_section.get("output"),
            file_config=logging_section.get("file_config"),
            format=logging_section.get("format"),
            components=logging_section.get("components"),
            filters=logging_section.get("filters"),
            structured=logging_section.get("structured"),
        )

    @staticmethod
    def parse_system_config(config: Dict[str, Any]) -> SystemConfig:
        """Parse system configuration."""
        system_section = config.get("system", {})

        return SystemConfig(
            reproducibility=system_section.get("reproducibility"),
            checkpointing=system_section.get("checkpointing"),
            progress=system_section.get("progress"),
            early_stopping=system_section.get("early_stopping"),
            error_handling=system_section.get("error_handling"),
            cleanup=system_section.get("cleanup"),
        )


class ConfigurationValidator:
    """Modular validator for CSPBench configurations."""

    @staticmethod
    def validate_batch_structure(config: Dict[str, Any]) -> str:
        """Validate and return batch structure type."""
        # Detect structure type
        task_type = config.get("task", {}).get("type")

        if task_type in ["execution", "optimization", "sensitivity"]:
            # New unified structure
            required_sections = ["metadados", "datasets", "algorithms", "task"]
            missing = [
                section for section in required_sections if section not in config
            ]

            if missing:
                # Check if it's batch_info instead of metadados (backward compatibility)
                if "batch_info" in config and "metadados" in missing:
                    missing.remove("metadados")

            if missing:
                raise BatchConfigurationError(f"Required sections missing: {missing}")

            return task_type

        elif "experiments" in config:
            # Legacy processing structure
            return "execution"

        else:
            raise BatchConfigurationError(
                "Configuration structure not recognized. "
                "Check if there's a 'task' section with valid 'type' or 'experiments' section"
            )

    @staticmethod
    def validate_optimization_config(opt_config: OptimizationConfig) -> None:
        """Validate optimization configuration."""
        if not opt_config.target_datasets:
            raise OptimizationConfigurationError("target_datasets list cannot be empty")

        if not opt_config.target_algorithm:
            raise OptimizationConfigurationError("target_algorithm field is required")

        if not opt_config.parameters:
            raise OptimizationConfigurationError("parameters section cannot be empty")

        # Validate parameters
        for param_name, param_config in opt_config.parameters.items():
            if not isinstance(param_config, dict):
                raise OptimizationConfigurationError(
                    f"Parameter configuration {param_name} must be an object"
                )

            param_type = param_config.get("type")
            if param_type not in ["int", "uniform", "categorical"]:
                raise OptimizationConfigurationError(
                    f"Invalid parameter type for {param_name}: {param_type}"
                )

    @staticmethod
    def validate_sensitivity_config(sens_config: SensitivityConfig) -> None:
        """Validate sensitivity analysis configuration."""
        if not sens_config.target_datasets:
            raise SensitivityConfigurationError("target_datasets list cannot be empty")

        if not sens_config.target_algorithm:
            raise SensitivityConfigurationError("target_algorithm field is required")

        if not sens_config.parameters:
            raise SensitivityConfigurationError("parameters section cannot be empty")

        # Validate analysis method
        valid_methods = ["morris", "sobol", "fast", "delta"]
        if sens_config.analysis_method not in valid_methods:
            raise SensitivityConfigurationError(
                f"Invalid analysis method: {sens_config.analysis_method}. "
                f"Valid methods: {valid_methods}"
            )

        # Validate parameters
        for param_name, param_config in sens_config.parameters.items():
            if not isinstance(param_config, dict):
                raise SensitivityConfigurationError(
                    f"Parameter configuration {param_name} must be an object"
                )

            param_type = param_config.get("type")
            if param_type not in ["integer", "float", "categorical"]:
                raise SensitivityConfigurationError(
                    f"Invalid parameter type for {param_name}: {param_type}"
                )
