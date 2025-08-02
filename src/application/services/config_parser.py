"""
CSPBench Configuration Parsing Module

Provides modular classes and functions for parsing and validating
batch, optimization, and sensitivity analysis configurations.

This module contains ONLY parameters documented in TEMPLATE.yaml.
"""

from dataclasses import dataclass, field
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
    """Batch metadata matching TEMPLATE.yaml Section 1."""

    name: str
    description: str
    author: str
    version: str
    creation_date: str
    tags: List[str]


@dataclass
class InfrastructureConfig:
    """Infrastructure configuration matching TEMPLATE.yaml Section 2."""

    history: Optional[Dict[str, Any]] = None
    result: Optional[Dict[str, Any]] = None


@dataclass
class DatasetConfig:
    """Dataset configuration matching TEMPLATE.yaml Section 3."""

    id: str
    name: str
    type: str  # "synthetic", "file", "entrez"
    parameters: Dict[str, Any]


@dataclass
class AlgorithmConfig:
    """Algorithm configuration matching TEMPLATE.yaml Section 4."""

    id: str
    name: str
    description: str
    algorithms: List[str]
    algorithm_params: Dict[str, Dict[str, Any]]


@dataclass
class TaskConfig:
    """Task configuration matching TEMPLATE.yaml Section 5."""

    type: str  # "execution", "optimization", "sensitivity"


@dataclass
class ExecutionConfig:
    """Execution configuration matching TEMPLATE.yaml Section 6A."""

    name: str
    datasets: List[str]
    algorithms: List[str]
    repetitions: int


@dataclass
class OptimizationConfig:
    """Optimization configuration matching TEMPLATE.yaml Section 6B."""

    name: str
    study_name: str
    direction: str
    trials: int
    timeout_per_trial: int
    repetitions: int
    datasets: List[str]
    algorithm: str
    parameters: Dict[str, Any]
    optuna_config: Optional[Dict[str, Any]] = None


@dataclass
class SensitivityConfig:
    """Sensitivity configuration matching TEMPLATE.yaml Section 6C."""

    name: str
    method: str  # "morris", "sobol", "fast", "delta"
    datasets: List[str]
    algorithm: str
    samples: int
    repetitions: int
    parameters: Dict[str, Any]
    output_metrics: List[str]
    morris: Optional[Dict[str, Any]] = None
    sobol: Optional[Dict[str, Any]] = None
    fast: Optional[Dict[str, Any]] = None


@dataclass
class ExportConfig:
    """Export configuration matching TEMPLATE.yaml Section 7."""

    enabled: bool = True
    destination: str = "outputs/{session}"
    formats: Optional[Dict[str, bool]] = None
    format_options: Optional[Dict[str, Any]] = None
    directory_structure: Optional[Dict[str, bool]] = None
    include: Optional[List[str]] = None


@dataclass
@dataclass
class PlotsConfig:
    """Plots configuration matching TEMPLATE.yaml Section 8."""

    # Output formats (required field)
    formats: List[str] = field(default_factory=lambda: ["png", "pdf"])

    # Optional fields with defaults
    enabled: bool = True
    # Common plots
    convergence: bool = True
    comparison: bool = True
    boxplots: bool = True
    scatter: bool = True
    heatmap: bool = True
    runtime: bool = True
    success_rate: bool = True
    # Optimization-specific plots
    optimization_history: bool = True
    parameter_importance: bool = True
    parallel_coordinate: bool = True
    # Sensitivity-specific plots
    sensitivity_indices: bool = True
    morris_trajectories: bool = True
    interaction_effects: bool = True

    # Legacy compatibility properties
    @property
    def plot_convergence(self) -> bool:
        return self.convergence

    @property
    def plot_comparison(self) -> bool:
        return self.comparison

    @property
    def plot_boxplots(self) -> bool:
        return self.boxplots

    @property
    def plot_runtime(self) -> bool:
        return self.runtime

    @property
    def plot_scatter(self) -> bool:
        return self.scatter

    @property
    def plot_heatmap(self) -> bool:
        return self.heatmap

    @property
    def plot_success_rate(self) -> bool:
        return self.success_rate

    @property
    def plot_optimization_history(self) -> bool:
        return self.optimization_history

    @property
    def plot_parameter_importance(self) -> bool:
        return self.parameter_importance

    @property
    def plot_parallel_coordinate(self) -> bool:
        return self.parallel_coordinate

    @property
    def plot_sensitivity_indices(self) -> bool:
        return self.sensitivity_indices

    @property
    def plot_morris_trajectories(self) -> bool:
        return self.morris_trajectories

    @property
    def plot_interaction_effects(self) -> bool:
        return self.interaction_effects


@dataclass
class MonitoringConfig:
    """Monitoring configuration matching TEMPLATE.yaml Section 9."""

    enabled: bool = True
    interface: str = "simple"  # "simple", "tui"
    update_interval: int = 3


@dataclass
class ResourcesConfig:
    """Resources configuration matching TEMPLATE.yaml Section 10."""

    cpu: Optional[Dict[str, Any]] = None
    memory: Optional[Dict[str, Any]] = None
    parallel: Optional[Dict[str, Any]] = None
    timeouts: Optional[Dict[str, Any]] = None


@dataclass
class LoggingConfig:
    """Logging configuration matching TEMPLATE.yaml Section 11."""

    level: str = "INFO"  # "DEBUG", "INFO", "WARNING", "ERROR"
    output: Optional[Dict[str, bool]] = None


@dataclass
class SystemConfig:
    """System configuration matching TEMPLATE.yaml Section 12."""

    global_seed: Optional[int] = None
    work_directory: Optional[str] = None
    force_cleanup: bool = False
    checkpointing: Optional[Dict[str, Any]] = None
    error_handling: Optional[Dict[str, Any]] = None
    progress_tracking: Optional[Dict[str, Any]] = None
    environment: Optional[Dict[str, Any]] = None
    reproducibility: Optional[Dict[str, Any]] = None


@dataclass
class BatchConfig:
    """Complete batch configuration combining all sections."""

    metadata: BatchMetadata
    infrastructure: Optional[InfrastructureConfig] = None
    datasets: List[DatasetConfig] = field(default_factory=list)
    algorithms: List[AlgorithmConfig] = field(default_factory=list)
    task: TaskConfig = None
    execution: Optional[Dict[str, Any]] = None
    optimization: Optional[Dict[str, Any]] = None
    sensitivity: Optional[Dict[str, Any]] = None
    export: ExportConfig = field(default_factory=ExportConfig)
    plots: PlotsConfig = field(default_factory=PlotsConfig)
    monitoring: MonitoringConfig = field(default_factory=MonitoringConfig)
    resources: ResourcesConfig = field(default_factory=ResourcesConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    system: SystemConfig = field(default_factory=SystemConfig)


class ConfigParser:
    """
    Configuration parser that handles ONLY parameters from TEMPLATE.yaml.

    This class ensures strict adherence to the standardized template structure.
    """

    @staticmethod
    def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
        """Load and validate a YAML configuration file."""
        try:
            with open(config_path, encoding="utf-8") as file:
                config = yaml.safe_load(file)
                if not config:
                    raise BatchConfigurationError(
                        f"Empty configuration file: {config_path}"
                    )
                return config
        except FileNotFoundError:
            raise BatchConfigurationError(
                f"Configuration file not found: {config_path}"
            )
        except yaml.YAMLError as e:
            raise BatchConfigurationError(f"Invalid YAML syntax: {e}")
        except Exception as e:
            raise BatchConfigurationError(f"Error loading configuration: {e}")

    @staticmethod
    def parse_metadata(config: Dict[str, Any]) -> BatchMetadata:
        """Parse metadata section (Section 1)."""
        metadata_section = config.get("metadata", {})

        required_fields = ["name", "description", "author", "version", "creation_date"]
        for field in required_fields:
            if field not in metadata_section:
                raise BatchConfigurationError(
                    f"Missing required metadata field: {field}"
                )

        return BatchMetadata(
            name=metadata_section["name"],
            description=metadata_section["description"],
            author=metadata_section["author"],
            version=metadata_section["version"],
            creation_date=metadata_section["creation_date"],
            tags=metadata_section.get("tags", []),
        )

    @staticmethod
    def parse_infrastructure(config: Dict[str, Any]) -> Optional[InfrastructureConfig]:
        """Parse infrastructure section (Section 2)."""
        infra_section = config.get("infrastructure")
        if not infra_section:
            return None

        return InfrastructureConfig(
            history=infra_section.get("history"), result=infra_section.get("result")
        )

    @staticmethod
    def parse_datasets(config: Dict[str, Any]) -> List[DatasetConfig]:
        """Parse datasets section (Section 3)."""
        datasets_section = config.get("datasets", [])

        datasets = []
        for dataset_data in datasets_section:
            # Validate required fields
            required_fields = ["id", "name", "type", "parameters"]
            for field in required_fields:
                if field not in dataset_data:
                    raise BatchConfigurationError(
                        f"Missing required dataset field: {field}"
                    )

            # Validate dataset type
            valid_types = ["synthetic", "file", "entrez"]
            if dataset_data["type"] not in valid_types:
                raise BatchConfigurationError(
                    f"Invalid dataset type: {dataset_data['type']}. Must be one of {valid_types}"
                )

            datasets.append(
                DatasetConfig(
                    id=dataset_data["id"],
                    name=dataset_data["name"],
                    type=dataset_data["type"],
                    parameters=dataset_data["parameters"],
                )
            )

        return datasets

    @staticmethod
    def parse_algorithms(config: Dict[str, Any]) -> List[AlgorithmConfig]:
        """Parse algorithms section (Section 4)."""
        algorithms_section = config.get("algorithms", [])

        algorithms = []
        for alg_data in algorithms_section:
            # Validate required fields
            required_fields = ["id", "name", "algorithms", "algorithm_params"]
            for field in required_fields:
                if field not in alg_data:
                    raise BatchConfigurationError(
                        f"Missing required algorithm field: {field}"
                    )

            algorithms.append(
                AlgorithmConfig(
                    id=alg_data["id"],
                    name=alg_data["name"],
                    description=alg_data.get("description", ""),
                    algorithms=alg_data["algorithms"],
                    algorithm_params=alg_data["algorithm_params"],
                )
            )

        return algorithms

    @staticmethod
    def parse_task(config: Dict[str, Any]) -> TaskConfig:
        """Parse task section (Section 5)."""
        task_section = config.get("task", {})

        task_type = task_section.get("type")
        if not task_type:
            raise BatchConfigurationError("Missing required task.type")

        valid_types = ["execution", "optimization", "sensitivity"]
        if task_type not in valid_types:
            raise BatchConfigurationError(
                f"Invalid task type: {task_type}. Must be one of {valid_types}"
            )

        return TaskConfig(type=task_type)

    @staticmethod
    def parse_execution(config: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Parse execution section (Section 6A)."""
        execution_section = config.get("execution")
        if not execution_section:
            return None

        executions = []
        for exec_data in execution_section.get("executions", []):
            # Validate required fields
            required_fields = ["name", "datasets", "algorithms", "repetitions"]
            for field in required_fields:
                if field not in exec_data:
                    raise BatchConfigurationError(
                        f"Missing required execution field: {field}"
                    )

            executions.append(
                ExecutionConfig(
                    name=exec_data["name"],
                    datasets=exec_data["datasets"],
                    algorithms=exec_data["algorithms"],
                    repetitions=exec_data["repetitions"],
                )
            )

        return {"executions": executions}

    @staticmethod
    def parse_optimization(config: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Parse optimization section (Section 6B)."""
        optimization_section = config.get("optimization")
        if not optimization_section:
            return None

        # Parse global optuna defaults
        method = optimization_section.get("method", "optuna")
        if method != "optuna":
            raise OptimizationConfigurationError(
                f"Unsupported optimization method: {method}"
            )

        optuna_defaults = optimization_section.get("optuna_defaults", {})

        # Parse individual optimizations
        optimizations = []
        for opt_data in optimization_section.get("optimizations", []):
            # Validate required fields
            required_fields = [
                "name",
                "study_name",
                "direction",
                "trials",
                "datasets",
                "algorithm",
                "parameters",
            ]
            for field in required_fields:
                if field not in opt_data:
                    raise OptimizationConfigurationError(
                        f"Missing required optimization field: {field}"
                    )

            # Validate direction
            if opt_data["direction"] not in ["minimize", "maximize"]:
                raise OptimizationConfigurationError(
                    f"Invalid direction: {opt_data['direction']}"
                )

            optimizations.append(
                OptimizationConfig(
                    name=opt_data["name"],
                    study_name=opt_data["study_name"],
                    direction=opt_data["direction"],
                    trials=opt_data["trials"],
                    timeout_per_trial=opt_data.get("timeout_per_trial", 300),
                    repetitions=opt_data.get("repetitions", 1),
                    datasets=opt_data["datasets"],
                    algorithm=opt_data["algorithm"],
                    parameters=opt_data["parameters"],
                    optuna_config=opt_data.get("optuna_config"),
                )
            )

        return {
            "method": method,
            "optuna_defaults": optuna_defaults,
            "optimizations": optimizations,
        }

    @staticmethod
    def parse_sensitivity(config: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Parse sensitivity section (Section 6C)."""
        sensitivity_section = config.get("sensitivity")
        if not sensitivity_section:
            return None

        # Parse global salib defaults
        method = sensitivity_section.get("method", "SALib")
        if method != "SALib":
            raise SensitivityConfigurationError(
                f"Unsupported sensitivity method: {method}"
            )

        salib_defaults = sensitivity_section.get("salib_defaults", {})

        # Parse individual analyses
        analyses = []
        for analysis_data in sensitivity_section.get("analyses", []):
            # Validate required fields
            required_fields = [
                "name",
                "method",
                "datasets",
                "algorithm",
                "parameters",
                "output_metrics",
            ]
            for field in required_fields:
                if field not in analysis_data:
                    raise SensitivityConfigurationError(
                        f"Missing required sensitivity field: {field}"
                    )

            # Validate method
            valid_methods = ["morris", "sobol", "fast", "delta"]
            if analysis_data["method"] not in valid_methods:
                raise SensitivityConfigurationError(
                    f"Invalid sensitivity method: {analysis_data['method']}"
                )

            analyses.append(
                SensitivityConfig(
                    name=analysis_data["name"],
                    method=analysis_data["method"],
                    datasets=analysis_data["datasets"],
                    algorithm=analysis_data["algorithm"],
                    samples=analysis_data.get(
                        "samples", salib_defaults.get("samples", 1000)
                    ),
                    repetitions=analysis_data.get("repetitions", 1),
                    parameters=analysis_data["parameters"],
                    output_metrics=analysis_data["output_metrics"],
                    morris=analysis_data.get("morris"),
                    sobol=analysis_data.get("sobol"),
                    fast=analysis_data.get("fast"),
                )
            )

        return {
            "method": method,
            "salib_defaults": salib_defaults,
            "analyses": analyses,
        }

    @staticmethod
    def parse_export(config: Dict[str, Any]) -> ExportConfig:
        """Parse export section - supports both new (output) and old (export) formats."""
        # Try new unified output configuration first
        output_section = config.get("output", {})
        if output_section:
            # Use new unified format
            results_config = output_section.get("results", {})
            structure_config = output_section.get("structure", {})

            return ExportConfig(
                enabled=results_config.get("enabled", True),
                destination=output_section.get("base_directory", "outputs/{session}"),
                formats=results_config.get("formats", {}),
                format_options=results_config.get("format_options", {}),
                directory_structure=structure_config,
                include=(
                    list(results_config.get("content", {}).keys())
                    if results_config.get("content")
                    else []
                ),
            )
        else:
            # Fallback to old export format
            export_section = config.get("export", {})
            return ExportConfig(
                enabled=export_section.get("enabled", True),
                destination=export_section.get("destination", "outputs/{session}"),
                formats=export_section.get("formats"),
                format_options=export_section.get("format_options"),
                directory_structure=export_section.get("directory_structure"),
                include=export_section.get("include"),
            )

    @staticmethod
    def parse_plots(config: Dict[str, Any]) -> PlotsConfig:
        """Parse plots section (Section 8)."""
        plots_section = config.get("plots", {})

        return PlotsConfig(
            enabled=plots_section.get("enabled", True),
            # Common plots
            convergence=plots_section.get("convergence", True),
            comparison=plots_section.get("comparison", True),
            boxplots=plots_section.get("boxplots", True),
            scatter=plots_section.get("scatter", True),
            heatmap=plots_section.get("heatmap", True),
            runtime=plots_section.get("runtime", True),
            success_rate=plots_section.get("success_rate", True),
            # Optimization-specific
            optimization_history=plots_section.get("optimization_history", True),
            parameter_importance=plots_section.get("parameter_importance", True),
            parallel_coordinate=plots_section.get("parallel_coordinate", True),
            # Sensitivity-specific
            sensitivity_indices=plots_section.get("sensitivity_indices", True),
            morris_trajectories=plots_section.get("morris_trajectories", True),
            interaction_effects=plots_section.get("interaction_effects", True),
            # Output formats
            formats=plots_section.get("formats", ["png", "pdf"]),
        )

    @staticmethod
    def parse_monitoring(config: Dict[str, Any]) -> MonitoringConfig:
        """Parse monitoring section (Section 9)."""
        monitoring_section = config.get("monitoring", {})

        interface = monitoring_section.get("interface", "simple")
        if interface not in ["simple", "tui"]:
            raise BatchConfigurationError(f"Invalid monitoring interface: {interface}")

        return MonitoringConfig(
            enabled=monitoring_section.get("enabled", True),
            interface=interface,
            update_interval=monitoring_section.get("update_interval", 3),
        )

    @staticmethod
    def parse_resources(config: Dict[str, Any]) -> ResourcesConfig:
        """Parse resources section (Section 10)."""
        resources_section = config.get("resources", {})

        return ResourcesConfig(
            cpu=resources_section.get("cpu"),
            memory=resources_section.get("memory"),
            parallel=resources_section.get("parallel"),
            timeouts=resources_section.get("timeouts"),
        )

    @staticmethod
    def parse_logging(config: Dict[str, Any]) -> LoggingConfig:
        """Parse logging section (Section 11)."""
        logging_section = config.get("logging", {})

        level = logging_section.get("level", "INFO")
        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR"]
        if level not in valid_levels:
            raise BatchConfigurationError(
                f"Invalid logging level: {level}. Must be one of {valid_levels}"
            )

        return LoggingConfig(level=level, output=logging_section.get("output"))

    @staticmethod
    def parse_system(config: Dict[str, Any]) -> SystemConfig:
        """Parse system section (Section 12)."""
        system_section = config.get("system", {})

        return SystemConfig(
            global_seed=system_section.get("global_seed"),
            work_directory=system_section.get("work_directory"),
            force_cleanup=system_section.get("force_cleanup", False),
            checkpointing=system_section.get("checkpointing"),
            error_handling=system_section.get("error_handling"),
            progress_tracking=system_section.get("progress_tracking"),
            environment=system_section.get("environment"),
            reproducibility=system_section.get("reproducibility"),
        )

    @classmethod
    def parse_config(cls, config_path: Union[str, Path]) -> BatchConfig:
        """
        Parse a complete batch configuration file.

        This method validates that ONLY parameters from TEMPLATE.yaml are used.
        """
        config = cls.load_config(config_path)

        # Parse all sections according to TEMPLATE.yaml structure
        metadata = cls.parse_metadata(config)
        infrastructure = cls.parse_infrastructure(config)
        datasets = cls.parse_datasets(config)
        algorithms = cls.parse_algorithms(config)
        task = cls.parse_task(config)
        execution = cls.parse_execution(config)
        optimization = cls.parse_optimization(config)
        sensitivity = cls.parse_sensitivity(config)
        export = cls.parse_export(config)
        plots = cls.parse_plots(config)
        monitoring = cls.parse_monitoring(config)
        resources = cls.parse_resources(config)
        logging = cls.parse_logging(config)
        system = cls.parse_system(config)

        # Validate task type consistency
        task_type = task.type
        if task_type == "execution" and not execution:
            raise BatchConfigurationError(
                "Task type 'execution' requires 'execution' section"
            )
        if task_type == "optimization" and not optimization:
            raise BatchConfigurationError(
                "Task type 'optimization' requires 'optimization' section"
            )
        if task_type == "sensitivity" and not sensitivity:
            raise BatchConfigurationError(
                "Task type 'sensitivity' requires 'sensitivity' section"
            )

        return BatchConfig(
            metadata=metadata,
            infrastructure=infrastructure,
            datasets=datasets,
            algorithms=algorithms,
            task=task,
            execution=execution,
            optimization=optimization,
            sensitivity=sensitivity,
            export=export,
            plots=plots,
            monitoring=monitoring,
            resources=resources,
            logging=logging,
            system=system,
        )


# Legacy function for backward compatibility
def load_batch_config(config_path: Union[str, Path]) -> BatchConfig:
    """Load and parse a batch configuration file (legacy function)."""
    return ConfigParser.parse_config(config_path)


def parse_batch_config(config_path: Union[str, Path]) -> BatchConfig:
    """Parse a batch configuration file."""
    return ConfigParser.parse_config(config_path)


def save_batch_config(config: BatchConfig, output_path: Union[str, Path]) -> None:
    """Save a batch configuration to a YAML file."""
    # Convert dataclass to dict for YAML serialization
    config_dict = {
        "metadata": {
            "name": config.metadata.name,
            "description": config.metadata.description,
            "author": config.metadata.author,
            "version": config.metadata.version,
            "creation_date": config.metadata.creation_date,
            "tags": config.metadata.tags,
        }
    }

    # Add optional sections only if they exist
    if config.infrastructure:
        config_dict["infrastructure"] = {
            "history": config.infrastructure.history,
            "result": config.infrastructure.result,
        }

    # Add datasets
    if config.datasets:
        config_dict["datasets"] = [
            {"id": ds.id, "name": ds.name, "type": ds.type, "parameters": ds.parameters}
            for ds in config.datasets
        ]

    # Add algorithms
    if config.algorithms:
        config_dict["algorithms"] = [
            {
                "id": alg.id,
                "name": alg.name,
                "description": alg.description,
                "algorithms": alg.algorithms,
                "algorithm_params": alg.algorithm_params,
            }
            for alg in config.algorithms
        ]

    # Add task
    config_dict["task"] = {"type": config.task.type}

    # Add task-specific configurations
    if config.execution:
        config_dict["execution"] = config.execution
    if config.optimization:
        config_dict["optimization"] = config.optimization
    if config.sensitivity:
        config_dict["sensitivity"] = config.sensitivity

    # Add all other sections with their current values
    config_dict["export"] = {
        "enabled": config.export.enabled,
        "destination": config.export.destination,
        "formats": config.export.formats,
        "format_options": config.export.format_options,
        "directory_structure": config.export.directory_structure,
        "include": config.export.include,
    }

    config_dict["plots"] = {
        "enabled": config.plots.enabled,
        "convergence": config.plots.convergence,
        "comparison": config.plots.comparison,
        "boxplots": config.plots.boxplots,
        "scatter": config.plots.scatter,
        "heatmap": config.plots.heatmap,
        "runtime": config.plots.runtime,
        "success_rate": config.plots.success_rate,
        "optimization_history": config.plots.optimization_history,
        "parameter_importance": config.plots.parameter_importance,
        "parallel_coordinate": config.plots.parallel_coordinate,
        "sensitivity_indices": config.plots.sensitivity_indices,
        "morris_trajectories": config.plots.morris_trajectories,
        "interaction_effects": config.plots.interaction_effects,
        "formats": config.plots.formats,
    }

    config_dict["monitoring"] = {
        "enabled": config.monitoring.enabled,
        "interface": config.monitoring.interface,
        "update_interval": config.monitoring.update_interval,
    }

    config_dict["resources"] = {
        "cpu": config.resources.cpu,
        "memory": config.resources.memory,
        "parallel": config.resources.parallel,
        "timeouts": config.resources.timeouts,
    }

    config_dict["logging"] = {
        "level": config.logging.level,
        "output": config.logging.output,
    }
    config_dict["system"] = {"reproducibility": config.system.reproducibility}

    # Write to YAML file
    with open(output_path, "w", encoding="utf-8") as file:
        yaml.dump(
            config_dict, file, default_flow_style=False, allow_unicode=True, indent=2
        )

    @staticmethod
    def parse_optimization_configs(
        config: Dict[str, Any],
    ) -> List["OptimizationConfig"]:
        """Extract optimization configurations for legacy compatibility."""
        optimization_section = config.get("optimization", {})
        if not optimization_section:
            raise OptimizationConfigurationError("optimization section not found")

        # Support new structure (optimizations) and legacy (direct fields)
        optimizations_list = optimization_section.get("optimizations")

        if optimizations_list is not None:
            # New structure: list of optimizations
            return [
                ConfigParser._parse_single_optimization(opt_config)
                for opt_config in optimizations_list
            ]
        else:
            # Legacy structure: single configuration
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

            return [ConfigParser._parse_single_optimization(legacy_config)]

    @staticmethod
    def _parse_single_optimization(opt_config: Dict[str, Any]) -> "OptimizationConfig":
        """Parse a single optimization configuration."""
        required_fields = ["algorithm", "parameters"]
        for field in required_fields:
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
    def parse_sensitivity_configs(config: Dict[str, Any]) -> List["SensitivityConfig"]:
        """Extract sensitivity analysis configurations for legacy compatibility."""
        sensitivity_section = config.get("sensitivity", {})
        if not sensitivity_section:
            raise SensitivityConfigurationError("sensitivity section not found")

        # Extract global configurations
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
            return [
                ConfigParser._parse_single_sensitivity(
                    sens_config, global_samples, global_repetitions
                )
                for sens_config in analyses_list
            ]
        else:
            # Legacy structure: single configuration
            legacy_config = sensitivity_section.copy()

            if "target_dataset" in legacy_config:
                legacy_config["datasets"] = [legacy_config.pop("target_dataset")]
            if "target_datasets" in legacy_config:
                legacy_config["datasets"] = legacy_config.pop("target_datasets")

            if "nome" not in legacy_config and "name" not in legacy_config:
                legacy_config["name"] = (
                    f"Analysis {legacy_config.get('target_algorithm', legacy_config.get('algorithm', 'Unknown'))}"
                )

            return [
                ConfigParser._parse_single_sensitivity(
                    legacy_config, global_samples, global_repetitions
                )
            ]

    @staticmethod
    def _parse_single_sensitivity(
        sens_config: Dict[str, Any],
        global_samples: int = 1000,
        global_repetitions: int = 3,
    ) -> "SensitivityConfig":
        """Parse a single sensitivity analysis configuration."""
        required_fields = ["algorithm", "parameters"]
        for field in required_fields:
            if field not in sens_config and f"target_{field}" not in sens_config:
                raise SensitivityConfigurationError(
                    f"Required field missing: {field} or target_{field}"
                )

        # Collect method-specific configuration
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
        """Parse batch resource configurations for legacy compatibility."""
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
