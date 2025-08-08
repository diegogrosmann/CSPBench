"""
CSPBench Batch Configuration Models

Domain models for batch configuration structures and related entities.
These models match the TEMPLATE.yaml sections and provide type safety for batch operations.
"""

from __future__ import annotations  # Enable forward references
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, TYPE_CHECKING

from src.domain.algorithms import CSPAlgorithm

from .enums import BatchType, DatasetType, ExportFormat, OptimizationMethod, SensitivityMethod


@dataclass
class BatchMetadata:
    """
    Batch metadata matching TEMPLATE.yaml Section 1.
    
    Contains essential information for batch identification and traceability.
    """
    name: str
    description: Optional[str] = None
    author: Optional[str] = None
    version: Optional[str] = None
    creation_date: Optional[str] = None
    tags: Optional[List[str]] = None


@dataclass
class DatasetConfig:
    """
    Dataset configuration matching TEMPLATE.yaml Section 2.
    
    Defines datasets available for use in batch executions.
    """
    id: str
    name: str
    type: DatasetType
    parameters: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if not self.parameters:
            raise ValueError("parameters cannot be empty")
        
        # Validações específicas por tipo de dataset
        if self.type == DatasetType.SYNTHETIC:
            required_params = {"n_strings", "string_length", "alphabet"}
            missing = required_params - set(self.parameters.keys())
            if missing:
                raise ValueError(f"Synthetic dataset missing required parameters: {missing}")
            
            # Validar tipos e valores
            if not isinstance(self.parameters.get("n_strings"), int) or self.parameters["n_strings"] <= 0:
                raise ValueError("n_strings must be positive integer")
            if not isinstance(self.parameters.get("string_length"), int) or self.parameters["string_length"] <= 0:
                raise ValueError("string_length must be positive integer")
            if not isinstance(self.parameters.get("alphabet"), str) or len(self.parameters["alphabet"]) == 0:
                raise ValueError("alphabet must be non-empty string")
        
        elif self.type == DatasetType.FILE:
            if "filepath" not in self.parameters:
                raise ValueError("File dataset requires 'filepath' parameter")
            if not isinstance(self.parameters["filepath"], str):
                raise ValueError("filepath must be string")
        
        elif self.type == DatasetType.ENTREZ:
            required_params = {"database", "query"}
            missing = required_params - set(self.parameters.keys())
            if missing:
                raise ValueError(f"Entrez dataset missing required parameters: {missing}")

@dataclass
class AlgorithmParamsConfig:
    """
    Single algorithm configuration with its specific parameters.
    
    Represents one algorithm and its parameter settings.
    """
    name: str  # Nome do algoritmo para lookup no registry
    parameters: Optional[Dict[str, Any]] = None
    algorithm_class: Optional[type['CSPAlgorithm']] = None

    def __post_init__(self):
        if not self.name and not self.algorithm_class:
            raise ValueError("Either 'name' or 'algorithm_class' must be provided")

    def create_algorithm(self, strings: List[str], alphabet: str) -> 'CSPAlgorithm':
        """
        Create an instance of the CSPAlgorithm with the configured parameters.
        
        Se algorithm_class não estiver definido, usa o AlgorithmRegistry para lookup por nome.
        
        Args:
            strings: List of dataset strings
            alphabet: Alphabet used in the dataset
            
        Returns:
            CSPAlgorithm: Configured algorithm instance
            
        Raises:
            ValueError: If algorithm name not found in registry or algorithm_class not defined
        """
        if self.algorithm_class is None:
            # Import aqui para evitar import circular
            # TODO: Implementar AlgorithmRegistry quando disponível
            raise ValueError(f"AlgorithmRegistry not yet implemented. Please provide algorithm_class directly.")
        
        params = self.parameters or {}
        return self.algorithm_class(strings=strings, alphabet=alphabet, **params)


@dataclass
class AlgorithmConfig:
    """
    Algorithm group configuration matching TEMPLATE.yaml Section 3.
    
    Defines a named group of algorithms with their configurations.
    This allows organizing related algorithm configurations together.
    """
    id: str
    name: str
    description: Optional[str] = None
    algorithms: List[AlgorithmParamsConfig] = field(default_factory=list)

    def __post_init__(self):
        if not self.algorithms or len(self.algorithms) < 1:
            raise ValueError("algorithms must contain at least one algorithm")


@dataclass
class TaskConfig:
    """
    Individual task configuration for execution, optimization, or sensitivity.

    Represents a specific task with its datasets, algorithms, and repetitions.
    Used within ExecutionConfig, OptimizationConfig, and SensitivityConfig.
    """
    id: str
    name: str
    repetitions: int
    datasets: List[DatasetConfig] = field(default_factory=list)
    algorithm_configs: List[AlgorithmConfig] = field(default_factory=list)

    def __post_init__(self):
        # Validar repetições
        if not isinstance(self.repetitions, int) or self.repetitions <= 0:
            raise ValueError("repetitions must be positive integer")
        
        if not self.datasets or len(self.datasets) < 1:
            raise ValueError("datasets must contain at least one DatasetConfig")
        if not self.algorithm_configs or len(self.algorithm_configs) < 1:
            raise ValueError("algorithm_configs must contain at least one AlgorithmConfig")


@dataclass
class ExecutionConfig:
    """
    Execution configuration matching TEMPLATE.yaml Section 5A.
    
    Contains a list of execution tasks to be performed.
    """
    type: BatchType
    tasks: List[TaskConfig] = field(default_factory=list)

    def __post_init__(self):
        if not self.tasks or len(self.tasks) < 1:
            raise ValueError("tasks must contain at least one TaskConfig")


@dataclass
class OptimizationTaskConfig(TaskConfig):
    """
    Optimization task configuration extending base TaskConfig.

    Adds optimization-specific parameters.
    """
    direction: Optional[str] = None
    n_trials: Optional[int] = None
    trial_timeout: Optional[int] = None
    search_space: Dict[str, Any] = field(default_factory=dict)
    extra_parameters: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        # Chamar validação da classe pai primeiro
        super().__post_init__()
        
        # Validações obrigatórias
        if self.direction is None:
            raise ValueError("direction is required")
        if self.n_trials is None:
            raise ValueError("n_trials is required")
        if self.trial_timeout is None:
            raise ValueError("trial_timeout is required")
        
        # Validações específicas de otimização
        if self.direction not in ["minimize", "maximize"]:
            raise ValueError("direction must be 'minimize' or 'maximize'")
        
        if not isinstance(self.n_trials, int) or self.n_trials <= 0:
            raise ValueError("n_trials must be positive integer")
        
        if not isinstance(self.trial_timeout, int) or self.trial_timeout <= 0:
            raise ValueError("trial_timeout must be positive integer")


@dataclass
class OptimizationConfig:
    """
    Optimization configuration matching TEMPLATE.yaml Section 5B.
    
    Configuration for hyperparameter optimization using Optuna.
    """
    method: OptimizationMethod = OptimizationMethod.OPTUNA
    defaults_parameters: Dict[str, Any] = field(default_factory=dict)
    tasks: List[OptimizationTaskConfig] = field(default_factory=list)

    def __post_init__(self):
        if not self.tasks or len(self.tasks) < 1:
            raise ValueError("tasks must contain at least one OptimizationTaskConfig")


@dataclass
class SensitivityTaskConfig(TaskConfig):
    """
    Sensitivity task configuration extending base TaskConfig.
    
    Adds sensitivity analysis-specific parameters.
    """
    method: Optional[SensitivityMethod] = None
    n_samples: Optional[int] = None
    analysis_parameters: Dict[str, Any] = field(default_factory=dict)
    output_metrics: List[str] = field(default_factory=list)
    extra_parameters: Optional[Dict[str, Any]] = None
    
    def __post_init__(self):
        super().__post_init__()
        
        # Validações específicas de sensitivity analysis
        if self.method is None:
            raise ValueError("method is required for sensitivity analysis")
        
        if self.n_samples is not None and (not isinstance(self.n_samples, int) or self.n_samples <= 0):
            raise ValueError("n_samples must be positive integer")
        
        # Validar métricas de saída se especificadas
        if self.output_metrics:
            for metric in self.output_metrics:
                if not isinstance(metric, str) or not metric.strip():
                    raise ValueError("All output_metrics must be non-empty strings")


@dataclass
class SensitivityConfig:
    """
    Sensitivity configuration matching TEMPLATE.yaml Section 5C.
    
    Configuration for sensitivity analysis using SALib.
    """
    method: SensitivityMethod = SensitivityMethod.SALIB  # Using enum for type safety
    defaults_parameters: Dict[str, Any] = field(default_factory=dict)
    tasks: List[SensitivityTaskConfig] = field(default_factory=list)

    def __post_init__(self):
        if not self.tasks or len(self.tasks) < 1:
            raise ValueError("tasks must contain at least one SensitivityTaskConfig")


@dataclass
class ResourceConfig:
    """
    Resource configuration for batch execution.
    
    Defines computational resource limits and constraints.
    """
    timeout_per_repetition: int
    timeout_total_batch: int

    max_memory_gb: Optional[float] = None
    max_cores: Optional[int] = None
    max_workers: Optional[int] = None
    max_internal_threads: Optional[int] = None

    def __post_init__(self):
        # Validar timeouts
        if not isinstance(self.timeout_per_repetition, int) or self.timeout_per_repetition <= 0:
            raise ValueError("timeout_per_repetition must be positive integer")
        
        if not isinstance(self.timeout_total_batch, int) or self.timeout_total_batch <= 0:
            raise ValueError("timeout_total_batch must be positive integer")
        
        # Validar recursos opcionais
        if self.max_memory_gb is not None:
            if not isinstance(self.max_memory_gb, (int, float)) or self.max_memory_gb <= 0:
                raise ValueError("max_memory_gb must be positive number")
        
        if self.max_cores is not None:
            if not isinstance(self.max_cores, int) or self.max_cores <= 0:
                raise ValueError("max_cores must be positive integer")
        
        if self.max_workers is not None:
            if not isinstance(self.max_workers, int) or self.max_workers <= 0:
                raise ValueError("max_workers must be positive integer")
        
        if self.max_internal_threads is not None:
            if not isinstance(self.max_internal_threads, int) or self.max_internal_threads <= 0:
                raise ValueError("max_internal_threads must be positive integer")


@dataclass
class ExportConfig:
    """
    Export configuration for results output.

    Defines how and where results should be exported.
    All fields are required.
    """
    enabled: bool
    formats: List[ExportFormat]
    include_metadata: bool

    def __post_init__(self):
        if not isinstance(self.enabled, bool):
            raise ValueError("enabled must be boolean")
        
        if not isinstance(self.include_metadata, bool):
            raise ValueError("include_metadata must be boolean")
        
        if not self.formats:
            raise ValueError("formats list cannot be empty")
        
        # Validar que todos os formatos são válidos
        for fmt in self.formats:
            if not isinstance(fmt, ExportFormat):
                raise ValueError(f"All formats must be ExportFormat enum values, got: {type(fmt)}")


@dataclass
class LoggingConfig:
    """
    Logging configuration for the batch execution.

    Defines logging levels and output settings.
    Supports section-level logging configuration.
    """
    enabled: bool
    level: str = "INFO"

    def __post_init__(self):
        if not isinstance(self.enabled, bool):
            raise ValueError("enabled must be boolean")
        
        valid_levels = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
        if self.level not in valid_levels:
            raise ValueError(f"level must be one of {valid_levels}, got: {self.level}")

@dataclass
class ReproducibilityConfig:
    """
    Reproducibility configuration for runtime behavior.

    Only includes the random_seed field to guarantee reproducibility.
    """
    random_seed: Optional[int] = None

@dataclass
class BatchConfig:
    """
    Complete batch configuration combining all sections.

    Root configuration object that contains all batch settings.
    """
    metadata: BatchMetadata
    datasets: List[DatasetConfig] = field(default_factory=list)
    algorithm_configs: List[AlgorithmConfig] = field(default_factory=list)
    
    # Execution type - exactly one must be defined
    execution: Optional[ExecutionConfig] = None
    optimization: Optional[OptimizationConfig] = None
    sensitivity: Optional[SensitivityConfig] = None
    
    resource_config: Optional[ResourceConfig] = None
    export_config: Optional[ExportConfig] = None
    logging_config: Optional[LoggingConfig] = None
    reproducibility: Optional[ReproducibilityConfig] = None

    def __post_init__(self):
        # Validar que pelo menos um tipo de execução está definido
        execution_types = [self.execution, self.optimization, self.sensitivity]
        defined_types = [t for t in execution_types if t is not None]
        
        if len(defined_types) == 0:
            raise ValueError("At least one execution type (execution, optimization, or sensitivity) must be defined")
        
        if len(defined_types) > 1:
            raise ValueError("Only one execution type can be defined at a time")
        
        # Validar datasets globais se definidos
        if self.datasets:
            for dataset in self.datasets:
                if not isinstance(dataset, DatasetConfig):
                    raise ValueError("All items in datasets must be DatasetConfig instances")
        
        # Validar algorithm_configs globais se definidos
        if self.algorithm_configs:
            for config in self.algorithm_configs:
                if not isinstance(config, AlgorithmConfig):
                    raise ValueError("All items in algorithm_configs must be AlgorithmConfig instances")

    @property
    def batch_type(self) -> BatchType:
        """Get the batch type based on which execution config is defined."""
        if self.execution is not None:
            return BatchType.EXECUTION
        elif self.optimization is not None:
            return BatchType.OPTIMIZATION
        elif self.sensitivity is not None:
            return BatchType.SENSITIVITY
        else:
            raise ValueError("No execution type defined")
