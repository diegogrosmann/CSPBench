"""
CSPBench Domain

Domain module containing the core business logic for CSPBench.
Implements algorithms, metrics, and data entities without external dependencies.
"""

from .algorithms import CSPAlgorithm, global_registry, register_algorithm
from .dataset import Dataset
from .distance import (
    DistanceCalculator,
    HammingDistanceCalculator,
    LevenshteinDistanceCalculator,
    create_distance_calculator,
)
from .work import WorkStatus, WorkItem

from .errors import (
    AlgorithmError,
    AlgorithmExecutionError,
    AlgorithmNotFoundError,
    AlgorithmParameterError,
    AlgorithmRegistrationError,
    AlgorithmTimeoutError,
    ApplicationError,
    BatchConfigurationError,
    BatchExecutionError,
    ConfigurationError,
    CSPBenchError,
    DatasetEmptyError,
    DatasetError,
    DatasetNotFoundError,
    DatasetValidationError,
    DomainError,
    ExecutionError,
    OptimizationConfigurationError,
    OptimizationExecutionError,
    SensitivityConfigurationError,
    SensitivityExecutionError,
)

__all__ = [
    # Algorithms
    "CSPAlgorithm",
    "register_algorithm",
    "global_registry",
    # Distance calculators
    "DistanceCalculator",
    "HammingDistanceCalculator",
    "LevenshteinDistanceCalculator",
    "create_distance_calculator",
    # Work entities
    "WorkStatus",
    "WorkItem",
    # Metrics and quality evaluation
    "diversity_metric",
    "consensus_strength",
    "solution_quality",
    "QualityEvaluator",
    # Datasets
    "Dataset",
    "SyntheticDatasetGenerator",
    # Errors
    "CSPBenchError",
    "DomainError",
    "ApplicationError",
    "DatasetError",
    "DatasetNotFoundError",
    "DatasetValidationError",
    "DatasetEmptyError",
    "AlgorithmError",
    "AlgorithmNotFoundError",
    "AlgorithmExecutionError",
    "AlgorithmRegistrationError",
    "AlgorithmTimeoutError",
    "AlgorithmParameterError",
    "ConfigurationError",
    "BatchConfigurationError",
    "OptimizationConfigurationError",
    "SensitivityConfigurationError",
    "ExecutionError",
    "BatchExecutionError",
    "OptimizationExecutionError",
    "SensitivityExecutionError",
]


# Lazy re-export to avoid circular import with application layer
def __getattr__(name):  # PEP 562
    if name == "SyntheticDatasetGenerator":
        from src.application.services.dataset_generator import (
            SyntheticDatasetGenerator as _SDG,
        )

        return _SDG
    raise AttributeError(name)
