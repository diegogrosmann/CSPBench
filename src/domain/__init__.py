"""
CSPBench Domain

Domain module containing the core business logic for CSPBench.
Implements algorithms, metrics, and data entities without external dependencies.
"""

from .algorithms import Algorithm, CSPAlgorithm, global_registry, register_algorithm
from .dataset import Dataset
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
from .metrics import (
    DistanceCalculator,
    QualityEvaluator,
    average_distance,
    consensus_strength,
    diversity_metric,
    hamming_distance,
    max_distance,
    median_distance,
    solution_quality,
)

__all__ = [
    # Algorithms
    "CSPAlgorithm",
    "Algorithm",
    "register_algorithm",
    "global_registry",
    # Metrics
    "hamming_distance",
    "max_distance",
    "average_distance",
    "median_distance",
    "diversity_metric",
    "consensus_strength",
    "solution_quality",
    "DistanceCalculator",
    "QualityEvaluator",
    # Dataset
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
