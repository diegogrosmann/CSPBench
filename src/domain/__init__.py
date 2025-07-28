"""
CSPBench Domain

Domain module containing the core business logic for CSPBench.
Implements algorithms, metrics, and data entities without external dependencies.
"""

from .algorithms import Algorithm, CSPAlgorithm, global_registry, register_algorithm
from .dataset import Dataset, SyntheticDatasetGenerator
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
    max_hamming,
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
    "max_hamming",
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
