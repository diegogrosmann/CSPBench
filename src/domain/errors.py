"""
CSPBench Domain Exceptions

Defines custom exceptions for domain and application errors.
"""


class CSPBenchError(Exception):
    """Base exception for all CSPBench exceptions."""

    pass


class DomainError(CSPBenchError):
    """Base exception for domain errors."""

    pass


class ApplicationError(CSPBenchError):
    """Base exception for application errors."""

    pass


# Dataset Exceptions
class DatasetError(DomainError):
    """Error related to datasets."""

    pass


class DatasetNotFoundError(DatasetError):
    """Dataset not found."""

    pass


class DatasetValidationError(DatasetError):
    """Dataset validation error."""

    pass


class DatasetEmptyError(DatasetError):
    """Empty dataset."""

    pass


# Algorithm Exceptions
class AlgorithmError(DomainError):
    """Error related to algorithms."""

    pass


class AlgorithmNotFoundError(AlgorithmError):
    """Algorithm not found."""

    pass


class AlgorithmExecutionError(AlgorithmError):
    """Error in algorithm execution."""

    pass


class AlgorithmRegistrationError(AlgorithmError):
    """Error in algorithm registration."""

    pass


class AlgorithmTimeoutError(AlgorithmError):
    """Timeout in algorithm execution."""

    pass


class AlgorithmParameterError(AlgorithmError):
    """Error in algorithm parameters."""

    pass


# Configuration Exceptions
class ConfigurationError(ApplicationError):
    """Configuration error."""

    pass


class BatchConfigurationError(ConfigurationError):
    """Error in batch configuration."""

    pass


class OptimizationConfigurationError(ConfigurationError):
    """Error in optimization configuration."""

    pass


class SensitivityConfigurationError(ConfigurationError):
    """Error in sensitivity analysis configuration."""

    pass


# Execution Exceptions
class ExecutionError(ApplicationError):
    """Execution error."""

    pass


class BatchExecutionError(ExecutionError):
    """Error in batch execution."""

    pass


class OptimizationExecutionError(ExecutionError):
    """Error in optimization execution."""

    pass


class SensitivityExecutionError(ExecutionError):
    """Error in sensitivity analysis execution."""

    pass


# Export Exceptions
class ExportError(ApplicationError):
    """Export error."""

    pass


class UnsupportedFormatError(ExportError):
    """Unsupported export format."""

    pass


class ExportDestinationError(ExportError):
    """Error in export destination."""

    pass


# Repository Exceptions
class RepositoryError(ApplicationError):
    """Repository error."""

    pass


class RepositoryConnectionError(RepositoryError):
    """Repository connection error."""

    pass


class RepositoryPermissionError(RepositoryError):
    """Repository permission error."""

    pass
