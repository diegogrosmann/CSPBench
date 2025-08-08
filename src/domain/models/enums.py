"""
CSPBench Domain Enums

Centralized enumerations for CSPBench domain models.
Following Clean Architecture principles with only Python standard library.
"""

from enum import Enum


class BatchType(Enum):
    """Types of batch execution supported by CSPBench."""
    EXECUTION = "execution"
    OPTIMIZATION = "optimization" 
    SENSITIVITY = "sensitivity"


class DatasetType(Enum):
    """Types of datasets supported by CSPBench."""
    SYNTHETIC = "synthetic"
    FILE = "file"
    ENTREZ = "entrez"


class Status(Enum):
    """Generic status for various operations (algorithms, tasks, batches, etc.)."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    TIMEOUT = "timeout"


class ExportFormat(Enum):
    """Supported export formats for results."""
    CSV = "csv"
    JSON = "json"
    PARQUET = "parquet"
    PICKLE = "pickle"


class OptimizationMethod(Enum):
    """Supported optimization methods for hyperparameter tuning."""
    OPTUNA = "optuna"

class SensitivityMethod(Enum):
    """Supported sensitivity analysis methods."""
    SALIB = "SALib"
    MORRIS = "morris"
    SOBOL = "sobol"
    FAST = "fast"
    DELTA = "delta"
    RBD_FAST = "rbd-fast"
    DGSM = "dgsm"
