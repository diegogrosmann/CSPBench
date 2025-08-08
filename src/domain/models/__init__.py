"""
Domain Models Package

Exports all domain models for clean imports.
"""

# Import batch models
from .batch import (
    BatchConfig,
    BatchMetadata,
    DatasetConfig,
    AlgorithmConfig,
    AlgorithmParamsConfig,
    TaskConfig,
    ExecutionConfig,
    OptimizationConfig,
    OptimizationTaskConfig,
    SensitivityConfig,
    SensitivityTaskConfig,
    ResourceConfig,
    ExportConfig,
    LoggingConfig,
    ReproducibilityConfig,
)

# Import result models
from .results import (
    AlgorithmResult,
    OptimizationResult,
    SensitivityResult,
    BatchResults,
)

# Import event models
from .events import (
    BaseEvent,
    ErrorEvent,
    BatchStartedEvent,
    BatchFinishedEvent,
    RepetitionEvent,
    RepetitionStartedEvent,
    RepetitionProgressEvent,
    RepetitionCallbackEvent,
    RepetitionFinishedEvent,
)

# Import enums
from .enums import (
    BatchType,
    DatasetType,
    Status,
    ExportFormat,
    OptimizationMethod,
    SensitivityMethod,
)

__all__ = [
    # Batch models
    "BatchConfig",
    "BatchMetadata",
    "DatasetConfig",
    "AlgorithmConfig",
    "AlgorithmParamsConfig",
    "TaskConfig",
    "ExecutionConfig",
    "OptimizationConfig",
    "OptimizationTaskConfig",
    "SensitivityConfig",
    "SensitivityTaskConfig",
    "ResourceConfig",
    "ExportConfig",
    "LoggingConfig",
    "ReproducibilityConfig",
    # Results models
    "AlgorithmResult",
    "OptimizationResult",
    "SensitivityResult",
    "BatchResults",
    # Event models
    "BaseEvent",
    "ErrorEvent",
    "BatchStartedEvent",
    "BatchFinishedEvent",
    "RepetitionEvent",
    "RepetitionStartedEvent",
    "RepetitionProgressEvent",
    "RepetitionCallbackEvent",
    "RepetitionFinishedEvent",
    # Enums
    "BatchType",
    "DatasetType",
    "Status",
    "ExportFormat",
    "OptimizationMethod",
    "SensitivityMethod",
]
