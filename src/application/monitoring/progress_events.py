"""Progress events for the monitoring system."""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Any, Dict, Optional


class TaskType(Enum):
    """Supported task types."""

    EXECUTION = "execution"
    OPTIMIZATION = "optimization"
    SENSITIVITY = "sensitivity"


class TaskStatus(Enum):
    """Task status."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ItemStatus(Enum):
    """Individual item status (repetition, trial, sample)."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ExecutionLevel(Enum):
    """Hierarchical execution levels."""

    EXECUTION = "execution"
    DATASET = "dataset"
    CONFIG = "config"
    ALGORITHM = "algorithm"
    REPETITION = "repetition"
    TRIAL = "trial"
    SAMPLE = "sample"


class ExecutionPhase(Enum):
    """Execution phases."""

    STARTING = "starting"
    RUNNING = "running"
    FINISHING = "finishing"
    COMPLETED = "completed"
    FAILED = "failed"


class UnifiedPhase(Enum):
    """Unified high-level phases for display across the system."""

    PROCESSING = "processing"
    OPTIMIZATION = "optimization"
    ANALYSIS = "analysis"


@dataclass
class ProgressEvent:
    """Base class for all progress events."""

    session_id: Optional[str] = None
    timestamp: Optional[datetime] = field(default_factory=datetime.now)


@dataclass
class TaskStartedEvent(ProgressEvent):
    """Event fired when a task starts."""

    task_type: TaskType = TaskType.EXECUTION
    task_name: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TaskFinishedEvent(ProgressEvent):
    """Event fired when a task finishes."""

    task_type: TaskType = TaskType.EXECUTION
    task_name: str = ""
    success: bool = True
    results: Dict[str, Any] = field(default_factory=dict)
    error_message: Optional[str] = None


@dataclass
class ExecutionStartedEvent(ProgressEvent):
    """Event fired when an execution block starts."""

    execution_name: str = ""
    total_items: int = 0
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ExecutionProgressEvent(ProgressEvent):
    """Event fired during execution progress."""

    execution_name: str = ""
    current_item: int = 0
    total_items: int = 0
    item_name: str = ""
    progress_percent: float = 0.0
    message: str = ""
    context: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ExecutionFinishedEvent(ProgressEvent):
    """Event fired when an execution block finishes."""

    execution_name: str = ""
    success: bool = True
    total_processed: int = 0
    results: Dict[str, Any] = field(default_factory=dict)
    error_message: Optional[str] = None


@dataclass
class AlgorithmProgressEvent(ProgressEvent):
    """Event fired during algorithm execution."""

    algorithm_name: str = ""
    progress_percent: float = 0.0
    message: str = ""
    item_id: Optional[str] = None
    context: Dict[str, Any] = field(default_factory=dict)


@dataclass
class AlgorithmFinishedEvent(ProgressEvent):
    """Event fired when an algorithm execution finishes."""

    algorithm_name: str = ""
    success: bool = True
    best_string: str = ""
    max_distance: int = -1
    execution_time: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)
    repetition_number: Optional[int] = None


@dataclass
class WarningEvent(ProgressEvent):
    """Event fired when an algorithm reports a warning."""

    algorithm_name: str = ""
    warning_message: str = ""
    item_id: Optional[str] = None
    repetition_number: Optional[int] = None
    execution_context: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ErrorEvent(ProgressEvent):
    """Event fired when an error occurs."""

    error_message: str = ""
    error_type: str = "generic"
    context: Dict[str, Any] = field(default_factory=dict)


@dataclass
class DisplayEvent(ProgressEvent):
    """Unified display event consumed by a single display for all phases.

    This keeps a minimal, phase-agnostic contract and carries extra data in payload.
    """

    phase: UnifiedPhase = UnifiedPhase.PROCESSING
    message: str = ""
    progress: float = 0.0  # 0..1 progress for the current unit
    dataset_id: Optional[str] = None
    algorithm_name: Optional[str] = None
    task_id: Optional[str] = None
    trial_no: Optional[int] = None
    rep_idx: Optional[int] = None
    payload: Dict[str, Any] = field(default_factory=dict)
