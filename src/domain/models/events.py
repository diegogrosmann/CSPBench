"""
CSPBench Event Models

Domain models for system events during batch execution.
These events are used for monitoring and progress tracking.
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, TYPE_CHECKING

from .enums import BatchType, Status

if TYPE_CHECKING:
    from .batch import BatchConfig


@dataclass
class BaseEvent:
    """
    Base class for all CSPBench events.
    
    All events in the system inherit from this base class.
    Includes optional status tracking for event lifecycle management.
    """
    session_id: Optional[str] = None
    timestamp: datetime = field(default_factory=datetime.now)

@dataclass
class ErrorEvent:
    """
    Event fired when an error occurs during execution.
    
    Provides detailed error information for debugging and monitoring.
    """
    error_message: str = ""
    error_type: str = ""
    context: Optional[Dict[str, Any]] = None
    stack_trace: Optional[str] = None


@dataclass
class BatchStartedEvent(BaseEvent):
    """
    Event fired when a batch starts execution.
    
    Contains comprehensive information about the batch being executed.
    """
    batch_config: Optional['BatchConfig'] = None
    status: Status = Status.RUNNING

    def __post_init__(self):
        if self.batch_config is None:
            raise ValueError("batch_config cannot be None")

@dataclass
class BatchFinishedEvent(BaseEvent):
    """
    Event fired when a batch finishes execution.
    
    Contains execution results and statistics.
    """
    success: Optional[bool] = None
    status: Optional[Status] = None
    total_duration: Optional[float] = None
    results_summary: Optional[Dict[str, Any]] = None
    errors: Optional[List[ErrorEvent]] = None
    
    def __post_init__(self):
        if self.success is None:
            raise ValueError("success must be provided")
        if self.status is None:
            raise ValueError("status must be provided")
        if self.total_duration is None:
            raise ValueError("total_duration must be provided")

    @staticmethod
    def completed(
        session_id: str,
        total_duration: float,
        results_summary: Optional[Dict[str, Any]] = None
    ) -> 'BatchFinishedEvent':
        """Create a BatchFinishedEvent for successful completion."""
        return BatchFinishedEvent(
            session_id=session_id,
            success=True,
            total_duration=total_duration,
            results_summary=results_summary or {},
            status=Status.COMPLETED
        )

    @staticmethod
    def failed(
        session_id: str,
        total_duration: float,
        errors: List[ErrorEvent],
        results_summary: Optional[Dict[str, Any]] = None
    ) -> 'BatchFinishedEvent':
        """Create a BatchFinishedEvent for failed execution."""
        return BatchFinishedEvent(
            session_id=session_id,
            success=False,
            total_duration=total_duration,
            errors=errors,
            results_summary=results_summary or {},
            status=Status.FAILED
        )

    @staticmethod
    def timeout(
        session_id: str,
        total_duration: float,
        results_summary: Optional[Dict[str, Any]] = None
    ) -> 'BatchFinishedEvent':
        """Create a BatchFinishedEvent for timeout scenario."""
        return BatchFinishedEvent(
            session_id=session_id,
            success=False,
            total_duration=total_duration,
            results_summary=results_summary or {},
            status=Status.TIMEOUT
        )

@dataclass
class RepetitionEvent(BaseEvent):
    """
    Base event for repetition-related events.
    """
    repetition_id: Optional[str] = None

    def __post_init__(self):
        if self.repetition_id is None:
            raise ValueError("repetition_id cannot be None")


@dataclass
class RepetitionStartedEvent(RepetitionEvent):
    """
    Event fired when a repetition starts execution.
    
    Contains initialization parameters and configuration details.
    """
    datasets_id: Optional[str] = None
    algorithm_config_id: Optional[str] = None
    algorithm_name: Optional[str] = None
    repetition_number: Optional[int] = None
    status: Status = Status.RUNNING

    def __post_init__(self):
        if self.datasets_id is None:
            raise ValueError("datasets_id cannot be None")
        if self.algorithm_config_id is None:
            raise ValueError("algorithm_config_id cannot be None")
        if self.algorithm_name is None:
            raise ValueError("algorithm_name cannot be None")
        if self.repetition_number is None:
            raise ValueError("repetition_number cannot be None")
    
    
@dataclass
class RepetitionProgressEvent(RepetitionEvent):
    """
    Event fired to report progress during repetition execution.
    
    Used for monitoring and progress tracking of long-running operations.
    """
    progress_percentage: Optional[float] = None
    extra_info: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        super().__post_init__()
        if self.progress_percentage is None:
            raise ValueError("progress_percentage must be provided")


@dataclass
class RepetitionCallbackEvent(RepetitionEvent):
    """
    Event fired for custom callback notifications during repetition execution.
    
    Allows algorithms to send custom events and data during processing.
    """
    callback_type: Optional[str] = None
    callback_data: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        super().__post_init__()
        if not self.callback_type:
            raise ValueError("callback_type cannot be empty")
        if not self.callback_data:
            raise ValueError("callback_data must be provided")


@dataclass
class RepetitionFinishedEvent(RepetitionEvent):
    """
    Event fired when a repetition finishes execution.
    
    Contains execution results, metrics, and performance statistics.
    """
    success: Optional[bool] = None
    status: Optional[Status] = None
    duration: Optional[float] = None
    results: Optional[Dict[str, Any]] = None
    metrics: Optional[Dict[str, float]] = None
    errors: Optional[List[ErrorEvent]] = None
    
    def __post_init__(self):
        super().__post_init__()
        if self.success is None:
            raise ValueError("success must be provided")
        if self.status is None:
            raise ValueError("status must be provided")
        if self.duration is None:
            raise ValueError("duration must be provided")

    @staticmethod
    def completed(
        session_id: str,
        repetition_id: str,
        duration: float,
        results: Optional[Dict[str, Any]] = None,
        metrics: Optional[Dict[str, float]] = None
    ) -> 'RepetitionFinishedEvent':
        """Create a RepetitionFinishedEvent for successful completion."""
        return RepetitionFinishedEvent(
            session_id=session_id,
            repetition_id=repetition_id,
            success=True,
            status=Status.COMPLETED,
            duration=duration,
            results=results or {},
            metrics=metrics or {}
        )

    @staticmethod
    def failed(
        session_id: str,
        repetition_id: str,
        duration: float,
        errors: List[ErrorEvent],
        results: Optional[Dict[str, Any]] = None,
        metrics: Optional[Dict[str, float]] = None,
        
    ) -> 'RepetitionFinishedEvent':
        """Create a RepetitionFinishedEvent for failed execution."""
        return RepetitionFinishedEvent(
            session_id=session_id,
            repetition_id=repetition_id,
            success=False,
            status=Status.FAILED,
            duration=duration,
            results=results or {},
            metrics=metrics or {},
            errors=errors
        )

    @staticmethod
    def timeout(
        session_id: str,
        repetition_id: str,
        duration: float,
        results: Optional[Dict[str, Any]] = None,
        metrics: Optional[Dict[str, float]] = None
    ) -> 'RepetitionFinishedEvent':
        """Create a RepetitionFinishedEvent for timeout scenario."""
        return RepetitionFinishedEvent(
            session_id=session_id,
            repetition_id=repetition_id,
            success=False,
            status=Status.TIMEOUT,
            duration=duration,
            results=results or {},
            metrics=metrics or {},
            errors=[]
        )

