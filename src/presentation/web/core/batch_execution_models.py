"""
Batch Execution models for CSPBench Web Interface.

This module provides comprehensive Pydantic models for batch execution management
including execution requests, status tracking, progress monitoring, and results
reporting for the CSPBench batch processing system.
"""

from datetime import datetime
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


class BatchExecutionRequest(BaseModel):
    """Request model for batch execution.

    This model handles requests to execute batch configurations with
    proper validation of file paths and monitoring options.

    Attributes:
        batch_file (str): Relative path to batch YAML file in batches/ directory.
        monitor_type (str): Type of monitoring ('none', 'log', 'websocket').
    """

    batch_file: str = Field(
        ..., description="Path to batch YAML file (relative to batches/ directory)"
    )
    monitor_type: str = Field(
        "log", description="Monitor type: 'none', 'log', 'websocket'"
    )

    @field_validator("batch_file")
    @classmethod
    def validate_batch_file(cls, v: str):
        """Validate batch file path.

        Ensures the batch file path is safe and points to a valid
        YAML configuration file within the allowed directory structure.

        Args:
            v (str): Batch file path to validate.

        Returns:
            str: Validated and sanitized file path.

        Raises:
            ValueError: If path is empty, contains unsafe patterns, or wrong extension.
        """
        if not v.strip():
            raise ValueError("Batch file path cannot be empty")

        # Basic path validation
        if ".." in v or v.startswith("/"):
            raise ValueError("Invalid batch file path. Use relative paths only.")

        if not v.endswith(".yaml") and not v.endswith(".yml"):
            raise ValueError("Batch file must have .yaml or .yml extension")

        return v.strip()

    @field_validator("monitor_type")
    @classmethod
    def validate_monitor_type(cls, v: str):
        """Validate monitor type.

        Ensures the monitoring type is one of the supported options
        for batch execution tracking.

        Args:
            v (str): Monitor type to validate.

        Returns:
            str: Validated monitor type.

        Raises:
            ValueError: If monitor type is not in the allowed list.
        """
        valid_types = ["none", "log", "websocket"]
        if v not in valid_types:
            raise ValueError(f"Invalid monitor type. Must be one of: {valid_types}")
        return v


class BatchExecutionResponse(BaseModel):
    """Response model for batch execution submission.

    This model provides information about a successfully submitted
    batch execution including tracking details and output location.

    Attributes:
        work_id (str): Unique identifier for tracking this execution.
        status (str): Initial execution status after submission.
        message (str): Human-readable submission status message.
        submitted_at (datetime): Timestamp when execution was submitted.
        output_path (Optional[str]): Directory where results will be stored.
    """

    work_id: str = Field(..., description="Unique work identifier")
    status: str = Field(..., description="Initial execution status")
    message: str = Field(..., description="Submission status message")
    submitted_at: datetime = Field(..., description="Submission timestamp")
    output_path: Optional[str] = Field(
        None, description="Path where results will be stored"
    )


class BatchStatusResponse(BaseModel):
    """Response model for batch execution status.

    This model provides current status information for batch executions
    including progress, timing, and error information.

    Attributes:
        work_id (str): Unique identifier for this execution.
        status (str): Current execution status.
        created_at (datetime): When the work was initially created.
        updated_at (datetime): When status was last updated.
        output_path (Optional[str]): Directory containing output files.
        progress (Optional[Dict]): Detailed progress information.
        error (Optional[str]): Error message if execution failed.
        config_name (Optional[str]): Name of the batch configuration used.
    """

    work_id: str = Field(..., description="Unique work identifier")
    status: str = Field(..., description="Current execution status")
    created_at: datetime = Field(..., description="Work creation timestamp")
    updated_at: datetime = Field(..., description="Last update timestamp")
    output_path: Optional[str] = Field(None, description="Output directory path")
    progress: Optional[Dict[str, Any]] = Field(
        None, description="Execution progress data"
    )
    error: Optional[str] = Field(None, description="Error message if status is error")
    config_name: Optional[str] = Field(None, description="Batch configuration name")


class BatchResultsResponse(BaseModel):
    """Response model for batch execution results.

    This model provides comprehensive information about completed
    batch executions including output files and execution statistics.

    Attributes:
        work_id (str): Unique identifier for this execution.
        status (str): Final execution status.
        output_path (Optional[str]): Directory containing result files.
        output_files (List[str]): List of generated output files.
        summary (Optional[Dict]): Execution summary and statistics.
        config_name (Optional[str]): Name of the batch configuration used.
        execution_time (Optional[float]): Total execution time in seconds.
    """

    work_id: str = Field(..., description="Unique work identifier")
    status: str = Field(..., description="Final execution status")
    output_path: Optional[str] = Field(None, description="Results directory path")
    output_files: List[str] = Field(
        default_factory=list, description="List of output files"
    )
    summary: Optional[Dict[str, Any]] = Field(None, description="Execution summary")
    config_name: Optional[str] = Field(None, description="Batch configuration name")
    execution_time: Optional[float] = Field(
        None, description="Total execution time in seconds"
    )


class BatchControlResponse(BaseModel):
    """Response model for batch control operations (cancel, pause, restart).

    This model provides feedback for control operations performed on
    batch executions such as cancellation, pausing, or restarting.

    Attributes:
        work_id (str): Unique identifier for the affected execution.
        operation (str): Type of control operation performed.
        success (bool): Whether the control operation was successful.
        message (str): Human-readable result of the operation.
        new_status (Optional[str]): New execution status after operation.
    """

    work_id: str = Field(..., description="Unique work identifier")
    operation: str = Field(..., description="Operation performed")
    success: bool = Field(..., description="Whether operation was successful")
    message: str = Field(..., description="Operation result message")
    new_status: Optional[str] = Field(None, description="New status after operation")


class BatchListResponse(BaseModel):
    """Response model for listing batch executions.

    This model provides a paginated list of batch executions with
    filtering and total count information.

    Attributes:
        executions (List[BatchStatusResponse]): List of batch execution statuses.
        total (int): Total number of executions in the system.
        filtered (Optional[int]): Number of executions after applying filters.
    """

    executions: List[BatchStatusResponse] = Field(
        ..., description="List of batch executions"
    )
    total: int = Field(..., description="Total number of executions")
    filtered: Optional[int] = Field(None, description="Number after filtering")


class BatchProgressUpdate(BaseModel):
    """Model for real-time progress updates via WebSocket.

    This model represents real-time progress updates sent through
    WebSocket connections to provide live execution monitoring.

    Attributes:
        work_id (str): Unique identifier for the execution being updated.
        timestamp (datetime): When this progress update was generated.
        event_type (str): Type of progress event (started, progress, completed, etc.).
        data (Dict[str, Any]): Structured progress data specific to event type.
        message (Optional[str]): Human-readable progress message.
    """

    work_id: str = Field(..., description="Unique work identifier")
    timestamp: datetime = Field(..., description="Update timestamp")
    event_type: str = Field(..., description="Type of progress event")
    data: Dict[str, Any] = Field(..., description="Progress data")
    message: Optional[str] = Field(None, description="Human-readable message")
