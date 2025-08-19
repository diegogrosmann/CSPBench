"""
Batch Execution models for CSPBench Web Interface.

This module provides Pydantic models for batch execution management including:
- Batch execution requests and responses
- Status tracking and monitoring
- Results and progress reporting
"""

from datetime import datetime
from typing import Dict, List, Optional, Any
from pydantic import BaseModel, Field, validator
import yaml


class BatchExecutionRequest(BaseModel):
    """Request model for batch execution."""
    
    batch_file: str = Field(..., description="Path to batch YAML file (relative to batches/ directory)")
    monitor_type: str = Field("log", description="Monitor type: 'none', 'log', 'websocket'")
    
    @validator('batch_file')
    def validate_batch_file(cls, v):
        """Validate batch file path."""
        if not v.strip():
            raise ValueError("Batch file path cannot be empty")
        
        # Basic path validation
        if '..' in v or v.startswith('/'):
            raise ValueError("Invalid batch file path. Use relative paths only.")
        
        if not v.endswith('.yaml') and not v.endswith('.yml'):
            raise ValueError("Batch file must have .yaml or .yml extension")
        
        return v.strip()
    
    @validator('monitor_type')
    def validate_monitor_type(cls, v):
        """Validate monitor type."""
        valid_types = ['none', 'log', 'websocket']
        if v not in valid_types:
            raise ValueError(f"Invalid monitor type. Must be one of: {valid_types}")
        return v


class BatchExecutionResponse(BaseModel):
    """Response model for batch execution submission."""
    
    work_id: str = Field(..., description="Unique work identifier")
    status: str = Field(..., description="Initial execution status")
    message: str = Field(..., description="Submission status message")
    submitted_at: datetime = Field(..., description="Submission timestamp")
    output_path: Optional[str] = Field(None, description="Path where results will be stored")


class BatchStatusResponse(BaseModel):
    """Response model for batch execution status."""
    
    work_id: str = Field(..., description="Unique work identifier")
    status: str = Field(..., description="Current execution status")
    created_at: datetime = Field(..., description="Work creation timestamp")
    updated_at: datetime = Field(..., description="Last update timestamp")
    output_path: Optional[str] = Field(None, description="Output directory path")
    progress: Optional[Dict[str, Any]] = Field(None, description="Execution progress data")
    error: Optional[str] = Field(None, description="Error message if status is error")
    config_name: Optional[str] = Field(None, description="Batch configuration name")


class BatchResultsResponse(BaseModel):
    """Response model for batch execution results."""
    
    work_id: str = Field(..., description="Unique work identifier")
    status: str = Field(..., description="Final execution status")
    output_path: Optional[str] = Field(None, description="Results directory path")
    output_files: List[str] = Field(default_factory=list, description="List of output files")
    summary: Optional[Dict[str, Any]] = Field(None, description="Execution summary")
    config_name: Optional[str] = Field(None, description="Batch configuration name")
    execution_time: Optional[float] = Field(None, description="Total execution time in seconds")


class BatchControlResponse(BaseModel):
    """Response model for batch control operations (cancel, pause, resume)."""
    
    work_id: str = Field(..., description="Unique work identifier")
    operation: str = Field(..., description="Operation performed")
    success: bool = Field(..., description="Whether operation was successful")
    message: str = Field(..., description="Operation result message")
    new_status: Optional[str] = Field(None, description="New status after operation")


class BatchListResponse(BaseModel):
    """Response model for listing batch executions."""
    
    executions: List[BatchStatusResponse] = Field(..., description="List of batch executions")
    total: int = Field(..., description="Total number of executions")
    filtered: Optional[int] = Field(None, description="Number after filtering")


class BatchProgressUpdate(BaseModel):
    """Model for real-time progress updates via WebSocket."""
    
    work_id: str = Field(..., description="Unique work identifier")
    timestamp: datetime = Field(..., description="Update timestamp")
    event_type: str = Field(..., description="Type of progress event")
    data: Dict[str, Any] = Field(..., description="Progress data")
    message: Optional[str] = Field(None, description="Human-readable message")
