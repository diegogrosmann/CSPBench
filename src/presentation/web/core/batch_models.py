"""
Batch models for CSPBench Web Interface.

This module provides Pydantic models for batch configuration management including:
- Batch file information
- Batch upload requests
- Batch operation responses
"""

from typing import List, Optional
from pydantic import BaseModel, Field, validator
import yaml


class BatchFileInfo(BaseModel):
    """Batch file information model."""

    name: str
    description: Optional[str] = None
    path: str
    size: str
    created: str
    modified: str
    is_template: bool = False
    metadata_name: Optional[str] = None
    metadata_description: Optional[str] = None
    metadata_author: Optional[str] = None
    metadata_version: Optional[str] = None
    metadata_creation_date: Optional[str] = None
    metadata_tags: List[str] = []


class BatchListResponse(BaseModel):
    """Batch files list response model."""

    files: List[BatchFileInfo]
    total: int


class BatchUploadRequest(BaseModel):
    """Request for uploading batch configuration content."""

    name: str = Field(..., description="Batch configuration name")
    content: str = Field(..., description="YAML batch configuration content")
    overwrite: bool = Field(False, description="Whether to overwrite existing file")

    @validator("content")
    def validate_yaml_content(cls, v):
        """Validate YAML content format."""
        if not v.strip():
            raise ValueError("Content cannot be empty")

        try:
            # Try to parse as YAML to validate format
            yaml.safe_load(v)
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML format: {str(e)}")

        return v

    @validator("name")
    def validate_name(cls, v):
        """Validate batch configuration name."""
        if not v.strip():
            raise ValueError("Name cannot be empty")

        # Remove potentially dangerous characters
        invalid_chars = ["/", "\\", "..", "<", ">", ":", '"', "|", "?", "*"]
        for char in invalid_chars:
            if char in v:
                raise ValueError(f"Name cannot contain '{char}'")

        return v.strip()


class OperationResponse(BaseModel):
    """Generic operation response model."""

    success: bool
    message: str
    data: Optional[dict] = None
