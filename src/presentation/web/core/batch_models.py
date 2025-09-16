"""
Batch models for CSPBench Web Interface.

This module provides Pydantic models for batch configuration management including
batch file information, upload requests, and operation responses. It handles
validation and serialization of batch-related data structures.
"""

from typing import List, Optional

import yaml
from pydantic import BaseModel, Field, field_validator


class BatchFileInfo(BaseModel):
    """Batch file information model.
    
    This model represents metadata and information about batch configuration
    files stored in the system, including file properties and embedded metadata.
    
    Attributes:
        name (str): Batch file name.
        description (Optional[str]): File description.
        path (str): Filesystem path to the batch file.
        size (str): Human-readable file size.
        created (str): Creation timestamp string.
        modified (str): Last modification timestamp string.
        is_template (bool): Whether this is a template file.
        metadata_name (Optional[str]): Name from batch metadata.
        metadata_description (Optional[str]): Description from batch metadata.
        metadata_author (Optional[str]): Author from batch metadata.
        metadata_version (Optional[str]): Version from batch metadata.
        metadata_creation_date (Optional[str]): Creation date from metadata.
        metadata_tags (List[str]): Tags from batch metadata.
    """

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
    """Batch files list response model.
    
    This model provides a response structure for listing batch files,
    including pagination information and file metadata.
    
    Attributes:
        files (List[BatchFileInfo]): List of batch file information objects.
        total (int): Total number of batch files available.
    """

    files: List[BatchFileInfo]
    total: int


class BatchUploadRequest(BaseModel):
    """Request for uploading batch configuration content.
    
    This model handles and validates requests to upload new batch
    configuration files to the system with proper YAML validation.
    
    Attributes:
        name (str): Name for the batch configuration file.
        content (str): YAML batch configuration content.
        overwrite (bool): Whether to overwrite existing files with same name.
    """

    name: str = Field(..., description="Batch configuration name")
    content: str = Field(..., description="YAML batch configuration content")
    overwrite: bool = Field(False, description="Whether to overwrite existing file")

    @field_validator("content")
    @classmethod
    def validate_yaml_content(cls, v: str):
        """Validate YAML content format.
        
        Ensures the uploaded content is valid YAML format that can
        be parsed successfully.
        
        Args:
            v (str): YAML content to validate.
            
        Returns:
            str: Validated YAML content.
            
        Raises:
            ValueError: If content is empty or invalid YAML format.
        """
        if not v.strip():
            raise ValueError("Content cannot be empty")

        try:
            # Try to parse as YAML to validate format
            yaml.safe_load(v)
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML format: {str(e)}")

        return v

    @field_validator("name")
    @classmethod
    def validate_name(cls, v: str):
        """Validate batch configuration name.
        
        Ensures the batch name is safe for filesystem operations
        by removing potentially dangerous characters.
        
        Args:
            v (str): Batch name to validate.
            
        Returns:
            str: Validated and sanitized batch name.
            
        Raises:
            ValueError: If name is empty or contains invalid characters.
        """
        if not v.strip():
            raise ValueError("Name cannot be empty")

        # Remove potentially dangerous characters
        invalid_chars = ["/", "\\", "..", "<", ">", ":", '"', "|", "?", "*"]
        for char in invalid_chars:
            if char in v:
                raise ValueError(f"Name cannot contain '{char}'")

        return v.strip()


class OperationResponse(BaseModel):
    """Generic operation response model.
    
    This model provides a standardized response format for various
    batch operations including success status, messages, and optional data.
    
    Attributes:
        success (bool): Whether the operation completed successfully.
        message (str): Human-readable operation result message.
        data (Optional[dict]): Additional response data if applicable.
    """

    success: bool
    message: str
    data: Optional[dict] = None
