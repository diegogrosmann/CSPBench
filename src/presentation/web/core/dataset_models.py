"""
Dataset-related API models for CSPBench Web Interface.

This module provides comprehensive Pydantic models for dataset management,
including creation, validation, metadata handling, and various dataset
generation methods (synthetic, file upload, NCBI download).
"""

from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


class DatasetType(str, Enum):
    """Types of datasets supported by the system.
    
    This enum defines the different methods available for creating
    datasets in the CSPBench system.
    """

    SYNTHETIC = "synthetic"  # Algorithmically generated datasets
    FILE = "file"           # User-uploaded FASTA files
    NCBI = "ncbi"          # Downloaded from NCBI databases


class SyntheticMethod(str, Enum):
    """Synthetic dataset generation methods.
    
    This enum defines the available algorithms for generating
    synthetic biological sequence datasets.
    """

    RANDOM = "random"        # Completely random sequences
    NOISE = "noise"         # Add noise to base sequence
    CLUSTERED = "clustered" # Generate clustered sequences
    MUTATIONS = "mutations" # Apply mutations to base sequence


class DatasetInfo(BaseModel):
    """Dataset information response model.
    
    This model provides comprehensive metadata and statistics
    about a dataset, including sequence statistics, creation
    information, and storage details.
    
    Attributes:
        id (str): Unique dataset identifier.
        name (str): Human-readable dataset name.
        type (DatasetType): Method used to create the dataset.
        size (int): Total number of sequences in dataset.
        min_length (int): Length of shortest sequence.
        max_length (int): Length of longest sequence.
        average_length (float): Mean sequence length.
        alphabet (str): Characters used in sequences.
        alphabet_size (int): Number of unique characters.
        diversity (float): Statistical diversity measure (0-1).
        uniform_lengths (bool): Whether all sequences have same length.
        file_path (Optional[str]): Filesystem path if persisted.
        file_size (Optional[str]): Human-readable file size.
        created_at (datetime): When dataset was created.
        generation_params (Optional[Dict]): Parameters used for generation.
    """

    id: str = Field(..., description="Dataset unique identifier")
    name: str = Field(..., description="Dataset name")
    type: DatasetType = Field(..., description="Dataset type")

    # Statistics
    size: int = Field(..., description="Number of sequences")
    min_length: int = Field(..., description="Minimum sequence length")
    max_length: int = Field(..., description="Maximum sequence length")
    average_length: float = Field(..., description="Average sequence length")
    alphabet: str = Field(..., description="Dataset alphabet")
    alphabet_size: int = Field(..., description="Size of alphabet")
    diversity: float = Field(..., description="Dataset diversity score")
    uniform_lengths: bool = Field(
        ..., description="Whether all sequences have same length"
    )

    # Metadata
    file_path: Optional[str] = Field(None, description="File path if persisted")
    file_size: Optional[str] = Field(None, description="Human readable file size")
    created_at: datetime = Field(..., description="Creation timestamp")
    generation_params: Optional[Dict[str, Any]] = Field(
        None, description="Parameters used for generation"
    )


class DatasetPreview(BaseModel):
    """Dataset preview with sample sequences.
    
    This model provides a preview of dataset contents including
    metadata and a sample of the actual sequences for user inspection.
    
    Attributes:
        info (DatasetInfo): Complete dataset metadata.
        sample_sequences (List[str]): First few sequences for preview.
        total_sequences (int): Total number of sequences in dataset.
    """

    info: DatasetInfo
    sample_sequences: List[str] = Field(
        ..., description="First few sequences for preview"
    )
    total_sequences: int = Field(..., description="Total number of sequences")


class DatasetUploadRequest(BaseModel):
    """Request for uploading dataset content.
    
    This model validates and handles user requests to upload
    FASTA-formatted sequence data to create new datasets.
    
    Attributes:
        name (str): Desired name for the new dataset.
        content (str): FASTA-formatted sequence content.
    """

    name: str = Field(..., description="Dataset name")
    content: str = Field(..., description="FASTA content")

    @field_validator("content")
    @classmethod
    def validate_fasta_content(cls, v: str):
        """Validate FASTA content format.
        
        Ensures the uploaded content follows valid FASTA format
        with proper headers and sequence data.
        
        Args:
            v (str): Content to validate.
            
        Returns:
            str: Validated content.
            
        Raises:
            ValueError: If content is empty or invalid FASTA format.
        """
        if not v.strip():
            raise ValueError("Content cannot be empty")

        lines = v.strip().split("\n")
        has_header = any(line.startswith(">") for line in lines)

        if not has_header:
            raise ValueError(
                "Content must be in FASTA format with headers starting with '>'"
            )

        return v


class SyntheticDatasetRequest(BaseModel):
    """Request for synthetic dataset generation.
    
    This model handles requests to generate synthetic biological
    sequence datasets using various algorithms and parameters.
    
    Attributes:
        name (str): Name for the generated dataset.
        method (SyntheticMethod): Algorithm to use for generation.
        n (int): Number of sequences to generate (3-1000).
        alphabet (str): Character set to use in sequences.
        seed (Optional[int]): Random seed for reproducible results.
        length (Optional[int]): Target sequence length (5-10000).
        base_sequence (Optional[str]): Template sequence for mutations/noise.
        noise_rate (Optional[float]): Probability of noise (0.0-1.0).
        mutation_rate (Optional[float]): Probability of mutations (0.0-1.0).
        mutation_types (Optional[List[str]]): Types of mutations to apply.
        num_clusters (Optional[int]): Number of sequence clusters.
        cluster_distance (Optional[float]): Distance between clusters (0.0-1.0).
        pad_char (Optional[str]): Character for sequence padding.
    """

    name: str = Field(..., description="Dataset name")
    method: SyntheticMethod = Field(..., description="Generation method")

    # Common parameters
    n: int = Field(20, ge=3, le=1000, description="Number of sequences")
    alphabet: str = Field("ACGT", description="Alphabet to use")
    seed: Optional[int] = Field(None, description="Random seed for reproducibility")

    # Method-specific parameters
    length: Optional[int] = Field(
        50, ge=5, le=10000, description="Sequence length (for random/clustered)"
    )
    base_sequence: Optional[str] = Field(
        None, description="Base sequence (for noise/mutations)"
    )
    noise_rate: Optional[float] = Field(
        0.1, ge=0.0, le=1.0, description="Noise rate (for noise method)"
    )
    mutation_rate: Optional[float] = Field(
        0.1, ge=0.0, le=1.0, description="Mutation rate (for mutations method)"
    )
    mutation_types: Optional[List[str]] = Field(
        ["substitution"], description="Types of mutations"
    )
    num_clusters: Optional[int] = Field(
        2, ge=1, description="Number of clusters (for clustered method)"
    )
    cluster_distance: Optional[float] = Field(
        0.2, ge=0.0, le=1.0, description="Distance between clusters"
    )
    pad_char: Optional[str] = Field("N", description="Padding character")


class NCBIDatasetRequest(BaseModel):
    """Request for NCBI dataset download.
    
    This model handles requests to download biological sequences
    from NCBI databases with filtering and processing options.
    
    Attributes:
        name (str): Name for the downloaded dataset.
        query (str): NCBI search query string.
        db (str): NCBI database to search (default: "nucleotide").
        max_sequences (Optional[int]): Maximum sequences to download (1-10000).
        min_length (Optional[int]): Minimum sequence length filter.
        max_length (Optional[int]): Maximum sequence length filter.
        uniform_policy (Optional[str]): Length uniformization strategy.
    """

    name: str = Field(..., description="Dataset name")
    query: str = Field(..., description="NCBI search query")
    db: str = Field("nucleotide", description="NCBI database")
    max_sequences: Optional[int] = Field(
        None, ge=1, le=10000, description="Maximum sequences to download"
    )
    min_length: Optional[int] = Field(
        None, ge=1, description="Minimum sequence length filter"
    )
    max_length: Optional[int] = Field(
        None, ge=1, description="Maximum sequence length filter"
    )
    uniform_policy: Optional[str] = Field(
        None, description="Uniformization policy (strict, pad, trim)"
    )

    @field_validator("max_sequences")
    @classmethod
    def validate_max_sequences(cls, v):
        """Validate max_sequences limit.
        
        Ensures the maximum sequences requested doesn't exceed
        system limits for performance and resource management.
        
        Args:
            v: Value to validate.
            
        Returns:
            Validated value.
            
        Raises:
            ValueError: If value exceeds maximum allowed limit.
        """
        if v is not None and v > 10000:
            raise ValueError(
                "Maximum sequences cannot exceed 10000 for performance reasons"
            )
        return v


class DatasetListResponse(BaseModel):
    """Response for dataset listing.
    
    This model provides a paginated list of datasets with
    metadata for display and navigation purposes.
    
    Attributes:
        datasets (List[DatasetInfo]): List of dataset information.
        total (int): Total number of datasets available.
    """

    datasets: List[DatasetInfo]
    total: int = Field(..., description="Total number of datasets")


class DatasetUpdateRequest(BaseModel):
    """Request for updating dataset metadata.
    
    This model handles requests to modify dataset properties
    such as name and other mutable metadata fields.
    
    Attributes:
        name (Optional[str]): New name for the dataset.
    """

    name: Optional[str] = Field(None, description="New dataset name")


class OperationResponse(BaseModel):
    """Generic operation response.
    
    This model provides a standardized response format for
    various dataset operations including success status,
    messages, and optional additional data.
    
    Attributes:
        success (bool): Whether the operation completed successfully.
        message (str): Human-readable result description.
        data (Optional[Dict]): Additional response data if needed.
    """

    success: bool = Field(..., description="Whether operation succeeded")
    message: str = Field(..., description="Operation result message")
    data: Optional[Dict[str, Any]] = Field(None, description="Additional response data")


class DatasetGenerationStatus(BaseModel):
    """Dataset generation status response.
    
    This model tracks the progress and status of dataset generation
    operations, particularly for long-running tasks like NCBI downloads.
    
    Attributes:
        status (str): Current status (pending, running, completed, failed).
        progress (int): Completion percentage (0-100).
        message (str): Current status message.
        dataset_id (Optional[str]): Generated dataset ID when completed.
        error (Optional[str]): Error description if generation failed.
    """

    status: str = Field(
        ..., description="Generation status (pending, running, completed, failed)"
    )
    progress: int = Field(0, ge=0, le=100, description="Progress percentage")
    message: str = Field("", description="Status message")
    dataset_id: Optional[str] = Field(
        None, description="Generated dataset ID if completed"
    )
    error: Optional[str] = Field(None, description="Error message if failed")
