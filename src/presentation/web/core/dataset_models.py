"""
Dataset-related API models for CSPBench Web Interface.
"""

from enum import Enum
from typing import Any, Dict, List, Optional, Union
from datetime import datetime

from pydantic import BaseModel, Field, validator


class DatasetType(str, Enum):
    """Types of datasets."""
    SYNTHETIC = "synthetic"
    FILE = "file"  
    NCBI = "ncbi"


class SyntheticMethod(str, Enum):
    """Synthetic dataset generation methods."""
    RANDOM = "random"
    NOISE = "noise"
    CLUSTERED = "clustered"
    MUTATIONS = "mutations"


class DatasetInfo(BaseModel):
    """Dataset information response."""
    
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
    uniform_lengths: bool = Field(..., description="Whether all sequences have same length")
    
    # Metadata
    file_path: Optional[str] = Field(None, description="File path if persisted")
    file_size: Optional[str] = Field(None, description="Human readable file size")
    created_at: datetime = Field(..., description="Creation timestamp")
    generation_params: Optional[Dict[str, Any]] = Field(None, description="Parameters used for generation")


class DatasetPreview(BaseModel):
    """Dataset preview with sample sequences."""
    
    info: DatasetInfo
    sample_sequences: List[str] = Field(..., description="First few sequences for preview")
    total_sequences: int = Field(..., description="Total number of sequences")


class DatasetUploadRequest(BaseModel):
    """Request for uploading dataset content."""
    
    name: str = Field(..., description="Dataset name")
    content: str = Field(..., description="FASTA content")
    
    @validator('content')
    def validate_fasta_content(cls, v):
        """Validate FASTA content format."""
        if not v.strip():
            raise ValueError("Content cannot be empty")
        
        lines = v.strip().split('\n')
        has_header = any(line.startswith('>') for line in lines)
        
        if not has_header:
            raise ValueError("Content must be in FASTA format with headers starting with '>'")
        
        return v


class SyntheticDatasetRequest(BaseModel):
    """Request for synthetic dataset generation."""
    
    name: str = Field(..., description="Dataset name")
    method: SyntheticMethod = Field(..., description="Generation method")
    
    # Common parameters
    n: int = Field(20, ge=3, le=1000, description="Number of sequences")
    alphabet: str = Field("ACGT", description="Alphabet to use")
    seed: Optional[int] = Field(None, description="Random seed for reproducibility")
    
    # Method-specific parameters
    length: Optional[int] = Field(50, ge=5, le=10000, description="Sequence length (for random/clustered)")
    base_sequence: Optional[str] = Field(None, description="Base sequence (for noise/mutations)")
    noise_rate: Optional[float] = Field(0.1, ge=0.0, le=1.0, description="Noise rate (for noise method)")
    mutation_rate: Optional[float] = Field(0.1, ge=0.0, le=1.0, description="Mutation rate (for mutations method)")
    mutation_types: Optional[List[str]] = Field(["substitution"], description="Types of mutations")
    num_clusters: Optional[int] = Field(2, ge=1, description="Number of clusters (for clustered method)")
    cluster_distance: Optional[float] = Field(0.2, ge=0.0, le=1.0, description="Distance between clusters")
    pad_char: Optional[str] = Field("N", description="Padding character")


class NCBIDatasetRequest(BaseModel):
    """Request for NCBI dataset download."""
    
    name: str = Field(..., description="Dataset name")
    query: str = Field(..., description="NCBI search query")
    db: str = Field("nucleotide", description="NCBI database")
    max_sequences: Optional[int] = Field(None, ge=1, le=10000, description="Maximum sequences to download")
    min_length: Optional[int] = Field(None, ge=1, description="Minimum sequence length filter")
    max_length: Optional[int] = Field(None, ge=1, description="Maximum sequence length filter") 
    uniform_policy: Optional[str] = Field(None, description="Uniformization policy (strict, pad, trim)")
    
    @validator('max_sequences')
    def validate_max_sequences(cls, v):
        """Validate max_sequences limit."""
        if v is not None and v > 10000:
            raise ValueError("Maximum sequences cannot exceed 10000 for performance reasons")
        return v


class DatasetListResponse(BaseModel):
    """Response for dataset listing."""
    
    datasets: List[DatasetInfo]
    total: int = Field(..., description="Total number of datasets")


class DatasetUpdateRequest(BaseModel):
    """Request for updating dataset metadata."""
    
    name: Optional[str] = Field(None, description="New dataset name")


class OperationResponse(BaseModel):
    """Generic operation response."""
    
    success: bool = Field(..., description="Whether operation succeeded")
    message: str = Field(..., description="Operation result message")
    data: Optional[Dict[str, Any]] = Field(None, description="Additional response data")


class DatasetGenerationStatus(BaseModel):
    """Dataset generation status response."""
    
    status: str = Field(..., description="Generation status (pending, running, completed, failed)")
    progress: int = Field(0, ge=0, le=100, description="Progress percentage")
    message: str = Field("", description="Status message")
    dataset_id: Optional[str] = Field(None, description="Generated dataset ID if completed")
    error: Optional[str] = Field(None, description="Error message if failed")
