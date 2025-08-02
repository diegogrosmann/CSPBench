"""
Pydantic models for web API requests and responses.
"""

from datetime import datetime
from typing import Dict, List, Optional

from pydantic import BaseModel, validator

from .security import SecurityValidator


class AlgorithmInfo(BaseModel):
    """Information about an available algorithm."""

    name: str
    description: str
    default_params: Dict
    is_deterministic: bool
    supports_internal_parallel: bool
    category: Optional[str] = None


class ExecutionRequest(BaseModel):
    """Request for algorithm execution."""

    algorithm: str
    dataset_content: Optional[str] = None
    dataset_name: Optional[str] = "uploaded_dataset.fasta"
    parameters: Dict = {}
    save_history: bool = False
    timeout: int = 300

    @validator("algorithm")
    def validate_algorithm(cls, v):
        if not v or not isinstance(v, str):
            raise ValueError("Algorithm name is required")
        return v

    @validator("dataset_name")
    def validate_dataset_name(cls, v):
        return SecurityValidator.sanitize_filename(v) if v else "uploaded_dataset.fasta"

    @validator("timeout")
    def validate_timeout(cls, v):
        if v < 10 or v > 3600:
            raise ValueError("Timeout must be between 10 and 3600 seconds")
        return v


class ExecutionResult(BaseModel):
    """Result of algorithm execution."""

    session_id: str
    status: str
    result: Optional[Dict] = None
    error: Optional[str] = None
    download_url: Optional[str] = None
    timestamp: str = datetime.now().isoformat()


class HealthCheck(BaseModel):
    """Health check response."""

    status: str
    timestamp: str
    components: Dict[str, bool]
    version: str = "0.1.0"


class DatasetGenerationRequest(BaseModel):
    """Request for synthetic dataset generation."""

    num_strings: int
    string_length: int
    alphabet: str = "ACGT"
    max_distance: Optional[int] = None
    generation_method: str = "random"
    seed: Optional[int] = None

    @validator("num_strings")
    def validate_num_strings(cls, v):
        if v < 1 or v > 10000:
            raise ValueError("Number of strings must be between 1 and 10000")
        return v

    @validator("string_length")
    def validate_string_length(cls, v):
        if v < 1 or v > 10000:
            raise ValueError("String length must be between 1 and 10000")
        return v

    @validator("alphabet")
    def validate_alphabet(cls, v):
        if not v or len(v) < 2:
            raise ValueError("Alphabet must have at least 2 characters")
        return v


class DatasetGenerationResult(BaseModel):
    """Result of dataset generation."""

    session_id: str
    status: str
    filename: Optional[str] = None
    download_url: Optional[str] = None
    error: Optional[str] = None
    metadata: Optional[Dict] = None
    dataset_info: Optional[Dict] = None
    sequences: Optional[List[str]] = None


class NCBIDatasetRequest(BaseModel):
    """Request for NCBI dataset generation."""

    query: str
    max_sequences: int = 100
    sequence_type: str = "nucleotide"
    email: str
    api_key: Optional[str] = None

    @validator("query")
    def validate_query(cls, v):
        if not v or len(v.strip()) < 3:
            raise ValueError("Query must be at least 3 characters long")
        return v.strip()

    @validator("max_sequences")
    def validate_max_sequences(cls, v):
        if v < 1 or v > 1000:
            raise ValueError("Max sequences must be between 1 and 1000")
        return v

    @validator("email")
    def validate_email(cls, v):
        import re

        if not re.match(r"^[^@]+@[^@]+\.[^@]+$", v):
            raise ValueError("Invalid email format")
        return v
