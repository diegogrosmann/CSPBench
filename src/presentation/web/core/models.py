"""
Pydantic models for web API requests and responses.
"""

from typing import Dict, List, Optional

from pydantic import BaseModel


class AlgorithmInfo(BaseModel):
    """Information about an available algorithm."""

    name: str
    description: str
    default_params: Dict
    is_deterministic: bool
    supports_internal_parallel: bool
    category: Optional[str] = None


class HealthCheck(BaseModel):
    """Health check response."""

    status: str
    timestamp: str
    components: Dict[str, bool]
    version: str = "0.1.0"