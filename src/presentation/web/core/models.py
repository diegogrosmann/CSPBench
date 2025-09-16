"""
Pydantic models for web API requests and responses.

This module defines the core data models used throughout the CSPBench
web interface for API communication and data validation.
"""

from typing import Dict, Optional

from pydantic import BaseModel


class AlgorithmInfo(BaseModel):
    """Information about an available algorithm.
    
    This model represents metadata and configuration information
    for algorithms available in the CSPBench system.
    
    Attributes:
        name (str): Algorithm identifier/name.
        description (str): Human-readable algorithm description.
        default_params (Dict): Default parameter values and their types.
        is_deterministic (bool): Whether algorithm produces consistent results.
        supports_internal_parallel (bool): Whether algorithm supports parallelization.
        category (Optional[str]): Algorithm category for organization.
    """

    name: str
    description: str
    default_params: Dict
    is_deterministic: bool
    supports_internal_parallel: bool
    category: Optional[str] = None


class HealthCheck(BaseModel):
    """Health check response model.
    
    This model represents the system health status including
    component availability and version information.
    
    Attributes:
        status (str): Overall system status (e.g., "healthy", "degraded").
        timestamp (str): ISO timestamp of the health check.
        components (Dict[str, bool]): Status of individual system components.
        version (str): Application version number.
    """

    status: str
    timestamp: str
    components: Dict[str, bool]
    version: str = "0.1.0"
