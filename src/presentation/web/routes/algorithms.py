"""
Algorithm-related API endpoints for CSPBench Web Interface.

This module provides REST API endpoints for algorithm discovery, metadata
retrieval, and parameter information. It integrates with the algorithm
registry to provide comprehensive algorithm information for the web interface.
"""

import logging
from typing import List

from fastapi import APIRouter, HTTPException

from src.domain.algorithms import global_registry

from ..core.models import AlgorithmInfo

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/algorithms", tags=["algorithms"])


@router.get("/", response_model=List[AlgorithmInfo])
async def get_algorithms():
    """Get list of available algorithms with comprehensive metadata.

    Retrieves all registered algorithms with their metadata including
    parameters, capabilities, and categorization information for
    display in the web interface.

    Returns:
        List[AlgorithmInfo]: Complete list of available algorithms with metadata.

    Raises:
        HTTPException: If algorithm registry access fails.

    Note:
        Gracefully handles metadata extraction errors for individual algorithms,
        providing basic information even when detailed metadata is unavailable.
        This ensures the interface remains functional even with problematic algorithms.
    """
    try:
        algorithms = []

        for name, algorithm_cls in global_registry.items():
            try:
                # Get metadata from class attributes instead of creating instance
                # This avoids the issue with required constructor arguments

                # Get default parameters from class attribute
                default_params = getattr(algorithm_cls, "default_params", {})

                # Get metadata from class attributes
                is_deterministic = getattr(algorithm_cls, "is_deterministic", True)
                supports_parallel = getattr(
                    algorithm_cls, "supports_internal_parallel", False
                )
                category = getattr(algorithm_cls, "category", "General")

                # Safely get description
                try:
                    description = algorithm_cls.__doc__ or f"{name} algorithm"
                except Exception:
                    description = f"{name} algorithm"

                algorithm_info = AlgorithmInfo(
                    name=name,
                    description=description,
                    default_params=default_params,
                    is_deterministic=is_deterministic,
                    supports_internal_parallel=supports_parallel,
                    category=category,
                )

                algorithms.append(algorithm_info)

            except Exception as e:
                logger.warning(f"Error getting metadata for algorithm {name}: {e}")
                # Add basic info even if metadata extraction fails
                # Use try/except for __doc__ access in case it also fails
                try:
                    description = algorithm_cls.__doc__ or f"{name} algorithm"
                except Exception:
                    description = f"{name} algorithm"

                algorithm_info = AlgorithmInfo(
                    name=name,
                    description=description,
                    default_params={},
                    is_deterministic=True,
                    supports_internal_parallel=False,
                    category="General",
                )
                algorithms.append(algorithm_info)

        return algorithms

    except Exception as e:
        logger.error(f"Error getting algorithms: {e}")
        raise HTTPException(status_code=500, detail="Failed to retrieve algorithms")


@router.get("/{algorithm_name}")
async def get_algorithm_info(algorithm_name: str):
    """Get detailed information about a specific algorithm.

    Retrieves comprehensive metadata and configuration information
    for a single algorithm including parameters, capabilities,
    and version information.

    Args:
        algorithm_name (str): Name/identifier of the algorithm to query.

    Returns:
        dict: Detailed algorithm information including all available metadata.

    Raises:
        HTTPException: If algorithm not found or information retrieval fails.

    Note:
        Provides extended information compared to the list endpoint,
        including version information and detailed parameter schemas
        for algorithm configuration interfaces.
    """
    try:
        algorithm_cls = global_registry.get(algorithm_name)
        if not algorithm_cls:
            raise HTTPException(
                status_code=404, detail=f"Algorithm '{algorithm_name}' not found"
            )

        # Get metadata from class attributes instead of creating instance
        # This avoids the issue with required constructor arguments

        return {
            "name": algorithm_name,
            "description": algorithm_cls.__doc__ or f"{algorithm_name} algorithm",
            "default_params": getattr(algorithm_cls, "default_params", {}),
            "is_deterministic": getattr(algorithm_cls, "is_deterministic", True),
            "supports_internal_parallel": getattr(
                algorithm_cls, "supports_internal_parallel", False
            ),
            "category": getattr(algorithm_cls, "category", "General"),
            "version": getattr(algorithm_cls, "__version__", "1.0.0"),
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting algorithm info for {algorithm_name}: {e}")
        raise HTTPException(
            status_code=500, detail="Failed to retrieve algorithm information"
        )
