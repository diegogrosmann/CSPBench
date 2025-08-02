"""
Algorithm-related API endpoints.
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
    """Get list of available algorithms with metadata."""
    try:
        algorithms = []

        for name, algorithm_cls in global_registry.items():
            try:
                # Create temporary instance to get metadata
                temp_instance = algorithm_cls()

                # Get default parameters
                default_params = getattr(temp_instance, "default_params", {})

                # Get metadata
                is_deterministic = getattr(temp_instance, "is_deterministic", True)
                supports_parallel = getattr(
                    temp_instance, "supports_internal_parallel", False
                )
                category = getattr(temp_instance, "category", "General")

                algorithm_info = AlgorithmInfo(
                    name=name,
                    description=algorithm_cls.__doc__ or f"{name} algorithm",
                    default_params=default_params,
                    is_deterministic=is_deterministic,
                    supports_internal_parallel=supports_parallel,
                    category=category,
                )

                algorithms.append(algorithm_info)

            except Exception as e:
                logger.warning(f"Error getting metadata for algorithm {name}: {e}")
                # Add basic info even if metadata extraction fails
                algorithm_info = AlgorithmInfo(
                    name=name,
                    description=algorithm_cls.__doc__ or f"{name} algorithm",
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
    """Get detailed information about a specific algorithm."""
    try:
        algorithm_cls = global_registry.get(algorithm_name)
        if not algorithm_cls:
            raise HTTPException(
                status_code=404, detail=f"Algorithm '{algorithm_name}' not found"
            )

        # Create temporary instance to get metadata
        temp_instance = algorithm_cls()

        return {
            "name": algorithm_name,
            "description": algorithm_cls.__doc__ or f"{algorithm_name} algorithm",
            "default_params": getattr(temp_instance, "default_params", {}),
            "is_deterministic": getattr(temp_instance, "is_deterministic", True),
            "supports_internal_parallel": getattr(
                temp_instance, "supports_internal_parallel", False
            ),
            "category": getattr(temp_instance, "category", "General"),
            "version": getattr(algorithm_cls, "__version__", "1.0.0"),
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting algorithm info for {algorithm_name}: {e}")
        raise HTTPException(
            status_code=500, detail="Failed to retrieve algorithm information"
        )
