"""
Health check and monitoring endpoints.
"""

import logging
from datetime import datetime

from fastapi import APIRouter

from src.domain.algorithms import global_registry

from ..core.config import web_config
from ..core.models import HealthCheck

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/health", tags=["health"])


@router.get("/")
async def simple_health():
    """Simple health check endpoint for basic monitoring."""
    try:
        algorithms_loaded = len(global_registry)
        return {
            "status": "healthy",
            "timestamp": datetime.now().isoformat(),
            "algorithms_loaded": algorithms_loaded,
            "version": "1.0.0",
        }
    except Exception as e:
        logger.error(f"Simple health check failed: {e}")
        return {
            "status": "error",
            "timestamp": datetime.now().isoformat(),
            "error": str(e),
            "version": "1.0.0",
        }


@router.get("/detailed", response_model=HealthCheck)
async def detailed_health_check():
    """Detailed health check following security guidelines."""
    try:
        # Check essential components
        algorithms_available = len(global_registry) > 0
        config_loaded = web_config.config is not None
        session_manager_available = web_config.get_session_manager() is not None
        experiment_service_available = web_config.get_experiment_service() is not None

        status = (
            "healthy"
            if all(
                [
                    algorithms_available,
                    config_loaded,
                    session_manager_available,
                    experiment_service_available,
                ]
            )
            else "degraded"
        )

        return HealthCheck(
            status=status,
            timestamp=datetime.now().isoformat(),
            components={
                "algorithms": algorithms_available,
                "configuration": config_loaded,
                "session_manager": session_manager_available,
                "experiment_service": experiment_service_available,
            },
        )
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return HealthCheck(
            status="degraded",
            timestamp=datetime.now().isoformat(),
            components={
                "algorithms": False,
                "configuration": False,
                "session_manager": False,
                "experiment_service": False,
            },
        )


@router.get("/metrics")
async def get_metrics():
    """Basic metrics endpoint (can be extended with Prometheus)."""
    try:
        return {
            "algorithms_count": len(global_registry),
            "status": "running",
            "uptime": "unknown",  # Would need startup tracking
            "timestamp": datetime.now().isoformat(),
        }
    except Exception as e:
        logger.error(f"Metrics collection failed: {e}")
        return {"error": "Metrics unavailable", "timestamp": datetime.now().isoformat()}
