"""
Health check and monitoring endpoints.
"""

import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from src.domain.algorithms import global_registry
from ..core.models import HealthCheck

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/health", tags=["health"])


class HealthStatus(BaseModel):
    """Health status model."""
    status: str
    version: str = "0.1.0"


def get_version():
    """Get application version from multiple sources."""
    try:
        # Try to get version from pyproject.toml
        project_root = Path(__file__).parent.parent.parent.parent.parent
        pyproject_path = project_root / "pyproject.toml"
        
        if pyproject_path.exists():
            with open(pyproject_path, 'r') as f:
                content = f.read()
                for line in content.split('\n'):
                    if line.strip().startswith('version ='):
                        version = line.split('=')[1].strip().strip('"\'')
                        return f"v{version}"
        
        # Fallback to environment variable or default
        return os.getenv("CSPBENCH_VERSION", "v0.1.0")
    except Exception:
        return "v0.1.0"


@router.get("/")
async def health_check():
    """Basic health check."""
    try:
        algorithms_loaded = len(global_registry)
        version = get_version()
        
        return HealthStatus(status="healthy", version=version)
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return HealthStatus(status="error", version=get_version())


@router.get("/detailed", response_model=HealthCheck)
async def detailed_health_check():
    """Detailed health check following security guidelines."""
    try:
        # Check essential components
        algorithms_available = len(global_registry) > 0

        status = "healthy" if algorithms_available else "degraded"

        return HealthCheck(
            status=status,
            algorithms=len(global_registry),
            config_loaded=True,
        )

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail={
                "status": "unhealthy",
                "algorithms": 0,
                "config_loaded": False,
                "error": str(e),
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


@router.get("/status")
async def get_system_status():
    """Get system status for UI display."""
    try:
        algorithms_count = len(global_registry)
        
        # Determine system status based on available components
        if algorithms_count > 0:
            status = "Online"
            status_type = "success"
        else:
            status = "Degraded"
            status_type = "warning"
            
        return {
            "status": status,
            "status_type": status_type,
            "algorithms_available": algorithms_count > 0,
            "timestamp": datetime.now().isoformat(),
        }
    except Exception as e:
        logger.error(f"Status check failed: {e}")
        return {
            "status": "Error",
            "status_type": "danger",
            "algorithms_available": False,
            "timestamp": datetime.now().isoformat(),
            "error": str(e),
        }


@router.get("/version")
async def get_system_version():
    """Get system version information."""
    try:
        version = get_version()
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
        
        return {
            "version": version,
            "python_version": python_version,
            "architecture": "Hexagonal",
            "timestamp": datetime.now().isoformat(),
        }
    except Exception as e:
        logger.error(f"Version check failed: {e}")
        return {
            "version": "Unknown",
            "python_version": "Unknown",
            "architecture": "Hexagonal",
            "timestamp": datetime.now().isoformat(),
            "error": str(e),
        }
