"""
Health check and monitoring endpoints.

This module provides health check endpoints for monitoring the CSPBench
web application status, component availability, and system metrics.
It supports both basic and detailed health checks for different monitoring needs.
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
    """Health status model.

    Simple health status response for basic health checks.

    Attributes:
        status (str): Overall system health status.
        version (str): Application version information.
    """

    status: str
    version: str = "0.1.0"


def get_version():
    """Get application version from multiple sources.

    Attempts to retrieve version information from various sources
    in order of preference: pyproject.toml, environment variable,
    or fallback to default version.

    Returns:
        str: Version string with 'v' prefix.

    Note:
        This function gracefully handles errors and always returns
        a valid version string, defaulting to "v0.1.0" if needed.
    """
    try:
        # Try to get version from pyproject.toml
        project_root = Path(__file__).parent.parent.parent.parent.parent
        pyproject_path = project_root / "pyproject.toml"

        if pyproject_path.exists():
            with open(pyproject_path) as f:
                content = f.read()
                for line in content.split("\n"):
                    if line.strip().startswith("version ="):
                        version = line.split("=")[1].strip().strip("\"'")
                        if version:  # Ensure version is not empty
                            return f"v{version}"

        # Fallback to environment variable or default
        env_version = os.getenv("CSPBENCH_VERSION")
        if env_version:
            return env_version

        return "v0.1.0"
    except Exception:
        return "v0.1.0"


@router.get("/")
async def health_check():
    """Basic health check endpoint.

    Provides a simple health status check with minimal resource usage.
    Suitable for load balancer health checks and basic monitoring.

    Returns:
        HealthStatus: Basic health status with version information.

    Note:
        This endpoint always returns a response, using "error" status
        if system components are unavailable rather than raising exceptions.
    """
    try:
        algorithms_loaded = len(global_registry)
        version = get_version()

        return HealthStatus(status="healthy", version=version)
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return HealthStatus(status="error", version=get_version())


@router.get("/detailed", response_model=HealthCheck)
async def detailed_health_check():
    """Detailed health check following security guidelines.

    Provides comprehensive health information including component
    status checks and detailed system information for monitoring
    and debugging purposes.

    Returns:
        HealthCheck: Detailed health status with component breakdown.

    Raises:
        HTTPException: If critical system components fail completely.
    """
    try:
        # Check essential components
        algorithms_available = len(global_registry) > 0

        status = "healthy" if algorithms_available else "degraded"

        return HealthCheck(
            status=status,
            timestamp=datetime.now().isoformat(),
            components={
                "algorithms": algorithms_available,
                "config": True,  # Assume config is loaded if we get here
            },
            version=get_version(),
        )

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail={
                "status": "unhealthy",
                "timestamp": datetime.now().isoformat(),
                "components": {
                    "algorithms": False,
                    "config": False,
                },
                "version": get_version(),
                "error": str(e),
            },
        )


@router.get("/metrics")
async def get_metrics():
    """Basic metrics endpoint for monitoring systems.

    Provides system metrics in a format suitable for monitoring
    systems like Prometheus. Can be extended with additional
    performance and operational metrics as needed.

    Returns:
        dict: System metrics including component counts and status.

    Note:
        This endpoint gracefully handles errors and returns
        error information rather than failing completely.
    """
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
    """Get system status for UI display.

    Provides system status information specifically formatted
    for display in user interfaces with appropriate status
    types and user-friendly messages.

    Returns:
        dict: UI-friendly system status with display information.

    Note:
        Status types (success, warning, danger) are designed
        for Bootstrap or similar UI framework integration.
    """
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
    """Get system version information.

    Provides comprehensive version information including
    application version, Python runtime, and architecture
    details for debugging and support purposes.

    Returns:
        dict: Detailed version and system information.

    Note:
        This endpoint is useful for support and debugging
        to identify exact system configuration and versions.
    """
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
