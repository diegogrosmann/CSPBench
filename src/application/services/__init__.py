"""
Application Services Module.

Implements use cases and application logic that coordinate between
domain entities and infrastructure services. This module provides
the main service layer for the CSPBench application.

Services:
    PipelineService: Core pipeline execution orchestration

Note:
    Uses PEP 562 lazy import to avoid heavy dependencies during tests.
    Services are loaded on-demand when first accessed.

Example:
    Using lazy import to access services:
    
    >>> from src.application.services import PipelineService
    >>> service = PipelineService()
"""

__all__: list[str] = ["PipelineService"]


def __getattr__(name):
    """
    PEP 562 lazy import to avoid heavy dependencies during tests.

    This function implements lazy loading of service modules to prevent
    importing heavy dependencies when they're not needed, particularly
    useful during testing.

    Args:
        name (str): Attribute name being accessed.

    Returns:
        type: Requested service class.

    Raises:
        AttributeError: If requested service doesn't exist.
    """
    if name == "PipelineService":
        from .pipeline_service import PipelineService as _PS

        return _PS
    raise AttributeError(name)
