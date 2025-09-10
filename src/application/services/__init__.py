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
"""

__all__: list[str] = ["PipelineService"]


def __getattr__(name):
    """
    PEP 562 lazy import to avoid heavy dependencies during tests.

    Args:
        name: Attribute name being accessed

    Returns:
        Requested service class

    Raises:
        AttributeError: If requested service doesn't exist
    """
    if name == "PipelineService":
        from .pipeline_service import PipelineService as _PS

        return _PS
    raise AttributeError(name)
