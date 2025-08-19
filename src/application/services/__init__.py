"""
Application Services

Implements use cases and application logic.
"""

__all__: list[str] = ["PipelineService"]


def __getattr__(name):  # PEP 562 lazy import to avoid heavy deps during tests
    if name == "PipelineService":
        from .pipeline_service import PipelineService as _PS

        return _PS
    raise AttributeError(name)
