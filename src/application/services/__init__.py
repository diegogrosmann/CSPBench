"""
Application Services

Implements use cases and application logic.
"""

__all__: list[str] = ["ExperimentService"]


def __getattr__(name):  # PEP 562 lazy import to avoid heavy deps during tests
    if name == "ExperimentService":
        from .experiment_service import ExperimentService as _ES

        return _ES
    raise AttributeError(name)
