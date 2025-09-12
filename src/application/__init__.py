"""
CSPBench Application Layer.

This module contains the application layer of the CSPBench system,
implementing use cases and application logic that orchestrate
domain entities and infrastructure services.

The application layer acts as a mediator between the presentation
layer (CLI/Web interfaces) and the domain/infrastructure layers,
ensuring proper separation of concerns and clean architecture principles.

Key Components:
    - Work management services
    - Execution orchestration
    - Configuration parsing and validation
    - Dataset generation and loading
    - Service coordination

Note:
    __all__ is intentionally empty to avoid premature imports of
    legacy services during testing phases.
"""

__all__: list[str] = []  # Avoid early import of legacy services
