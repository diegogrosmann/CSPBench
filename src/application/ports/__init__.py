"""
Application Ports Module.

Defines interfaces (ports) that must be implemented by infrastructure
layer following the hexagonal architecture pattern. These ports
provide clean separation between application logic and infrastructure
concerns.

Key Interfaces:
- AlgorithmRegistry: Algorithm management and discovery
- ExportPort: Result export and formatting
- ExecutorPort: Execution orchestration
- ExecutorInterface: Execution control interface

These interfaces allow the application layer to remain independent
of specific infrastructure implementations, facilitating testing
and component replacement.
"""

from .repositories import (
    AlgorithmRegistry,
    ExportPort,
)

__all__ = [
    "AlgorithmRegistry",
    "ExportPort",
    "ExecutorPort",
    "ExecutorInterface",
]
