"""
Infrastructure Orchestrators Module

Implementations for algorithm execution and coordination.
"""

from .executors import Executor

# Backward compatibility shim for legacy tests/imports
LocalExecutor = Executor
from .optimization_orchestrator import OptimizationOrchestrator
from .optimization_report_generator import OptimizationReportGenerator
from .sensitivity_orchestrator import SensitivityOrchestrator

__all__ = [
    "Executor",
    "LocalExecutor",
    "OptimizationOrchestrator",
    "OptimizationReportGenerator",
    "SensitivityOrchestrator",
]
