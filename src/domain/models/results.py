"""
CSPBench Result Models

Domain models for execution results and outcomes.
These models store the results of algorithm executions, optimizations, and sensitivity analyses.
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional

from .enums import BatchType


@dataclass
class AlgorithmResult:
    """
    Result from a single algorithm execution.
    
    Contains all information about an algorithm run including performance metrics.
    """
    algorithm_name: str
    dataset_id: str
    repetition: int
    success: bool
    execution_time: float
    memory_usage_mb: Optional[float] = None
    solution_quality: Optional[Dict[str, Any]] = None
    solution_value: Optional[float] = None
    solution_string: Optional[str] = None
    iterations: Optional[int] = None
    convergence_info: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class OptimizationResult:
    """
    Result from a single optimization trial.
    
    Contains information about hyperparameter optimization trials.
    """
    trial_number: int
    algorithm_name: str
    dataset_id: str
    parameters: Dict[str, Any]
    objective_value: float
    success: bool
    execution_time: float
    trial_metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SensitivityResult:
    """
    Result from sensitivity analysis.
    
    Contains sensitivity indices and analysis results.
    """
    method: str
    algorithm_name: str
    dataset_id: str
    parameters_analyzed: List[str]
    sensitivity_indices: Dict[str, Any]
    samples_used: int
    confidence_intervals: Optional[Dict[str, Any]] = None
    analysis_metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class BatchResults:
    """
    Consolidated results from a complete batch execution.
    
    Root results object containing all execution outcomes and statistics.
    """
    batch_name: str
    batch_type: BatchType
    session_id: str
    start_time: datetime
    end_time: datetime
    total_duration_seconds: float
    successful_runs: int
    failed_runs: int
    total_runs: int
    algorithm_results: List[AlgorithmResult] = field(default_factory=list)
    optimization_results: List[OptimizationResult] = field(default_factory=list)
    sensitivity_results: List[SensitivityResult] = field(default_factory=list)
    summary_statistics: Dict[str, Any] = field(default_factory=dict)
    export_paths: Dict[str, str] = field(default_factory=dict)
    resource_usage: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
