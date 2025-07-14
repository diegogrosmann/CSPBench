"""Interfaces e tipos para sistema de monitoramento."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Dict, Any, Optional, List


class TaskType(Enum):
    """Tipos de tarefa suportados."""
    EXECUTION = "execution"
    OPTIMIZATION = "optimization"
    SENSITIVITY = "sensitivity"


@dataclass
class ExecutionData:
    """Dados específicos de execução."""
    current_execution: str = ""
    total_executions: int = 1
    completed_executions: int = 0
    current_algorithm: str = ""
    total_algorithms: int = 0
    completed_algorithms: int = 0
    algorithm_progress: Dict[str, float] = field(default_factory=dict)
    algorithm_results: Dict[str, Any] = field(default_factory=dict)
    algorithm_callback_info: Dict[str, str] = field(default_factory=dict)
    best_distance: Optional[int] = None
    current_dataset: str = ""
    total_datasets: int = 1
    current_dataset_index: int = 1
    current_task_info: str = ""


@dataclass
class OptimizationData:
    """Dados específicos de otimização."""
    current_optimization: str = ""
    total_optimizations: int = 1
    completed_optimizations: int = 0
    study_name: str = ""
    n_trials: int = 0
    completed_trials: int = 0
    current_trial: Optional[int] = None
    best_value: Optional[float] = None
    best_params: Optional[Dict[str, Any]] = None
    trial_callback_info: str = ""
    current_dataset: str = ""
    total_datasets: int = 1
    current_dataset_index: int = 1
    current_task_info: str = ""
    sampler_name: str = "TPESampler"


@dataclass
class SensitivityData:
    """Dados específicos de análise de sensibilidade."""
    current_analysis: str = ""
    total_analyses: int = 1
    completed_analyses: int = 0
    analysis_method: str = ""
    n_samples: int = 0
    completed_samples: int = 0
    parameters: List[str] = field(default_factory=list)
    current_parameter: str = ""
    parameter_progress: Dict[str, float] = field(default_factory=dict)
    sensitivity_results: Dict[str, Dict[str, float]] = field(default_factory=dict)
    callback_info: str = ""
    current_dataset: str = ""
    total_datasets: int = 1
    current_dataset_index: int = 1
    current_task_info: str = ""
    method_details: str = ""


@dataclass
class MonitoringData:
    """Dados gerais de monitoramento."""
    task_type: TaskType
    batch_name: str
    start_time: datetime
    current_time: datetime
    session_id: str = ""
    
    # Específico por tipo
    execution_data: Optional[ExecutionData] = None
    optimization_data: Optional[OptimizationData] = None
    sensitivity_data: Optional[SensitivityData] = None


class MonitoringInterface(ABC):
    """Interface comum para sistemas de monitoramento."""
    
    @abstractmethod
    def start_monitoring(self, task_type: TaskType, config: Dict[str, Any]) -> None:
        """Inicia o monitoramento para um tipo de tarefa."""
        pass
    
    @abstractmethod
    def update_progress(self, progress_data: MonitoringData) -> None:
        """Atualiza informações de progresso."""
        pass
    
    @abstractmethod
    def finish_monitoring(self, results: Dict[str, Any]) -> None:
        """Finaliza o monitoramento."""
        pass
    
    @abstractmethod
    def show_error(self, error: str) -> None:
        """Exibe erro."""
        pass
    
    @abstractmethod
    def close(self) -> None:
        """Fecha e limpa recursos do monitor."""
        pass
