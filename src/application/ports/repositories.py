"""
Application Ports - Infrastructure Interfaces.

This module defines the interfaces (ports) that must be implemented
by the infrastructure layer following the hexagonal architecture pattern.

These ports decouple the application layer from specific infrastructure
implementations, allowing for easy testing and component replacement.

Key Interfaces:
    - AlgorithmRegistry: Algorithm management and discovery
    - ExportPort: Result export and formatting
    - ExecutionEngine: Task execution orchestration

The interfaces defined here serve as contracts between the application
layer and infrastructure implementations, ensuring loose coupling and
high testability.
"""

from abc import ABC, abstractmethod
from typing import (
    Any,
    Dict,
    List,
    Protocol,
    runtime_checkable,
)

from src.domain import CSPAlgorithm, Dataset
from src.domain.config import (
    AlgParams,
    CSPBenchConfig,
    TaskConfig,
)
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)


class AlgorithmRegistry(Protocol):
    """
    Port for algorithm registry.
    
    Defines the interface for algorithm management, including registration,
    discovery, and metadata retrieval. This port abstracts the algorithm
    storage and retrieval mechanism from the application layer.
    """

    def get_algorithm(self, name: str) -> type[CSPAlgorithm]:
        """
        Get algorithm class by name.

        Args:
            name (str): Algorithm name.

        Returns:
            type[CSPAlgorithm]: Algorithm class.
            
        Raises:
            KeyError: If algorithm not found.
        """
        ...

    def list_algorithms(self) -> List[str]:
        """
        List available algorithms.

        Returns:
            List[str]: Available algorithm names.
        """
        ...

    def register_algorithm(self, algorithm_class: type[CSPAlgorithm]) -> None:
        """
        Register new algorithm.

        Args:
            algorithm_class (type[CSPAlgorithm]): Algorithm class to register.
        """
        ...

    def algorithm_exists(self, name: str) -> bool:
        """
        Check if algorithm exists.

        Args:
            name (str): Algorithm name.

        Returns:
            bool: True if exists, False otherwise.
        """
        ...

    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """
        Get algorithm metadata.

        Args:
            name (str): Algorithm name.

        Returns:
            Dict[str, Any]: Algorithm metadata including description, parameters, etc.
        """
        ...


@runtime_checkable
class ExportPort(Protocol):
    """
    Port for result export.
    
    Defines the interface for exporting results in various formats
    and managing export destinations. This port abstracts the export
    mechanism from the application layer.
    """

    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """
        Export results in specific format.

        Args:
            results (Dict[str, Any]): Result data to export.
            format_type (str): Export format (csv, json, xlsx, etc.).
            destination (str): Export destination path.

        Returns:
            str: Path of exported file.
        """
        ...

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """
        Export batch results.

        Args:
            batch_results (List[Dict[str, Any]]): List of results to export.
            format_type (str): Export format.
            destination (str): Export destination path.

        Returns:
            str: Path of exported file.
        """
        ...

    def get_supported_formats(self) -> List[str]:
        """
        List supported formats.

        Returns:
            List[str]: Available export formats.
        """
        ...

    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """
        Export optimization results.

        Args:
            optimization_data (Dict[str, Any]): Optimization data to export.
            destination (str): Export destination path.

        Returns:
            str: Path of exported file.
        """
        ...


@runtime_checkable
class ExecutionEngine(Protocol):
    """
    Port for pipeline task execution engines.
    
    Defines the interface for executing individual tasks within
    the pipeline execution framework. This port abstracts the
    task execution mechanism from the application layer.
    """

    def __init__(
        self,
        combination_store: CombinationScopedPersistence,
        execution_controller: ExecutionController,
        batch_config: CSPBenchConfig,
    ):
        """
        Initialize execution engine.
        
        Args:
            combination_store (CombinationScopedPersistence): Store for combination data.
            execution_controller (ExecutionController): Execution control interface.
            batch_config (CSPBenchConfig): Batch configuration.
        """
        ...

    def run(
        self,
        task: TaskConfig,
        dataset_obj: Dataset,
        alg: AlgParams,
    ) -> BaseStatus:
        """
        Execute a specific task.

        Args:
            task (TaskConfig): Task configuration to execute.
            dataset_obj (Dataset): Dataset object for the task.
            alg (AlgParams): Algorithm parameters.

        Returns:
            BaseStatus: Execution status result.
        """
        ...


class AbstractAlgorithmRegistry(ABC):
    """
    Abstract base class for algorithm registry.
    
    Provides a concrete base class that can be inherited by
    infrastructure implementations of the AlgorithmRegistry port.
    """

    @abstractmethod
    def get_algorithm(self, name: str) -> type[CSPAlgorithm]:
        """
        Get algorithm class by name.
        
        Args:
            name (str): Algorithm name.
            
        Returns:
            type[CSPAlgorithm]: Algorithm class.
        """
        pass

    @abstractmethod
    def list_algorithms(self) -> List[str]:
        """
        List available algorithms.
        
        Returns:
            List[str]: Available algorithm names.
        """
        pass

    @abstractmethod
    def register_algorithm(self, algorithm_class: type[CSPAlgorithm]) -> None:
        """
        Register new algorithm.
        
        Args:
            algorithm_class (type[CSPAlgorithm]): Algorithm class to register.
        """
        pass

    @abstractmethod
    def algorithm_exists(self, name: str) -> bool:
        """
        Check if algorithm exists.
        
        Args:
            name (str): Algorithm name.
            
        Returns:
            bool: True if exists, False otherwise.
        """
        pass

    @abstractmethod
    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """
        Get algorithm metadata.
        
        Args:
            name (str): Algorithm name.
            
        Returns:
            Dict[str, Any]: Algorithm metadata.
        """
        pass


class AbstractExportPort(ABC):
    """
    Abstract base class for result export.
    
    Provides a concrete base class that can be inherited by
    infrastructure implementations of the ExportPort.
    """

    @abstractmethod
    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """
        Export results in specific format.
        
        Args:
            results (Dict[str, Any]): Result data to export.
            format_type (str): Export format.
            destination (str): Export destination path.
            
        Returns:
            str: Path of exported file.
        """
        pass

    @abstractmethod
    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """
        Export batch results.
        
        Args:
            batch_results (List[Dict[str, Any]]): List of results to export.
            format_type (str): Export format.
            destination (str): Export destination path.
            
        Returns:
            str: Path of exported file.
        """
        pass

    @abstractmethod
    def get_supported_formats(self) -> List[str]:
        """
        List supported formats.
        
        Returns:
            List[str]: Available export formats.
        """
        pass

    @abstractmethod
    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """
        Export optimization results.
        
        Args:
            optimization_data (Dict[str, Any]): Optimization data to export.
            destination (str): Export destination path.
            
        Returns:
            str: Path of exported file.
        """
        pass


class AbstractExecutionEngine(ABC):
    """
    Abstract base class for task execution engines.
    
    Provides a concrete base class that can be inherited by
    infrastructure implementations of the ExecutionEngine port.
    """

    def __init__(
        self,
        combination_store: CombinationScopedPersistence,
        execution_controller: ExecutionController,
        batch_config: CSPBenchConfig,
    ):
        """
        Initialize execution engine.
        
        Args:
            combination_store (CombinationScopedPersistence): Store for combination data.
            execution_controller (ExecutionController): Execution control interface.
            batch_config (CSPBenchConfig): Batch configuration.
        """
        self._combination_store = combination_store
        self._execution_controller = execution_controller
        self._batch_config = batch_config

    @abstractmethod
    def run(
        self,
        task: TaskConfig,
        dataset_obj: Dataset,
        alg: AlgParams,
    ) -> BaseStatus:
        """
        Execute a specific task.
        
        Args:
            task (TaskConfig): Task configuration to execute.
            dataset_obj (Dataset): Dataset object for the task.
            alg (AlgParams): Algorithm parameters.
            
        Returns:
            BaseStatus: Execution status result.
        """
        pass
