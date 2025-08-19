"""
Portas da Aplicação - Interfaces para Infraestrutura

Este módulo define as interfaces (ports) que devem ser implementadas
pela camada de infraestrutura seguindo o padrão de arquitetura hexagonal.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Protocol, runtime_checkable, Callable

from src.domain import CSPAlgorithm, Dataset
from src.domain.config import AlgParams, ResourcesConfig, SystemConfig
from src.infrastructure.monitoring.monitor_interface import Monitor


class AlgorithmRegistry(Protocol):
    """Port para registry de algoritmos."""

    def get_algorithm(self, name: str) -> type[CSPAlgorithm]:
        """
        Obtém classe de algoritmo por nome.

        Args:
            name: Nome do algoritmo

        Returns:
            type[CSPAlgorithm]: Classe do algoritmo
        """
        ...

    def list_algorithms(self) -> List[str]:
        """
        Lista algoritmos disponíveis.

        Returns:
            List[str]: Nomes dos algoritmos disponíveis
        """
        ...

    def register_algorithm(self, algorithm_class: type[CSPAlgorithm]) -> None:
        """
        Registra novo algoritmo.

        Args:
            algorithm_class: Classe do algoritmo a ser registrada
        """
        ...

    def algorithm_exists(self, name: str) -> bool:
        """
        Verifica se algoritmo existe.

        Args:
            name: Nome do algoritmo

        Returns:
            bool: True se existe
        """
        ...

    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """
        Obtém metadados de algoritmo.

        Args:
            name: Nome do algoritmo

        Returns:
            Dict[str, Any]: Metadados do algoritmo
        """
        ...


@runtime_checkable
class ExportPort(Protocol):
    """Port para exportação de resultados."""

    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """
        Exporta resultados em formato específico.

        Args:
            results: Dados dos resultados
            format_type: Formato de exportação (csv, json, xlsx, etc.)
            destination: Destino da exportação

        Returns:
            str: Caminho do arquivo exportado
        """
        ...

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """
        Exporta resultados de batch.

        Args:
            batch_results: Lista de resultados
            format_type: Formato de exportação
            destination: Destino da exportação

        Returns:
            str: Caminho do arquivo exportado
        """
        ...

    def get_supported_formats(self) -> List[str]:
        """
        Lista formatos suportados.

        Returns:
            List[str]: Formatos disponíveis
        """
        ...

    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """
        Exporta resultados de otimização.

        Args:
            optimization_data: Dados da otimização
            destination: Destino da exportação

        Returns:
            str: Caminho do arquivo exportado
        """
        ...


@runtime_checkable
class ExecutionEngine(Protocol):
    """Port para engines de execução de tarefas do pipeline."""

    def run(
        self,
        task: Any,  # ExperimentTask | OptimizationTask | SensitivityTask
        dataset_obj: Dataset,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        monitor: Monitor | None = None,
        system_config: SystemConfig | None = None,
        check_control: Callable[[], str] | None = None,
        store: Any = None,
    ) -> Dict[str, Any]:
        """
        Executa uma tarefa específica.

        Args:
            task: Tarefa a ser executada
            dataset_obj: Dataset objeto
            alg: Parâmetros do algoritmo
            resources: Configuração de recursos
            monitor: Monitor para logs (padrão: NoOpMonitor)
            system_config: Configuração do sistema
            check_control: Função de controle pause/cancel (opcional)
            store: Store para persistência (opcional)

        Returns:
            Dict[str, Any]: Resultados da execução
        """
        ...


class AbstractAlgorithmRegistry(ABC):
    """Interface ABC para registry de algoritmos."""

    @abstractmethod
    def get_algorithm(self, name: str) -> type[CSPAlgorithm]:
        """Obtém classe de algoritmo por nome."""
        pass

    @abstractmethod
    def list_algorithms(self) -> List[str]:
        """Lista algoritmos disponíveis."""
        pass

    @abstractmethod
    def register_algorithm(self, algorithm_class: type[CSPAlgorithm]) -> None:
        """Registra novo algoritmo."""
        pass

    @abstractmethod
    def algorithm_exists(self, name: str) -> bool:
        """Verifica se algoritmo existe."""
        pass

    @abstractmethod
    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """Obtém metadados de algoritmo."""
        pass


class AbstractExportPort(ABC):
    """Interface ABC para exportação de resultados."""

    @abstractmethod
    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """Exporta resultados em formato específico."""
        pass

    @abstractmethod
    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """Exporta resultados de batch."""
        pass

    @abstractmethod
    def get_supported_formats(self) -> List[str]:
        """Lista formatos suportados."""
        pass

    @abstractmethod
    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """Exporta resultados de otimização."""
        pass


class AbstractStore(ABC):
    """Interface ABC para persistência de execuções."""

    @abstractmethod
    def record_execution_start(
        self,
        *,
        unit_id: str,
        task_id: str,
        dataset_id: str,
        algorithm: str,
        mode: str,
        repetition: int | None = None,
        trial: int | None = None,
        sample: int | None = None,
        params: Dict[str, Any] | None = None,
    ) -> None:
        """Registra início de execução."""
        pass

    @abstractmethod
    def record_execution_end(
        self,
        unit_id: str,
        status: str,
        result: Dict[str, Any] | None,
        objective: float | None,
    ) -> None:
        """Registra fim de execução."""
        pass

    @abstractmethod
    def get_completed_repetitions(
        self, task_id: str, dataset_id: str, algorithm: str
    ) -> List[Dict[str, Any]]:
        """Obtém repetições já completadas para resume."""
        pass


class AbstractExecutionEngine(ABC):
    """Interface ABC para engines de execução de tarefas."""

    @abstractmethod
    def run(
        self,
        task: Any,
        dataset_obj: Dataset,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        monitor: Monitor | None = None,
        system_config: SystemConfig | None = None,
        check_control: Callable[[], str] | None = None,
        store: AbstractStore | None = None,
    ) -> Dict[str, Any]:
        """Executa uma tarefa específica."""
        pass
