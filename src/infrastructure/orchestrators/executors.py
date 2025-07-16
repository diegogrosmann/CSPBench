"""
Executor Principal - Roteador de Tarefas

Implementa interface padronizada delegando execução para orquestradores especializados.
"""

import time
from typing import Any, Dict, List

from src.application.ports import ExecutorInterface
from src.domain import Dataset
from src.domain.errors import AlgorithmExecutionError
from src.infrastructure.logging_config import get_logger


class Executor(ExecutorInterface):
    """Executor principal que delega tarefas para orquestradores especializados."""

    def __init__(self):
        """Inicializa executor como roteador puro."""
        self._logger = get_logger(__name__)
        self._current_batch_config = None

    def set_batch_config(self, batch_config: Dict[str, Any]) -> None:
        """Define configuração do batch atual (delegado para orquestradores)."""
        self._current_batch_config = batch_config
        self._logger.debug(f"Configuração de batch definida: {type(batch_config)}")

    def execute_batch(
        self, batch_config: Dict[str, Any], monitoring_service=None
    ) -> List[Dict[str, Any]]:
        """
        Executa batch de algoritmos delegando para ExecutionOrchestrator.

        Args:
            batch_config: Configuração do batch
            monitoring_service: Serviço de monitoramento opcional

        Returns:
            List[Dict[str, Any]]: Lista de resultados da execução
        """
        try:
            from src.infrastructure import (
                DomainAlgorithmRegistry,
                FileDatasetRepository,
            )
            from src.infrastructure.orchestrators.execution_orchestrator import (
                ExecutionOrchestrator,
            )

            # Configurar dependências
            algorithm_registry = DomainAlgorithmRegistry()
            dataset_repository = FileDatasetRepository("./datasets")

            # Criar orquestrador de execução
            orchestrator = ExecutionOrchestrator(
                algorithm_registry=algorithm_registry,
                dataset_repository=dataset_repository,
                monitoring_service=monitoring_service,
            )

            # Delegar execução
            return orchestrator.execute_batch(batch_config, monitoring_service)

        except Exception as e:
            self._logger.error(f"Erro na execução do batch: {e}")
            raise AlgorithmExecutionError(f"Erro na execução do batch: {e}") from e

    def execute_optimization(
        self,
        algorithm_name: str,
        dataset: Dataset,
        optimization_config: Dict[str, Any],
        monitoring_service=None,
        config_index: int = 1,
        total_configs: int = 1,
        dataset_index: int = 1,
        total_datasets: int = 1,
        dataset_name: str | None = None,
    ) -> Dict[str, Any]:
        """
        Executa otimização de hiperparâmetros delegando para OptimizationOrchestrator.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para otimização
            optimization_config: Configuração da otimização
            monitoring_service: Serviço de monitoramento opcional
            config_index: Índice da configuração atual
            total_configs: Total de configurações
            dataset_index: Índice do dataset atual
            total_datasets: Total de datasets
            dataset_name: Nome original do dataset (opcional)

        Returns:
            Dict[str, Any]: Resultado da otimização
        """
        try:
            from src.infrastructure import (
                DomainAlgorithmRegistry,
                FileDatasetRepository,
            )
            from src.infrastructure.orchestrators.optimization_orchestrator import (
                OptimizationOrchestrator,
            )

            # Configurar dependências
            algorithm_registry = DomainAlgorithmRegistry()
            dataset_repository = FileDatasetRepository("./datasets")

            # Usar nome original do dataset se fornecido, senão usar temporário
            if dataset_name:
                dataset_id = dataset_name
            else:
                dataset_id = f"temp_optimization_{int(time.time())}"

            # Salvar o dataset temporariamente no repositório para o orquestrador
            dataset_repository.save(dataset, dataset_id)

            # Criar configuração completa para o orquestrador
            full_config = {
                "algorithm": algorithm_name,
                "dataset": dataset_id,
                "optimization": optimization_config,
                "export": optimization_config.get("export", {"enabled": True}),
                "monitoring": optimization_config.get("monitoring", {}),
                "resources": optimization_config.get("resources", {}),
                "plots": optimization_config.get("plots", {}),
            }

            # Criar orquestrador
            orchestrator = OptimizationOrchestrator(
                algorithm_registry=algorithm_registry,
                dataset_repository=dataset_repository,
                config=full_config,
                monitoring_service=monitoring_service,
                config_index=config_index,
                total_configs=total_configs,
                dataset_index=dataset_index,
                total_datasets=total_datasets,
            )

            # Executar otimização
            results = orchestrator.run_optimization()

            # Limpar dataset temporário apenas se era temporário
            if not dataset_name:
                try:
                    dataset_repository.delete(dataset_id)
                except:
                    pass

            return results

        except Exception as e:
            self._logger.error(f"Erro na otimização: {e}")
            raise AlgorithmExecutionError(f"Erro na otimização: {e}") from e

    def execute_sensitivity_analysis(
        self,
        algorithm_name: str,
        dataset: Dataset,
        sensitivity_config: Dict[str, Any],
        monitoring_service=None,
    ) -> Dict[str, Any]:
        """
        Executa análise de sensibilidade delegando para SensitivityOrchestrator.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para análise
            sensitivity_config: Configuração da análise
            monitoring_service: Serviço de monitoramento opcional

        Returns:
            Dict[str, Any]: Resultado da análise de sensibilidade
        """
        try:
            from src.infrastructure import DomainAlgorithmRegistry
            from src.infrastructure.orchestrators.sensitivity_orchestrator import (
                SensitivityOrchestrator,
            )

            # Configurar registry
            algorithm_registry = DomainAlgorithmRegistry()

            # Criar orquestrador de sensibilidade
            orchestrator = SensitivityOrchestrator(
                algorithm_registry, self, monitoring_service=monitoring_service
            )

            # Executar análise
            results = orchestrator.execute_sensitivity_analysis(
                algorithm_name, dataset, sensitivity_config
            )

            return results

        except Exception as e:
            self._logger.error(f"Erro na análise de sensibilidade: {e}")
            raise AlgorithmExecutionError(
                f"Erro na análise de sensibilidade: {e}"
            ) from e

    def get_execution_status(self, execution_id: str) -> str:
        """
        Obtém status de uma execução específica.

        Args:
            execution_id: ID da execução

        Returns:
            str: Status da execução (running, completed, failed, cancelled)
        """
        # Por enquanto, retornamos um status genérico
        # TODO: Implementar sistema de tracking de execuções
        return "completed"

    def cancel_execution(self, execution_id: str) -> bool:
        """
        Cancela uma execução em andamento.

        Args:
            execution_id: ID da execução a cancelar

        Returns:
            bool: True se cancelamento foi bem-sucedido
        """
        # Por enquanto, retornamos False
        # TODO: Implementar sistema de cancelamento
        return False
