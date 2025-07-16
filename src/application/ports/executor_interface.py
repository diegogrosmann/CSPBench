"""
Interface Padronizada para Executores

Define o contrato comum que todos os executores devem implementar,
garantindo consistência na execução de diferentes tipos de tarefas.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

from src.domain import Dataset


class ExecutorInterface(ABC):
    """Interface padronizada para executores de algoritmos CSP."""

    @abstractmethod
    def execute_batch(
        self,
        batch_config: Dict[str, Any],
        monitoring_service=None,
    ) -> List[Dict[str, Any]]:
        """
        Executa batch de algoritmos (incluindo execuções únicas).

        Args:
            batch_config: Configuração do batch
            monitoring_service: Serviço de monitoramento opcional

        Returns:
            List[Dict[str, Any]]: Lista de resultados da execução
        """
        pass

    @abstractmethod
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
    ) -> Dict[str, Any]:
        """
        Executa otimização de hiperparâmetros.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para otimização
            optimization_config: Configuração da otimização
            monitoring_service: Serviço de monitoramento opcional
            config_index: Índice da configuração atual
            total_configs: Total de configurações
            dataset_index: Índice do dataset atual
            total_datasets: Total de datasets

        Returns:
            Dict[str, Any]: Resultado da otimização contendo:
                - best_params: melhores parâmetros encontrados
                - best_value: melhor valor alcançado
                - n_trials: número de trials executados
                - optimization_history: histórico da otimização
        """
        pass

    @abstractmethod
    def execute_sensitivity_analysis(
        self,
        algorithm_name: str,
        dataset: Dataset,
        sensitivity_config: Dict[str, Any],
        monitoring_service=None,
    ) -> Dict[str, Any]:
        """
        Executa análise de sensibilidade.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para análise
            sensitivity_config: Configuração da análise
            monitoring_service: Serviço de monitoramento opcional

        Returns:
            Dict[str, Any]: Resultado da análise contendo:
                - sensitivity_indices: índices de sensibilidade
                - parameter_rankings: ranking dos parâmetros
                - analysis_method: método utilizado
                - samples_executed: número de amostras executadas
        """
        pass

    @abstractmethod
    def get_execution_status(self, execution_id: str) -> str:
        """
        Obtém status de uma execução específica.

        Args:
            execution_id: ID da execução

        Returns:
            str: Status da execução (running, completed, failed, cancelled)
        """
        pass

    @abstractmethod
    def cancel_execution(self, execution_id: str) -> bool:
        """
        Cancela uma execução em andamento.

        Args:
            execution_id: ID da execução a cancelar

        Returns:
            bool: True se cancelamento foi bem-sucedido
        """
        pass

    def set_batch_config(self, batch_config: Dict[str, Any]) -> None:
        """
        Define configuração do batch atual (método opcional).

        Args:
            batch_config: Configuração do batch
        """
        pass
