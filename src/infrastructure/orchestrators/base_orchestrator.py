"""
Base para Orquestradores

Fornece funcionalidade comum para todos os orquestradores,
incluindo integração padronizada com sistema de monitoramento.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

from src.infrastructure.logging_config import get_logger


class BaseOrchestrator(ABC):
    """Base abstrata para todos os orquestradores do sistema."""

    def __init__(self, monitoring_service=None):
        """
        Inicializa orquestrador base.

        Args:
            monitoring_service: Serviço de monitoramento opcional
        """
        self.monitoring_service = monitoring_service
        self._logger = get_logger(__name__)

    @abstractmethod
    def execute(self, **kwargs) -> Dict[str, Any]:
        """
        Executa a orquestração específica.

        Returns:
            Dict[str, Any]: Resultado da orquestração
        """
        pass

    def _report_progress(self, progress: float, message: str) -> None:
        """
        Reporta progresso ao sistema de monitoramento.

        Args:
            progress: Progresso em percentual (0-100)
            message: Mensagem descritiva do progresso
        """
        if self.monitoring_service:
            try:
                self.monitoring_service.report_progress(progress, message)
            except Exception as e:
                self._logger.warning(f"Erro ao reportar progresso: {e}")

    def _update_task_data(self, task_type: str, **data) -> None:
        """
        Atualiza dados específicos da tarefa no monitoramento.

        Args:
            task_type: Tipo da tarefa (execution, optimization, sensitivity)
            **data: Dados específicos da tarefa
        """
        if not self.monitoring_service:
            return

        try:
            if task_type == "execution":
                self.monitoring_service.update_execution_data(**data)
            elif task_type == "optimization":
                self.monitoring_service.update_optimization_data(**data)
            elif task_type == "sensitivity":
                self.monitoring_service.update_sensitivity_data(**data)
            else:
                self._logger.warning(f"Tipo de tarefa desconhecido: {task_type}")
        except Exception as e:
            self._logger.warning(f"Erro ao atualizar dados da tarefa: {e}")

    def _log_error(self, error: Exception, context: str = "") -> None:
        """
        Registra erro com contexto apropriado.

        Args:
            error: Exceção ocorrida
            context: Contexto adicional do erro
        """
        error_msg = f"Erro na orquestração"
        if context:
            error_msg += f" ({context})"
        error_msg += f": {error}"

        self._logger.error(error_msg)

    def _log_info(self, message: str) -> None:
        """
        Registra informação no log.

        Args:
            message: Mensagem a registrar
        """
        self._logger.info(message)

    def _log_debug(self, message: str) -> None:
        """
        Registra informação de debug no log.

        Args:
            message: Mensagem a registrar
        """
        self._logger.debug(message)
