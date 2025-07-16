"""Serviço básico de monitoramento compatível com a nova interface."""

from typing import Any, Dict, Optional

from src.infrastructure.logging_config import get_logger
from src.presentation.monitoring.interfaces import TaskType
from src.presentation.monitoring.monitor_factory import MonitorFactory


class BasicMonitoringService:
    """Serviço básico de monitoramento para integração com o sistema de execução."""

    def __init__(self, config: Dict[str, Any]):
        """
        Inicializa o serviço de monitoramento.

        Args:
            config: Configuração do batch/sistema
        """
        self.config = config
        self.logger = get_logger(__name__)
        self.monitor = None
        self.is_active = False

    def start_monitoring(
        self,
        task_type: TaskType,
        batch_name: str,
        batch_config: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Inicia o monitoramento.

        Args:
            task_type: Tipo da tarefa (EXECUTION, OPTIMIZATION, SENSITIVITY)
            batch_name: Nome do batch
            batch_config: Configuração do batch (opcional)
        """
        try:
            # Cria monitor baseado na configuração
            self.monitor = MonitorFactory.create_monitor(self.config)

            if not self.monitor:
                self.logger.info("Monitoramento desabilitado na configuração")
                return

            # Inicia o monitoramento na interface
            self.monitor.start_task(task_type, batch_name, batch_config or {})
            self.is_active = True

            self.logger.info(f"Monitoramento iniciado para tarefa: {task_type.value}")

        except Exception as e:
            self.logger.error(f"Erro ao iniciar monitoramento: {e}")
            self.monitor = None

    def show_error(self, error: str) -> None:
        """
        Exibe erro no monitor.

        Args:
            error: Mensagem de erro
        """
        if self.monitor:
            try:
                self.monitor.show_error(error)
            except Exception as e:
                self.logger.error(f"Erro ao exibir erro no monitor: {e}")

    def finish_monitoring(self, results: Optional[Dict[str, Any]] = None) -> None:
        """
        Finaliza o monitoramento.

        Args:
            results: Resultados finais (opcional)
        """
        if not self.is_active:
            return

        try:
            if self.monitor:
                success = bool(results is not None and results.get("results"))
                self.monitor.finish_task(success=success, final_results=results)

            self.is_active = False
            self.logger.info("Monitoramento finalizado")

        except Exception as e:
            self.logger.error(f"Erro ao finalizar monitoramento: {e}")

    def close(self) -> None:
        """Fecha o serviço de monitoramento."""
        if self.monitor:
            try:
                self.monitor.close()
            except Exception as e:
                self.logger.error(f"Erro ao fechar monitor: {e}")

        self.is_active = False
        self.monitor = None

    def update_item(
        self, item_id: str, progress: float, message: str = "", context=None
    ) -> None:
        """
        Atualiza um item individual.

        Args:
            item_id: ID único do item
            progress: Progresso (0.0 a 100.0)
            message: Mensagem de status
            context: Contexto hierárquico (opcional)
        """
        if self.monitor:
            try:
                self.monitor.update_item(item_id, progress, message, context)
            except Exception as e:
                self.logger.error(f"Erro ao atualizar item {item_id}: {e}")

    def start_item(
        self,
        item_id: str,
        item_type: str = "repetition",
        context=None,
        metadata=None,
    ) -> None:
        """
        Inicia um item individual.

        Args:
            item_id: ID único do item
            item_type: Tipo do item (repetition, trial, sample, etc.)
            context: Contexto hierárquico (opcional)
            metadata: Metadados opcionais
        """
        if self.monitor:
            try:
                self.monitor.start_item(item_id, item_type, context, metadata)
            except Exception as e:
                self.logger.error(f"Erro ao iniciar item {item_id}: {e}")

    def update_hierarchy(
        self,
        level,
        level_id: str,
        progress: float,
        message: str = "",
        data=None,
    ) -> None:
        """
        Atualiza progresso hierárquico.

        Args:
            level: Nível hierárquico (ExecutionLevel)
            level_id: ID do nível
            progress: Progresso (0.0 a 100.0)
            message: Mensagem de status
            data: Dados adicionais específicos do nível
        """
        if self.monitor:
            try:
                self.monitor.update_hierarchy(level, level_id, progress, message, data)
            except Exception as e:
                self.logger.error(f"Erro ao atualizar hierarquia {level}: {e}")

    def finish_item(
        self,
        item_id: str,
        success: bool = True,
        result=None,
        error: Optional[str] = None,
    ) -> None:
        """
        Finaliza um item individual.

        Args:
            item_id: ID único do item
            success: Se foi executado com sucesso
            result: Resultado da execução
            error: Mensagem de erro se falhou
        """
        if self.monitor:
            try:
                self.monitor.finish_item(item_id, success, result, error)
            except Exception as e:
                self.logger.error(f"Erro ao finalizar item {item_id}: {e}")

    def algorithm_callback(
        self,
        algorithm_name: str,
        progress: float,
        message: str = "",
        item_id: Optional[str] = None,
    ) -> None:
        """
        Callback direto do algoritmo durante execução.

        Args:
            algorithm_name: Nome do algoritmo
            progress: Progresso (0.0 a 100.0)
            message: Mensagem de status
            item_id: ID único do item (opcional)
        """
        if self.monitor:
            try:
                self.monitor.algorithm_callback(
                    algorithm_name, progress, message, item_id
                )
            except Exception as e:
                self.logger.error(
                    f"Erro no callback do algoritmo {algorithm_name}: {e}"
                )
