"""Basic monitoring service compatible with the new interface."""

from typing import Any, Dict, Optional

from src.infrastructure.logging_config import get_logger
from src.presentation.monitoring.interfaces import TaskType
from src.presentation.monitoring.monitor_factory import MonitorFactory


class BasicMonitoringService:
    """Basic monitoring service for integration with execution system."""

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize monitoring service.

        Args:
            config: Batch/system configuration
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
        Start monitoring.

        Args:
            task_type: Task type (EXECUTION, OPTIMIZATION, SENSITIVITY)
            batch_name: Batch name
            batch_config: Batch configuration (optional)
        """
        try:
            # Create monitor based on configuration
            self.monitor = MonitorFactory.create_monitor(self.config)

            if not self.monitor:
                self.logger.info("Monitoring disabled in configuration")
                return

            # Start monitoring in interface
            self.monitor.start_task(task_type, batch_name, batch_config or {})
            self.is_active = True

            self.logger.info(f"Monitoring started for task: {task_type.value}")

        except Exception as e:
            self.logger.error(f"Error starting monitoring: {e}")
            self.monitor = None

    def show_error(self, error: str) -> None:
        """
        Display error in monitor.

        Args:
            error: Error message
        """
        if self.monitor:
            try:
                self.monitor.show_error(error)
            except Exception as e:
                self.logger.error(f"Error displaying error in monitor: {e}")

    def finish_monitoring(self, results: Optional[Dict[str, Any]] = None) -> None:
        """
        Finish monitoring.

        Args:
            results: Final results (optional)
        """
        if not self.is_active:
            return

        try:
            if self.monitor:
                success = bool(results is not None and results.get("results"))
                self.monitor.finish_task(success=success, final_results=results)

            self.is_active = False
            self.logger.info("Monitoring finished")

        except Exception as e:
            self.logger.error(f"Error finishing monitoring: {e}")

    def close(self) -> None:
        """Close monitoring service."""
        if self.monitor:
            try:
                self.monitor.close()
            except Exception as e:
                self.logger.error(f"Error closing monitor: {e}")

        self.is_active = False
        self.monitor = None

    def update_item(
        self, item_id: str, progress: float, message: str = "", context=None
    ) -> None:
        """
        Update an individual item.

        Args:
            item_id: Unique item ID
            progress: Progress (0.0 to 100.0)
            message: Status message
            context: Hierarchical context (optional)
        """
        if self.monitor:
            try:
                self.monitor.update_item(item_id, progress, message, context)
            except Exception as e:
                self.logger.error(f"Error updating item {item_id}: {e}")

    def start_item(
        self,
        item_id: str,
        item_type: str = "repetition",
        context=None,
        metadata=None,
    ) -> None:
        """
        Start an individual item.

        Args:
            item_id: Unique item ID
            item_type: Item type (repetition, trial, sample, etc.)
            context: Hierarchical context (optional)
            metadata: Optional metadata
        """
        if self.monitor:
            try:
                self.monitor.start_item(item_id, item_type, context, metadata)
            except Exception as e:
                self.logger.error(f"Error starting item {item_id}: {e}")

    def update_hierarchy(
        self,
        level,
        level_id: str,
        progress: float,
        message: str = "",
        data=None,
    ) -> None:
        """
        Update hierarchical progress.

        Args:
            level: Hierarchical level (ExecutionLevel)
            level_id: Level ID
            progress: Progress (0.0 to 100.0)
            message: Status message
            data: Additional level-specific data
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
