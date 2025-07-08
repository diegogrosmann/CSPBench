"""
ExecutionController - Camada de orquestração para execução de tarefas.

Responsável por registrar listeners (IConsole) e encaminhar callbacks
vindos do SchedulerExecutor, isolando a UI do controle de execução.
"""

import logging
from typing import Any, Dict, List, Optional, Protocol

from .execution_console import IConsole
from .executor import IExecutor, TaskHandle
from .task_result import TaskResult

logger = logging.getLogger(__name__)


class IExecutionListener(Protocol):
    """
    Interface para listeners de eventos de execução.

    Define os métodos que devem ser implementados por listeners
    que desejam receber eventos de execução de tarefas.
    """

    def on_task_start(self, task_id: str, algorithm_name: str, slot: Any) -> None:
        """
        Chamado quando uma tarefa inicia.

        Args:
            task_id: ID da tarefa
            algorithm_name: Nome do algoritmo
            slot: Slot da tarefa
        """
        ...

    def on_task_progress(self, task_id: str, message: str, slot: Any) -> None:
        """
        Chamado para atualizar progresso da tarefa.

        Args:
            task_id: ID da tarefa
            message: Mensagem de progresso
            slot: Slot da tarefa
        """
        ...

    def on_task_finish(self, task_id: str, result: TaskResult, slot: Any) -> None:
        """
        Chamado quando uma tarefa termina.

        Args:
            task_id: ID da tarefa
            result: Resultado da tarefa
            slot: Slot da tarefa
        """
        ...


class ExecutionController:
    """
    Controlador de execução que gerencia listeners e eventos.

    Recebe executor.submit(...) e apenas ouve eventos, sem decidir
    quantas tarefas iniciar. Encaminha eventos para listeners registrados.
    """

    def __init__(self, executor: IExecutor):
        """
        Inicializa o controlador de execução.

        Args:
            executor: Executor a ser usado para submissão de tarefas
        """
        self.executor = executor
        self.listeners: List[IConsole] = []
        self.active_tasks: Dict[str, Dict[str, Any]] = {}
        logger.info("ExecutionController inicializado")

    def add_listener(self, listener: IConsole) -> None:
        """
        Adiciona um listener para eventos de execução.

        Args:
            listener: Listener que implementa IConsole
        """
        self.listeners.append(listener)
        logger.debug("Listener adicionado: %s", type(listener).__name__)

    def remove_listener(self, listener: IConsole) -> None:
        """
        Remove um listener.

        Args:
            listener: Listener a ser removido
        """
        if listener in self.listeners:
            self.listeners.remove(listener)
            logger.debug("Listener removido: %s", type(listener).__name__)

    def submit_task(
        self, task_id: str, algorithm_instance: Any, meta: Dict[str, Any]
    ) -> TaskHandle:
        """
        Submete uma tarefa para execução.

        Args:
            task_id: ID único da tarefa
            algorithm_instance: Instância do algoritmo a ser executado
            meta: Metadados da tarefa

        Returns:
            TaskHandle: Handle da tarefa submetida
        """
        # Registrar tarefa
        self.active_tasks[task_id] = {"meta": meta, "handle": None, "status": "pending"}

        # Notificar listeners sobre início da tarefa
        self._notify_task_start(task_id, meta)

        # Submeter para o executor
        handle = self.executor.submit(algorithm_instance)
        self.active_tasks[task_id]["handle"] = handle
        self.active_tasks[task_id]["status"] = "running"

        logger.info("Tarefa %s submetida para execução", task_id)
        return handle

    def update_task_progress(
        self, task_id: str, progress: float, message: str = ""
    ) -> None:
        """
        Atualiza o progresso de uma tarefa.

        Args:
            task_id: ID da tarefa
            progress: Progresso (0.0 a 100.0)
            message: Mensagem opcional
        """
        if task_id not in self.active_tasks:
            logger.warning(
                f"Tentativa de atualizar progresso de tarefa inexistente: {task_id}"
            )
            return

        self._notify_task_progress(task_id, progress, message)

    def finish_task(self, task_id: str, result: TaskResult) -> None:
        """
        Finaliza uma tarefa.

        Args:
            task_id: ID da tarefa
            result: Resultado da execução
        """
        if task_id not in self.active_tasks:
            logger.warning("Tentativa de finalizar tarefa inexistente: %s", task_id)
            return

        self.active_tasks[task_id]["status"] = "completed"
        self._notify_task_finish(task_id, result)

        # Remover tarefa das ativas
        del self.active_tasks[task_id]
        logger.info("Tarefa %s finalizada", task_id)

    def get_task_status(self, task_id: str) -> Optional[str]:
        """
        Obtém o status de uma tarefa.

        Args:
            task_id: ID da tarefa

        Returns:
            Optional[str]: Status da tarefa ou None se não encontrada
        """
        if task_id in self.active_tasks:
            return self.active_tasks[task_id]["status"]
        return None

    def get_active_tasks(self) -> Dict[str, Dict[str, Any]]:
        """
        Retorna todas as tarefas ativas.

        Returns:
            Dict: Dicionário com tarefas ativas
        """
        return self.active_tasks.copy()

    def cleanup(self) -> None:
        """
        Limpa recursos e notifica listeners.
        """
        # Notificar listeners sobre limpeza
        for listener in self.listeners:
            try:
                listener.cleanup()
            except Exception as e:
                logger.error(
                    "Erro ao limpar listener %s: %s", type(listener).__name__, e
                )

        # Limpar estado
        self.active_tasks.clear()
        self.listeners.clear()
        logger.info("ExecutionController finalizado")

    def _notify_task_start(self, task_id: str, meta: Dict[str, Any]) -> None:
        """
        Notifica listeners sobre início de tarefa.

        Args:
            task_id: ID da tarefa
            meta: Metadados da tarefa
        """
        for listener in self.listeners:
            try:
                listener.on_task_start(task_id, meta)
            except Exception as e:
                logger.error(
                    f"Erro ao notificar listener {type(listener).__name__}: {e}"
                )

    def _notify_task_progress(
        self, task_id: str, progress: float, message: str = ""
    ) -> None:
        """
        Notifica listeners sobre progresso de tarefa.

        Args:
            task_id: ID da tarefa
            progress: Progresso (0.0 a 100.0)
            message: Mensagem opcional
        """
        for listener in self.listeners:
            try:
                listener.on_task_progress(task_id, progress, message)
            except Exception as e:
                logger.error(
                    f"Erro ao notificar listener {type(listener).__name__}: {e}"
                )

    def _notify_task_finish(self, task_id: str, result: TaskResult) -> None:
        """
        Notifica listeners sobre fim de tarefa.

        Args:
            task_id: ID da tarefa
            result: Resultado da execução
        """
        for listener in self.listeners:
            try:
                listener.on_task_finish(task_id, result)
            except Exception as e:
                logger.error(
                    f"Erro ao notificar listener {type(listener).__name__}: {e}"
                )


def create_execution_controller(executor: IExecutor) -> ExecutionController:
    """
    Cria um controlador de execução.

    Args:
        executor: Executor a ser usado

    Returns:
        ExecutionController: Controlador criado
    """
    return ExecutionController(executor)
