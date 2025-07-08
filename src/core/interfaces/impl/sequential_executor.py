"""
Executor sequencial para algoritmos.

Implementa a interface IExecutor para execução sequencial de algoritmos,
processando uma tarefa por vez.
"""

import logging
from typing import Union

from ..algorithm import IAlgorithm, Result
from ..executor import IExecutor, TaskHandle, TaskStatus

logger = logging.getLogger(__name__)


class SequentialExecutor(IExecutor):
    """
    Executor sequencial para algoritmos.

    Executa algoritmos um por vez, bloqueando até a conclusão de cada tarefa.
    """

    def __init__(self, timeout: int = 300, **kwargs):
        """
        Inicializa o executor sequencial.

        Args:
            timeout: Timeout por tarefa em segundos
            **kwargs: Argumentos adicionais (ignorados)
        """
        self.timeout = timeout
        self._shutdown = False

        # Importação tardia para evitar importação circular
        from ...scheduler.executor import SchedulerExecutor

        self._scheduler_executor = SchedulerExecutor(start_delay=0.0, timeout=timeout)

        logger.info("SequentialExecutor inicializado com timeout=%ss", timeout)

    def __enter__(self):
        """Context manager entry."""
        self._scheduler_executor.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self._scheduler_executor.__exit__(exc_type, exc_val, exc_tb)

    def submit(
        self, algorithm_instance: IAlgorithm, *, priority: int = 0
    ) -> TaskHandle:
        """
        Submete um algoritmo para execução sequencial.

        Args:
            algorithm_instance: Instância do algoritmo
            priority: Prioridade da tarefa (ignorada em execução sequencial)

        Returns:
            TaskHandle: Handle da tarefa
        """
        if self._shutdown:
            raise RuntimeError("Executor foi encerrado")

        # Usar o SchedulerExecutor diretamente
        handle = self._scheduler_executor.submit(algorithm_instance, priority=priority)

        logger.info("Tarefa %s submetida para execução sequencial", handle.task_id)
        return handle

    def poll(self, handle: TaskHandle) -> TaskStatus:
        """
        Verifica o status de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            TaskStatus: Status atual da tarefa
        """
        return self._scheduler_executor.poll(handle)

    def result(self, handle: TaskHandle) -> Union[Result, Exception]:
        """
        Obtém o resultado de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            Result ou Exception: Resultado da execução
        """
        return self._scheduler_executor.result(handle)

    def cancel(self, handle: TaskHandle) -> bool:
        """
        Cancela uma tarefa em execução.

        Args:
            handle: Handle da tarefa

        Returns:
            bool: True se a tarefa foi cancelada com sucesso
        """
        return self._scheduler_executor.cancel(handle)

    def shutdown(self, wait: bool = True) -> None:
        """
        Encerra o executor e limpa recursos.

        Args:
            wait: Se deve aguardar conclusão das tarefas pendentes
        """
        self._shutdown = True
        self._scheduler_executor.shutdown(wait=wait)
        logger.info("SequentialExecutor encerrado")

    def wait_for_all(
        self, handles: list[TaskHandle], timeout: int | None = None
    ) -> None:
        """
        Aguarda a conclusão de todas as tarefas.

        Args:
            handles: Lista de handles das tarefas
            timeout: Timeout em segundos
        """
        self._scheduler_executor.wait_for_all(handles, timeout=timeout)
