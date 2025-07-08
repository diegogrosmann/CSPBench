"""
Interface padronizada para executores de algoritmos.

Define o protocolo IExecutor que permite execução sequencial ou paralela
de algoritmos com a mesma API.
"""

from enum import Enum
from typing import Any, Optional, Protocol, Union
from uuid import uuid4

from .algorithm import IAlgorithm, Result


class TaskStatus(Enum):
    """Status de uma tarefa em execução."""

    QUEUED = "queued"  # Tarefa em fila, aguardando execução
    RUNNING = "running"  # Tarefa executando
    DONE = "done"  # Tarefa concluída
    ERROR = "error"  # Tarefa com erro
    TIMEOUT = "timeout"  # Tarefa com timeout


class TaskHandle:
    """
    Handle para uma tarefa em execução.

    Permite verificar status e obter resultados de tarefas assíncronas.
    """

    def __init__(self, task_id: str):
        self.task_id = task_id
        self._status = TaskStatus.RUNNING
        self._result: Optional[Union[Result, Exception]] = None
        self._error: Optional[Exception] = None

    def __str__(self) -> str:
        return f"TaskHandle({self.task_id}, {self._status.value})"

    def __repr__(self) -> str:
        return self.__str__()


class IExecutor(Protocol):
    """
    Interface padronizada para executores de algoritmos.

    Define uma API comum que funciona tanto para execução sequencial
    quanto paralela, permitindo troca transparente de implementações.
    """

    def submit(
        self, algorithm_instance: IAlgorithm, *, priority: int = 0
    ) -> TaskHandle:
        """
        Submete um algoritmo para execução.

        Args:
            algorithm_instance: Instância do algoritmo a ser executado
            priority: Prioridade da tarefa (0 = normal, valores menores = maior prioridade)

        Returns:
            TaskHandle: Handle para acompanhar o progresso da tarefa
        """
        ...

    def poll(self, handle: TaskHandle) -> TaskStatus:
        """
        Verifica o status atual de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            TaskStatus: Status atual da tarefa
        """
        ...

    def result(self, handle: TaskHandle) -> Union[Result, Exception]:
        """
        Obtém o resultado de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            Result ou Exception: Resultado da execução ou exceção em caso de erro

        Raises:
            RuntimeError: Se a tarefa ainda não foi finalizada
        """
        ...

    def cancel(self, handle: TaskHandle) -> bool:
        """
        Cancela uma tarefa em execução.

        Args:
            handle: Handle da tarefa

        Returns:
            bool: True se a tarefa foi cancelada com sucesso
        """
        ...

    def shutdown(self, wait: bool = True) -> None:
        """
        Encerra o executor e limpa recursos.

        Args:
            wait: Se True, aguarda todas as tarefas terminarem
        """
        ...


def create_task_handle() -> TaskHandle:
    """
    Cria um novo handle de tarefa com ID único.

    Returns:
        TaskHandle: Novo handle com ID único
    """
    task_id = f"task_{uuid4().hex[:8]}"
    return TaskHandle(task_id)
