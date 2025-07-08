"""
Interface IConsole para isolamento da camada de UI.

Define o protocolo padrão para consoles que recebem eventos de execução
de tarefas, permitindo diferentes implementações (SimpleConsole, CursesConsole).
"""

from typing import Any, Dict, Protocol

from ..interfaces.task_result import TaskResult


class IConsole(Protocol):
    """
    Interface padronizada para consoles de execução.

    Define métodos para receber eventos de tarefas durante a execução,
    permitindo diferentes implementações de UI sem acoplamento.
    """

    def on_task_start(self, task_id: str, meta: Dict[str, Any]) -> None:
        """
        Chamado quando uma tarefa é iniciada.

        Args:
            task_id: ID único da tarefa
            meta: Metadados da tarefa (algoritmo, parâmetros, etc.)
        """
        ...

    def on_task_progress(self, task_id: str, pct: float, msg: str = "") -> None:
        """
        Chamado para atualizar o progresso de uma tarefa.

        Args:
            task_id: ID único da tarefa
            pct: Porcentagem de progresso (0.0 a 100.0)
            msg: Mensagem opcional de progresso
        """
        ...

    def on_task_finish(self, task_id: str, result: TaskResult) -> None:
        """
        Chamado quando uma tarefa é finalizada.

        Args:
            task_id: ID único da tarefa
            result: Resultado da execução da tarefa
        """
        ...

    def cleanup(self) -> None:
        """
        Limpa recursos do console.

        Chamado no final da execução para garantir limpeza adequada.
        """
        ...
