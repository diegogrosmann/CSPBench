"""
Listener que conecta ExecutionController com IConsole.

Este módulo fornece um listener que traduz eventos do ExecutionController
para chamadas na interface IConsole, permitindo atualização visual das tarefas.
"""

import logging

from src.core.interfaces.console import IConsole, TaskSlot
from src.core.interfaces.execution_controller import IExecutionListener
from src.core.interfaces.task_result import TaskResult

logger = logging.getLogger(__name__)


class ConsoleExecutionListener(IExecutionListener):
    """
    Listener que conecta ExecutionController com IConsole.

    Este listener recebe eventos do ExecutionController e os traduz
    para chamadas apropriadas na interface IConsole.
    """

    def __init__(self, console: IConsole):
        """
        Inicializa o listener.

        Args:
            console: Interface de console para output
        """
        self.console = console
        self.task_slots = {}
        logger.info("ConsoleExecutionListener inicializado")

    def on_task_start(self, task_id: str, algorithm_name: str, slot: TaskSlot) -> None:
        """
        Chamado quando uma tarefa inicia.

        Args:
            task_id: ID da tarefa
            algorithm_name: Nome do algoritmo
            slot: Slot da tarefa
        """
        try:
            self.task_slots[task_id] = slot
            self.console.show_task_start(algorithm_name, slot)
            logger.debug(
                "Tarefa %s (%s) iniciada no slot %s",
                task_id,
                algorithm_name,
                slot.slot_id,
            )
        except Exception as e:
            logger.error("Erro ao mostrar início da tarefa %s: %s", task_id, e)

    def on_task_progress(self, task_id: str, message: str, slot: TaskSlot) -> None:
        """
        Chamado para atualizar progresso da tarefa.

        Args:
            task_id: ID da tarefa
            message: Mensagem de progresso
            slot: Slot da tarefa
        """
        try:
            if task_id in self.task_slots:
                self.console.update_progress(slot, message)
                logger.debug("Progresso da tarefa %s: %s", task_id, message)
        except Exception as e:
            logger.error("Erro ao atualizar progresso da tarefa %s: %s", task_id, e)

    def on_task_finish(self, task_id: str, result: TaskResult, slot: TaskSlot) -> None:
        """
        Chamado quando uma tarefa termina.

        Args:
            task_id: ID da tarefa
            result: Resultado da tarefa
            slot: Slot da tarefa
        """
        try:
            if task_id in self.task_slots:
                success = result.success

                if success:
                    distance = result.distance or 0
                    time_taken = result.time or 0.0
                    message = (
                        f"Concluído - Distância: {distance}, Tempo: {time_taken:.2f}s"
                    )
                else:
                    message = "Falhou"

                self.console.show_task_end(slot, success, message, result.error)

                logger.info("Tarefa %s concluída - Sucesso: %s", task_id, success)

                # Limpar slot
                del self.task_slots[task_id]

        except Exception as e:
            logger.error("Erro ao mostrar conclusão da tarefa %s: %s", task_id, e)


def create_console_listener(console: IConsole) -> ConsoleExecutionListener:
    """
    Cria um listener de console.

    Args:
        console: Interface de console

    Returns:
        Listener configurado
    """
    return ConsoleExecutionListener(console)
