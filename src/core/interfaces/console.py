"""
Interface padronizada para consoles.

Define o protocolo IConsole que permite abstração da camada de UI,
com fallback automático entre CursesConsole e SimpleConsole.
"""

from typing import Optional, Protocol

from src.utils.curses_console import BaseConsoleManager, create_console_manager


class TaskSlot:
    """
    Representa um slot de tarefa no console.

    Permite rastrear e atualizar o progresso de tarefas individuais
    em uma interface visual organizada.
    """

    def __init__(self, slot_id: int, name: str):
        self.slot_id = slot_id
        self.name = name
        self.active = False
        self.current_message = ""

    def __str__(self) -> str:
        return f"TaskSlot({self.slot_id}, {self.name}, active={self.active})"

    def __repr__(self) -> str:
        return self.__str__()


class IConsole(Protocol):
    """
    Interface padronizada para consoles.

    Define uma camada fina sobre CursesConsole com fallback automático
    para SimpleConsole, permitindo UI consistente em diferentes ambientes.
    """

    def show_task_start(self, name: str, slot: TaskSlot) -> None:
        """
        Mostra o início de uma tarefa em um slot específico.

        Args:
            name: Nome da tarefa
            slot: Slot onde a tarefa será exibida
        """
        ...

    def update_progress(self, slot: TaskSlot, message: str) -> None:
        """
        Atualiza o progresso de uma tarefa em um slot.

        Args:
            slot: Slot da tarefa
            message: Mensagem de progresso
        """
        ...

    def show_task_end(
        self, slot: TaskSlot, success: bool, message: str, error: Optional[str] = None
    ) -> None:
        """
        Mostra o fim de uma tarefa.

        Args:
            slot: Slot da tarefa
            success: True se a tarefa foi bem-sucedida
            message: Mensagem final
            error: Mensagem de erro (opcional)
        """
        ...

    def print(self, *args, **kwargs) -> None:
        """
        Imprime uma mensagem no console.

        Args:
            *args: Argumentos para impressão
            **kwargs: Argumentos nomeados para impressão
        """
        ...

    def clear(self) -> None:
        """Limpa a tela do console."""
        ...

    def cleanup(self) -> None:
        """Limpa recursos do console."""
        ...


class ConsoleAdapter:
    """
    Adaptador que implementa IConsole sobre BaseConsoleManager.

    Fornece a interface padronizada com fallback automático e
    gerenciamento de slots de tarefas.
    """

    def __init__(self, force_simple: bool = False):
        self._console_manager = create_console_manager(force_simple)
        self._slots: dict[int, TaskSlot] = {}
        self._next_slot_id = 1

    def allocate_slot(self, name: str) -> TaskSlot:
        """
        Aloca um novo slot para uma tarefa.

        Args:
            name: Nome da tarefa

        Returns:
            TaskSlot: Novo slot alocado
        """
        slot_id = self._next_slot_id
        self._next_slot_id += 1

        slot = TaskSlot(slot_id, name)
        self._slots[slot_id] = slot
        return slot

    def release_slot(self, slot: TaskSlot) -> None:
        """
        Libera um slot de tarefa.

        Args:
            slot: Slot a ser liberado
        """
        if slot.slot_id in self._slots:
            del self._slots[slot.slot_id]

    def show_task_start(self, name: str, slot: TaskSlot) -> None:
        """Mostra o início de uma tarefa."""
        slot.active = True
        slot.current_message = f"Iniciando {name}..."
        self._console_manager.print(f"🔄 [{slot.slot_id}] {slot.current_message}")

    def update_progress(self, slot: TaskSlot, message: str) -> None:
        """Atualiza o progresso de uma tarefa."""
        if slot.active:
            slot.current_message = message
            self._console_manager.print(f"⏳ [{slot.slot_id}] {message}")

    def show_task_end(
        self, slot: TaskSlot, success: bool, message: str, error: Optional[str] = None
    ) -> None:
        """Mostra o fim de uma tarefa."""
        slot.active = False
        slot.current_message = message

        if success:
            icon = "✅"
            display_message = message
        elif error and "timeout" in error.lower():
            icon = "⏰"
            display_message = f"{message} (Timeout)"
        else:
            icon = "❌"
            display_message = f"{message} - {error}" if error else message

        self._console_manager.print(f"{icon} [{slot.slot_id}] {display_message}")

    def print(self, *args, **kwargs) -> None:
        """Imprime uma mensagem no console."""
        self._console_manager.print(*args, **kwargs)

    def clear(self) -> None:
        """Limpa a tela do console."""
        self._console_manager.clear()

    def cleanup(self) -> None:
        """Limpa recursos do console."""
        self._console_manager.cleanup()

    def show_status(self, message: str) -> None:
        """
        Mostra uma mensagem de status geral.

        Args:
            message: Mensagem de status
        """
        if hasattr(self._console_manager, "show_status"):
            self._console_manager.show_status(message)
        else:
            self.print(f"Status: {message}")


def create_console(force_simple: bool = False) -> ConsoleAdapter:
    """
    Cria uma instância de console com fallback automático.

    Args:
        force_simple: Se True, força uso do console simples

    Returns:
        ConsoleAdapter: Instância do console adaptado
    """
    return ConsoleAdapter(force_simple)
