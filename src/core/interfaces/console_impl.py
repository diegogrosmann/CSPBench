"""
Implementações básicas da interface IConsole.

Este módulo fornece implementações concretas de IConsole para diferentes
ambientes de execução, incluindo console simples e integração com curses.
"""

from typing import Optional

from src.core.interfaces.console import IConsole, TaskSlot


class SimpleConsole:
    """
    Implementação simples de IConsole que usa prints padrão.

    Esta implementação fornece uma saída básica para ambientes
    onde curses não está disponível ou não é desejado.
    """

    def __init__(self):
        """Inicializa o console simples."""
        self.active_tasks = {}

    def show_task_start(self, name: str, slot: TaskSlot) -> None:
        """
        Mostra o início de uma tarefa.

        Args:
            name: Nome da tarefa
            slot: Slot da tarefa
        """
        print(f"🚀 Iniciando: {name} [Slot {slot.slot_id}]")
        self.active_tasks[slot.slot_id] = {"name": name, "slot": slot, "active": True}
        slot.active = True
        slot.current_message = "Iniciando..."

    def update_progress(self, slot: TaskSlot, message: str) -> None:
        """
        Atualiza o progresso de uma tarefa.

        Args:
            slot: Slot da tarefa
            message: Mensagem de progresso
        """
        if slot.slot_id in self.active_tasks:
            slot.current_message = message
            print(f"⏳ {slot.name}: {message}")

    def show_task_end(
        self, slot: TaskSlot, success: bool, message: str, error: Optional[str] = None
    ) -> None:
        """
        Mostra o fim de uma tarefa.

        Args:
            slot: Slot da tarefa
            success: True se a tarefa foi bem-sucedida
            message: Mensagem final
            error: Mensagem de erro opcional
        """
        status_icon = "✅" if success else "❌"
        print(f"{status_icon} {slot.name}: {message}")

        if error:
            print(f"   Erro: {error}")

        slot.active = False
        slot.current_message = message

        if slot.slot_id in self.active_tasks:
            del self.active_tasks[slot.slot_id]

    def print(self, *args, **kwargs) -> None:
        """
        Imprime uma mensagem no console.

        Args:
            *args: Argumentos para print
            **kwargs: Argumentos nomeados para print
        """
        print(*args, **kwargs)

    def clear(self) -> None:
        """Limpa a tela do console."""
        # Em console simples, apenas adiciona algumas linhas
        print("\n" * 3)

    def cleanup(self) -> None:
        """Limpa recursos do console."""
        # Console simples não precisa de cleanup
        pass


class CursesConsoleAdapter:
    """
    Adaptador para integrar CursesInterface com IConsole.

    Esta classe adapta a interface curses existente para funcionar
    com o novo sistema de IConsole, mantendo compatibilidade.
    """

    def __init__(self, curses_interface=None):
        """
        Inicializa o adaptador.

        Args:
            curses_interface: Interface curses opcional
        """
        self.curses_interface = curses_interface
        self.active_slots = {}

    def show_task_start(self, name: str, slot: TaskSlot) -> None:
        """
        Mostra o início de uma tarefa no curses.

        Args:
            name: Nome da tarefa
            slot: Slot da tarefa
        """
        self.active_slots[slot.slot_id] = slot
        slot.active = True
        slot.current_message = "Iniciando..."

        if self.curses_interface:
            try:
                # Assumir que a interface curses tem um método para mostrar
                # tarefas
                if hasattr(self.curses_interface, "show_task_start"):
                    self.curses_interface.show_task_start(name, slot)
                elif hasattr(self.curses_interface, "update_task_status"):
                    self.curses_interface.update_task_status(
                        slot.slot_id, name, "Iniciando..."
                    )
            except Exception as e:
                # Fallback para console simples
                print(f"🚀 Iniciando: {name} [Slot {slot.slot_id}]")

    def update_progress(self, slot: TaskSlot, message: str) -> None:
        """
        Atualiza o progresso de uma tarefa no curses.

        Args:
            slot: Slot da tarefa
            message: Mensagem de progresso
        """
        slot.current_message = message

        if self.curses_interface:
            try:
                if hasattr(self.curses_interface, "update_progress"):
                    self.curses_interface.update_progress(slot, message)
                elif hasattr(self.curses_interface, "update_task_status"):
                    self.curses_interface.update_task_status(
                        slot.slot_id, slot.name, message
                    )
            except Exception as e:
                # Fallback para console simples
                print(f"⏳ {slot.name}: {message}")

    def show_task_end(
        self, slot: TaskSlot, success: bool, message: str, error: Optional[str] = None
    ) -> None:
        """
        Mostra o fim de uma tarefa no curses.

        Args:
            slot: Slot da tarefa
            success: True se a tarefa foi bem-sucedida
            message: Mensagem final
            error: Mensagem de erro opcional
        """
        slot.active = False
        slot.current_message = message

        if self.curses_interface:
            try:
                if hasattr(self.curses_interface, "show_task_end"):
                    self.curses_interface.show_task_end(slot, success, message, error)
                elif hasattr(self.curses_interface, "update_task_status"):
                    status_msg = message if success else f"ERRO: {error or message}"
                    self.curses_interface.update_task_status(
                        slot.slot_id, slot.name, status_msg
                    )
            except Exception as e:
                # Fallback para console simples
                status_icon = "✅" if success else "❌"
                print(f"{status_icon} {slot.name}: {message}")
                if error:
                    print(f"   Erro: {error}")

        if slot.slot_id in self.active_slots:
            del self.active_slots[slot.slot_id]

    def print(self, *args, **kwargs) -> None:
        """
        Imprime uma mensagem no console.

        Args:
            *args: Argumentos para print
            **kwargs: Argumentos nomeados para print
        """
        if self.curses_interface:
            try:
                if hasattr(self.curses_interface, "print"):
                    self.curses_interface.print(*args, **kwargs)
                elif hasattr(self.curses_interface, "add_message"):
                    message = " ".join(str(arg) for arg in args)
                    self.curses_interface.add_message(message)
                else:
                    # Fallback
                    print(*args, **kwargs)
            except Exception as e:
                print(*args, **kwargs)
        else:
            print(*args, **kwargs)

    def clear(self) -> None:
        """Limpa a tela do console."""
        if self.curses_interface:
            try:
                if hasattr(self.curses_interface, "clear"):
                    self.curses_interface.clear()
            except Exception as e:
                pass

    def cleanup(self) -> None:
        """Limpa recursos do console."""
        if self.curses_interface:
            try:
                if hasattr(self.curses_interface, "cleanup"):
                    self.curses_interface.cleanup()
            except Exception as e:
                pass


def create_console(use_curses: bool = True, curses_interface=None) -> IConsole:
    """
    Cria uma instância de console apropriada.

    Args:
        use_curses: Se True, tenta usar curses; caso contrário, usa console simples
        curses_interface: Interface curses opcional

    Returns:
        Instância de console
    """
    if use_curses and curses_interface:
        return CursesConsoleAdapter(curses_interface)
    else:
        return SimpleConsole()
