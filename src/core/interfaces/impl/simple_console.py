"""
Implementação de console simples para IConsole.

Console baseado em prints padrão que implementa a interface IConsole
para ambientes onde curses não está disponível.
"""

import logging
from typing import Optional

from ..console import IConsole, TaskSlot

logger = logging.getLogger(__name__)


class SimpleConsole:
    """
    Console simples baseado em prints.

    Implementa IConsole usando prints padrão do Python,
    adequado para ambientes sem suporte a curses.
    """

    def __init__(self):
        """Inicializa o console simples."""
        self._active_slots = {}
        logger.debug("SimpleConsole inicializado")

    def show_task_start(self, name: str, slot: TaskSlot) -> None:
        """
        Mostra o início de uma tarefa.

        Args:
            name: Nome da tarefa
            slot: Slot da tarefa
        """
        slot.active = True
        slot.current_message = "Iniciando..."
        self._active_slots[slot.slot_id] = slot

        print(f"[{slot.slot_id:2d}] {name}: Iniciando...")

    def update_progress(self, slot: TaskSlot, message: str) -> None:
        """
        Atualiza o progresso de uma tarefa.

        Args:
            slot: Slot da tarefa
            message: Mensagem de progresso
        """
        slot.current_message = message

        # Atualizar apenas se a mensagem mudou significativamente
        if self._should_update_progress(slot, message):
            print(f"[{slot.slot_id:2d}] {slot.name}: {message}")

    def show_task_end(
        self, slot: TaskSlot, success: bool, message: str, error: Optional[str] = None
    ) -> None:
        """
        Mostra o fim de uma tarefa.

        Args:
            slot: Slot da tarefa
            success: True se bem-sucedida
            message: Mensagem final
            error: Mensagem de erro (opcional)
        """
        slot.active = False
        slot.current_message = message

        status = "✅" if success else "❌"
        print(f"[{slot.slot_id:2d}] {slot.name}: {status} {message}")

        if error:
            print(f"[{slot.slot_id:2d}] {slot.name}: Erro: {error}")

        # Remover do rastreamento
        if slot.slot_id in self._active_slots:
            del self._active_slots[slot.slot_id]

    def print(self, *args, **kwargs) -> None:
        """
        Imprime uma mensagem no console.

        Args:
            *args: Argumentos para print
            **kwargs: Argumentos nomeados para print
        """
        print(*args, **kwargs)

    def clear(self) -> None:
        """Limpa a tela (não implementado para console simples)."""
        # Console simples não limpa a tela
        pass

    def cleanup(self) -> None:
        """Limpa recursos do console."""
        self._active_slots.clear()
        logger.debug("SimpleConsole limpo")

    def _should_update_progress(self, slot: TaskSlot, message: str) -> bool:
        """
        Determina se deve atualizar o progresso.

        Args:
            slot: Slot da tarefa
            message: Nova mensagem

        Returns:
            bool: True se deve atualizar
        """
        # Atualizar se a mensagem mudou
        if slot.current_message != message:
            return True

        # Atualizar mensagens de tempo periodicamente
        if "executando" in message.lower() or "s)" in message.lower():
            return True

        return False
