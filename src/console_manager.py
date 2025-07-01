# src/console_manager.py
"""
Gerenciador de console para centralizar e sincronizar a saída,
evitando que múltiplas threads corrompam a exibição.

Classes:
    ConsoleManager: Singleton thread-safe para saída no console.

Atributos:
    console (ConsoleManager): Instância global para uso em toda a aplicação.
"""
import sys
import threading
from typing import Optional

class ConsoleManager:
    """
    Gerencia a saída no console de forma thread-safe, evitando conflitos
    entre múltiplas threads e centralizando prints e animações.

    Métodos:
        print(message, end, flush): Imprime mensagem thread-safe.
        print_warning(prefix, message): Imprime aviso formatado.
        print_inline(message, flush): Imprime mensagem na mesma linha.
    """
    _instance: Optional['ConsoleManager'] = None
    _lock = threading.Lock()

    def __new__(cls) -> 'ConsoleManager':
        """
        Garante singleton thread-safe.

        Returns:
            ConsoleManager: Instância única da classe.
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super(ConsoleManager, cls).__new__(cls)
        return cls._instance

    def print(self, message: str, end: str = "\n", flush: bool = False) -> None:
        """
        Imprime mensagem no console de forma thread-safe.

        Args:
            message (str): Mensagem a ser impressa.
            end (str): Sufixo (default: '\n').
            flush (bool): Força flush imediato.
        """
        with self._lock:
            # Limpar a linha antes de imprimir uma nova mensagem completa
            if end == "\n":
                sys.stdout.write(f"\r{' ' * 70}\r")
            
            sys.stdout.write(message + end)
            if flush:
                sys.stdout.flush()

    def print_inline(self, message: str, flush: bool = True) -> None:
        """
        Imprime mensagem na mesma linha (inline), thread-safe.

        Args:
            message (str): Mensagem a ser impressa.
            flush (bool): Força flush imediato.
        """
        with self._lock:
            sys.stdout.write(f"\r{message}")
            if flush:
                sys.stdout.flush()

    def print_warning(self, prefix: str, message: str) -> None:
        """
        Imprime aviso formatado no console.

        Args:
            prefix (str): Prefixo do aviso (não utilizado).
            message (str): Mensagem do aviso.
        """
        self.print(f"  AVISO: {message}")

console = ConsoleManager()
