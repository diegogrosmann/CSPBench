# src/console_manager.py
"""
Gerenciador de console para centralizar e sincronizar a saída,
evitando que múltiplas threads corrompam a exibição.
"""
import sys
import threading
from typing import Optional

class ConsoleManager:
    _instance: Optional['ConsoleManager'] = None
    _lock = threading.Lock()

    def __new__(cls) -> 'ConsoleManager':
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super(ConsoleManager, cls).__new__(cls)
                    cls._instance.spinner_active = False
                    cls._instance.spinner_prefix = ""
                    cls._instance.spinner_animation = None
        return cls._instance

    def set_spinner(self, prefix: str, animation) -> None:
        with self._lock:
            self.spinner_prefix = prefix
            self.spinner_animation = animation

    def start_spinner(self) -> None:
        with self._lock:
            self.spinner_active = True

    def stop_spinner(self) -> None:
        with self._lock:
            self.spinner_active = False
            # Limpar linha quando parar spinner
            sys.stdout.write(f"\r{' ' * 70}")
            sys.stdout.write("\r")
            sys.stdout.flush()

    def print(self, message: str, end: str = "\n", flush: bool = False) -> None:
        with self._lock:
            sys.stdout.write(message + end)
            if flush:
                sys.stdout.flush()

    def print_warning(self, prefix: str, message: str) -> None:
        self.print(f"  AVISO: {message}")

console = ConsoleManager()
console = ConsoleManager()
