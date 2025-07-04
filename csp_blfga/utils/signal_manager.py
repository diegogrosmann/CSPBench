"""
Sistema centralizado de handlers de sinais para o CSP-BLFGA.

Classes:
    SignalManager: Gerenciador centralizado de sinais do sistema.

Funções:
    setup_signal_handlers(...): Configuração padrão de handlers.
    graceful_shutdown(...): Encerramento gracioso com limpeza.
"""

import logging
import signal
import sys
import threading
from collections.abc import Callable

logger = logging.getLogger(__name__)


class SignalManager:
    """
    Gerenciador centralizado de sinais do sistema.

    Permite registrar callbacks personalizados para diferentes sinais
    e garante encerramento gracioso do sistema.
    """

    def __init__(self):
        self.shutdown_callbacks: list[Callable[[], None]] = []
        self.interrupted = False
        self._lock = threading.Lock()

    def register_shutdown_callback(self, callback: Callable[[], None]):
        """
        Registra um callback para ser executado durante o shutdown.

        Args:
            callback: Função a ser chamada durante o encerramento
        """
        with self._lock:
            self.shutdown_callbacks.append(callback)

    def remove_shutdown_callback(self, callback: Callable[[], None]):
        """
        Remove um callback de shutdown.

        Args:
            callback: Função a ser removida
        """
        with self._lock:
            if callback in self.shutdown_callbacks:
                self.shutdown_callbacks.remove(callback)

    def signal_handler(self, signum: int, frame):
        """
        Handler centralizado para sinais de interrupção.

        Args:
            signum: Número do sinal recebido
            frame: Frame atual de execução
        """
        signal_name = signal.Signals(signum).name
        logger.info(f"Recebido sinal {signal_name} ({signum})")

        with self._lock:
            if self.interrupted:
                # Segunda interrupção - forçar saída
                print("\n⚠️ Segunda interrupção detectada. Forçando saída...")
                sys.exit(1)

            self.interrupted = True

        print(f"\n🛑 Recebido {signal_name}. Iniciando encerramento gracioso...")

        # Executar callbacks de shutdown
        self._execute_shutdown_callbacks()

        print("✅ Encerramento gracioso completo.")
        sys.exit(0)

    def _execute_shutdown_callbacks(self):
        """Executa todos os callbacks de shutdown registrados."""
        with self._lock:
            for callback in self.shutdown_callbacks:
                try:
                    callback()
                except Exception as e:
                    logger.error(f"Erro em callback de shutdown: {e}")

    def setup_handlers(self):
        """Configura os handlers padrão de sinais."""
        signal.signal(signal.SIGINT, self.signal_handler)  # Ctrl+C
        signal.signal(signal.SIGTERM, self.signal_handler)  # Terminação

        logger.info("Handlers de sinais configurados")

    def is_interrupted(self) -> bool:
        """Verifica se o sistema foi interrompido."""
        with self._lock:
            return self.interrupted


# Instância global do gerenciador
_signal_manager = SignalManager()


def get_signal_manager() -> SignalManager:
    """Retorna a instância global do gerenciador de sinais."""
    return _signal_manager


def setup_signal_handlers():
    """Configuração padrão de handlers de sinais."""
    _signal_manager.setup_handlers()


def register_shutdown_callback(callback: Callable[[], None]):
    """
    Registra um callback para ser executado durante o shutdown.

    Args:
        callback: Função a ser chamada durante o encerramento
    """
    _signal_manager.register_shutdown_callback(callback)


def graceful_shutdown(message: str | None = None):
    """
    Executa um encerramento gracioso do sistema.

    Args:
        message: Mensagem opcional para exibir
    """
    if message:
        print(f"\n{message}")

    _signal_manager._execute_shutdown_callbacks()
    print("✅ Encerramento gracioso completo.")
    sys.exit(0)


def is_interrupted() -> bool:
    """Verifica se o sistema foi interrompido por sinal."""
    return _signal_manager.is_interrupted()
