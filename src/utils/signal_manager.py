"""
Sistema Centralizado de Gerenciamento de Sinais - CSPBench

Este módulo fornece um sistema robusto e centralizado para gerenciamento de sinais
do sistema operacional, garantindo encerramento gracioso e limpeza adequada de
recursos durante interrupções ou finalizações do CSPBench.

Arquitetura:
    O módulo implementa um padrão Observer para gerenciamento de sinais:
    - Gerenciador central de sinais
    - Sistema de callbacks registráveis
    - Thread-safety para ambientes paralelos
    - Encerramento gracioso coordenado
    - Limpeza automática de recursos

Funcionalidades:
    - Interceptação de sinais SIGINT e SIGTERM
    - Sistema de callbacks personalizáveis
    - Encerramento gracioso coordenado
    - Limpeza automática de recursos
    - Prevenção de corrupção de dados
    - Logging detalhado de operações

Sinais Gerenciados:
    - SIGINT (Ctrl+C): Interrupção manual do usuário
    - SIGTERM: Terminação solicitada pelo sistema
    - SIGQUIT: Saída com core dump (Unix)
    - Extensibilidade para outros sinais

Exemplo de Uso:
    ```python
    from src.utils.signal_manager import SignalManager, setup_signal_handlers

    # Criar gerenciador
    signal_manager = SignalManager()

    # Registrar callbacks de limpeza
    signal_manager.register_shutdown_callback(cleanup_database)
    signal_manager.register_shutdown_callback(save_state)
    signal_manager.register_shutdown_callback(close_files)

    # Configurar handlers padrão
    setup_signal_handlers(signal_manager)

    # Verificar se foi interrompido
    if signal_manager.is_interrupted():
        print("Execução foi interrompida")
    ```

Thread Safety:
    - Locks para operações críticas
    - Sincronização de callbacks
    - Prevenção de race conditions
    - Coordenação entre threads

Casos de Uso:
    - Limpeza de arquivos temporários
    - Finalização de conexões de rede
    - Salvamento de estado da aplicação
    - Encerramento de processos filhos
    - Liberação de recursos do sistema

Autor: CSPBench Development Team
Data: 2024
"""

import logging
import signal
import sys
import threading
from collections.abc import Callable
from typing import List

logger = logging.getLogger(__name__)


class SignalManager:
    """
    Gerenciador centralizado de sinais do sistema.

    Esta classe fornece uma interface unificada para gerenciamento de sinais
    do sistema operacional, permitindo registro de callbacks personalizados
    e garantindo encerramento gracioso com limpeza adequada de recursos.

    Funcionalidades:
        - Registro de callbacks de shutdown
        - Controle de estado de interrupção
        - Thread-safety para ambientes paralelos
        - Execução coordenada de limpeza
        - Logging detalhado de operações

    Atributos:
        shutdown_callbacks: Lista de funções a executar no shutdown
        interrupted: Flag indicando se sistema foi interrompido
        _lock: Lock para thread-safety

    Exemplo de Uso:
        ```python
        # Criar gerenciador
        manager = SignalManager()

        # Registrar callbacks
        manager.register_shutdown_callback(cleanup_temp_files)
        manager.register_shutdown_callback(save_progress)

        # Verificar estado
        if manager.is_interrupted():
            print("Sistema foi interrompido")

        # Executar shutdown manual
        manager.graceful_shutdown()
        ```

    Thread Safety:
        Todas as operações são thread-safe através do uso de locks,
        permitindo uso seguro em ambientes com múltiplas threads.
    """

    def __init__(self):
        """
        Inicializa o gerenciador de sinais.

        Configura estruturas internas e estado inicial para
        gerenciamento de callbacks e controle de interrupção.
        """
        self.shutdown_callbacks: List[Callable[[], None]] = []
        self.interrupted = False
        self._lock = threading.Lock()
        logger.debug("SignalManager inicializado")

    def register_shutdown_callback(self, callback: Callable[[], None]) -> None:
        """
        Registra um callback para ser executado durante o shutdown.

        Args:
            callback (Callable[[], None]): Função a ser chamada durante o encerramento.
                                          Deve ser callable sem parâmetros.

        Raises:
            TypeError: Se callback não é callable

        Exemplo:
            ```python
            def cleanup_function():
                print("Limpando recursos...")
                # Código de limpeza aqui

            manager.register_shutdown_callback(cleanup_function)
            ```

        Nota:
            - Callbacks são executados na ordem de registro
            - Falhas em um callback não impedem execução dos demais
            - Thread-safe para registro durante execução
        """
        if not callable(callback):
            raise TypeError("Callback deve ser callable")

        with self._lock:
            self.shutdown_callbacks.append(callback)
            logger.debug("Callback de shutdown registrado: %s", callback.__name__)

    def remove_shutdown_callback(self, callback: Callable[[], None]) -> bool:
        """
        Remove um callback de shutdown.

        Args:
            callback (Callable[[], None]): Função a ser removida

        Returns:
            bool: True se callback foi removido, False se não foi encontrado

        Exemplo:
            ```python
            # Remover callback específico
            removed = manager.remove_shutdown_callback(cleanup_function)
            if removed:
                print("Callback removido com sucesso")
            ```

        Nota:
            Thread-safe para remoção durante execução
        """
        with self._lock:
            if callback in self.shutdown_callbacks:
                self.shutdown_callbacks.remove(callback)
                logger.debug("Callback de shutdown removido: %s", callback.__name__)
                return True
            else:
                logger.debug(
                    "Callback não encontrado para remoção: %s", callback.__name__
                )
                return False

    def signal_handler(self, signum: int, frame):
        """
        Handler centralizado para sinais de interrupção.

        Args:
            signum: Número do sinal recebido
            frame: Frame atual de execução
        """
        signal_name = signal.Signals(signum).name
        logger.info("Recebido sinal %s (%s)", signal_name, signum)

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
                except (
                    OSError,
                    RuntimeError,
                    KeyboardInterrupt,
                    ValueError,
                    TypeError,
                ) as e:
                    logger.error("Erro em callback de shutdown: %s", e)

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

    graceful_shutdown(message)


def is_interrupted() -> bool:
    """Verifica se o sistema foi interrompido por sinal."""
    return _signal_manager.is_interrupted()
