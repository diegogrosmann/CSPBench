"""
Console Manager baseado em curses para CSP-BLFGA.

Este módulo implementa uma interface de console avançada usando curses
para melhor experiência do usuário em terminais compatíveis.

Classes:
    BaseConsoleManager: Interface base para console managers
    CursesConsole: Console manager baseado em curses
    SimpleConsole: Console manager simples (fallback)

Funções:
    detect_tty_support(): Detecta se o terminal suporta curses
    create_console_manager(): Factory para criar console apropriado
"""

import curses
import os
import sys
import threading
from abc import ABC, abstractmethod


class BaseConsoleManager(ABC):
    """
    Interface base para console managers.

    Define os métodos que devem ser implementados por todos os
    console managers para garantir compatibilidade.
    """

    @abstractmethod
    def print(self, *args, **kwargs):
        """Imprime uma mensagem no console."""
        raise NotImplementedError

    @abstractmethod
    def clear(self):
        """Limpa a tela."""
        raise NotImplementedError

    @abstractmethod
    def set_title(self, title: str):
        """Define título da janela."""
        raise NotImplementedError

    @abstractmethod
    def show_progress(self, message: str, current: int, total: int):
        """Mostra progresso."""
        raise NotImplementedError

    @abstractmethod
    def cleanup(self):
        """Limpa recursos do console."""
        raise NotImplementedError

    @abstractmethod
    def hide_progress(self):
        """Esconde progresso."""
        raise NotImplementedError

    @abstractmethod
    def show_status(self, message: str):
        """Mostra status."""
        raise NotImplementedError


def detect_tty_support() -> bool:
    """
    Detecta se o terminal suporta curses/TTY.

    Returns:
        bool: True se suporta curses, False caso contrário
    """
    # Verificar se está rodando em um TTY
    if not sys.stdout.isatty():
        return False

    # Verificar variáveis de ambiente que indicam terminal não-interativo
    if os.environ.get("TERM") == "dumb":
        return False

    if os.environ.get("CI") == "true":
        return False

    # Tentar inicializar curses
    try:
        _ = curses.initscr()  # Variável temporária para evitar unused variable
        curses.endwin()
        return True
    except Exception:  # pylint: disable=broad-except
        return False


class SimpleConsole(BaseConsoleManager):
    """
    Console manager simples que funciona como fallback.

    Este é usado quando curses não está disponível ou
    quando o terminal não suporta funcionalidades avançadas.
    """

    def __init__(self):
        self.lock = threading.Lock()

    def print(self, *args, **kwargs):
        """Print simples thread-safe."""
        with self.lock:
            print(*args, **kwargs)

    def clear(self):
        """Limpa a tela (no-op para console simples)."""

    def set_title(self, title: str):
        """Define título (no-op para console simples)."""

    def show_progress(self, message: str, current: int, total: int):
        """Mostra progresso simples."""
        if total > 0:
            percentage = (current / total) * 100
            self.print(
                f"\r{message}: {percentage:.1f}% ({current}/{total})",
                end="",
                flush=True,
            )

    def hide_progress(self):
        """Esconde progresso (no-op para console simples)."""

    def show_status(self, message: str):
        """Mostra status."""
        self.print(f"Status: {message}")

    def cleanup(self):
        """Limpa recursos (no-op para console simples)."""


class CursesConsole(BaseConsoleManager):
    """
    Console manager avançado baseado em curses.

    Fornece funcionalidades como:
    - Janelas múltiplas
    - Barras de progresso
    - Status em tempo real
    - Cores e formatação
    """

    def __init__(self):
        self.stdscr = None
        self.initialized = False
        self.lock = threading.Lock()
        self.progress_win = None
        self.status_win = None
        self.main_win = None
        self.current_y = 0
        self.max_y = 0
        self.max_x = 0
        self._init_curses()

    def _init_curses(self):
        """Inicializa curses."""
        try:
            self.stdscr = curses.initscr()
            curses.noecho()
            curses.cbreak()
            self.stdscr.keypad(True)

            # Configurar cores se disponível
            if curses.has_colors():
                curses.start_color()
                curses.init_pair(1, curses.COLOR_WHITE, curses.COLOR_BLACK)
                curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)
                curses.init_pair(3, curses.COLOR_YELLOW, curses.COLOR_BLACK)
                curses.init_pair(4, curses.COLOR_RED, curses.COLOR_BLACK)
                curses.init_pair(5, curses.COLOR_BLUE, curses.COLOR_BLACK)

            self.max_y, self.max_x = self.stdscr.getmaxyx()

            # Criar janelas
            self._create_windows()

            self.initialized = True

        except Exception:  # pylint: disable=broad-except
            # Se falhar, usar console simples
            self.initialized = False
            if self.stdscr:
                curses.endwin()
                self.stdscr = None

    def _create_windows(self):
        """Cria as janelas do curses."""
        if not self.stdscr:
            return

        # Janela principal (maior parte da tela)
        main_height = self.max_y - 4
        self.main_win = curses.newwin(main_height, self.max_x, 0, 0)
        self.main_win.scrollok(True)

        # Janela de progresso (2 linhas)
        self.progress_win = curses.newwin(2, self.max_x, main_height, 0)

        # Janela de status (2 linhas)
        self.status_win = curses.newwin(2, self.max_x, main_height + 2, 0)

        # Desenhar bordas
        self.progress_win.box()
        self.status_win.box()

        self.stdscr.refresh()

    def print(self, *args, **kwargs):
        """Print com suporte a curses."""
        if not self.initialized or not self.main_win:
            # Fallback para print normal
            print(*args, **kwargs)
            return

        with self.lock:
            try:
                # Converter argumentos para string
                message = " ".join(str(arg) for arg in args)

                # Adicionar à janela principal
                self.main_win.addstr(f"{message}\n")
                self.main_win.refresh()

            except Exception:  # pylint: disable=broad-except
                # Fallback para print normal em caso de erro
                print(*args, **kwargs)

    def clear(self):
        """Limpa a tela."""
        if not self.initialized:
            os.system("clear" if os.name == "posix" else "cls")
            return

        with self.lock:
            try:
                if self.main_win:
                    self.main_win.clear()
                    self.main_win.refresh()
            except Exception:  # pylint: disable=broad-except
                ...

    def set_title(self, title: str):
        """Define título da janela."""
        if not self.initialized or not self.stdscr:
            return

        with self.lock:
            try:
                # Mostrar título na primeira linha
                self.stdscr.addstr(0, 2, f" {title} ", curses.color_pair(2))
                self.stdscr.refresh()
            except Exception:  # pylint: disable=broad-except
                ...

    def show_progress(self, message: str, current: int, total: int):
        """Mostra barra de progresso."""
        if not self.initialized or not self.progress_win:
            # Fallback para progresso simples
            if total > 0:
                percentage = (current / total) * 100
                print(
                    f"\r{message}: {percentage:.1f}% ({current}/{total})",
                    end="",
                    flush=True,
                )
            return

        with self.lock:
            try:
                self.progress_win.clear()
                self.progress_win.box()

                # Calcular progresso
                if total > 0:
                    percentage = (current / total) * 100
                    bar_width = self.max_x - 20
                    filled_width = int((current / total) * bar_width)

                    # Desenhar barra
                    progress_text = f"{message}: {percentage:.1f}%"
                    self.progress_win.addstr(1, 2, progress_text[: self.max_x - 4])

                    # Barra visual
                    bar_y = 1
                    bar_x = 2
                    for i in range(bar_width):
                        if i < filled_width:
                            self.progress_win.addch(
                                bar_y,
                                bar_x + len(progress_text) + 2 + i,
                                "█",
                                curses.color_pair(2),
                            )
                        else:
                            self.progress_win.addch(
                                bar_y,
                                bar_x + len(progress_text) + 2 + i,
                                "░",
                                curses.color_pair(1),
                            )

                self.progress_win.refresh()

            except Exception:  # pylint: disable=broad-except
                ...

    def hide_progress(self):
        """Esconde barra de progresso."""
        if not self.initialized or not self.progress_win:
            return

        with self.lock:
            try:
                self.progress_win.clear()
                self.progress_win.box()
                self.progress_win.refresh()
            except Exception:  # pylint: disable=broad-except
                ...

    def show_status(self, message: str):
        """Mostra status na linha de status."""
        if not self.initialized or not self.status_win:
            self.print(f"Status: {message}")
            return

        with self.lock:
            try:
                self.status_win.clear()
                self.status_win.box()
                self.status_win.addstr(1, 2, f"Status: {message}"[: self.max_x - 4])
                self.status_win.refresh()
            except Exception:  # pylint: disable=broad-except
                ...

    def cleanup(self):
        """Limpa recursos do curses."""
        if self.initialized and self.stdscr:
            try:
                curses.endwin()
            except Exception:  # pylint: disable=broad-except
                ...

    def __del__(self):
        """Destrutor para garantir limpeza."""
        self.cleanup()


def create_console_manager(force_simple: bool = False) -> BaseConsoleManager:
    """
    Factory para criar o console manager apropriado.

    Args:
        force_simple: Se True, força uso do console simples

    Returns:
        BaseConsoleManager: Instância do console apropriado
    """
    if force_simple or not detect_tty_support():
        return SimpleConsole()
    else:
        # Tentar criar curses console, fallback para simples
        try:
            console = CursesConsole()
            if console.initialized:
                return console
            else:
                return SimpleConsole()
        except Exception:  # pylint: disable=broad-except
            return SimpleConsole()
