"""
Factory para criação de consoles.

Fornece funções para criar instâncias de consoles com fallback automático
e detecção de ambiente.
"""

import logging
import os
from typing import Union

from ..console import IConsole
from .simple_console import SimpleConsole
from .curses_console import CursesConsole


logger = logging.getLogger(__name__)


def create_console(prefer_curses: bool = True) -> IConsole:
    """
    Cria uma instância de console com fallback automático.

    Args:
        prefer_curses: Se True, tenta usar CursesConsole primeiro

    Returns:
        IConsole: Instância do console criado
    """
    # Verificar se estamos em um ambiente que suporta curses
    if prefer_curses and _can_use_curses():
        try:
            console = CursesConsole()
            logger.info("Console curses criado com sucesso")
            return console
        except Exception as e:
            logger.warning("Falha ao criar console curses: %s", e)

    # Fallback para console simples
    console = SimpleConsole()
    logger.info("Console simples criado")
    return console


def create_auto_console() -> IConsole:
    """
    Cria automaticamente o melhor console disponível.

    Returns:
        IConsole: Instância do console criado
    """
    return create_console(prefer_curses=True)


def create_simple_console() -> IConsole:
    """
    Cria explicitamente um console simples.

    Returns:
        IConsole: Instância do SimpleConsole
    """
    return SimpleConsole()


def create_curses_console() -> IConsole:
    """
    Cria explicitamente um console curses.

    Returns:
        IConsole: Instância do CursesConsole

    Raises:
        RuntimeError: Se curses não estiver disponível
    """
    if not _can_use_curses():
        raise RuntimeError("Curses não está disponível neste ambiente")

    return CursesConsole()


def _can_use_curses() -> bool:
    """
    Verifica se curses pode ser usado no ambiente atual.

    Returns:
        bool: True se curses estiver disponível
    """
    # Verificar se estamos em um terminal
    if not os.isatty(0):
        return False

    # Verificar se a variável TERM está definida
    if not os.environ.get("TERM"):
        return False

    # Verificar se não estamos em um ambiente que não suporta curses
    term = os.environ.get("TERM", "").lower()
    if term in ("dumb", "unknown"):
        return False

    # Verificar se curses está disponível
    try:
        import curses

        # Teste básico para verificar se curses funciona
        curses.setupterm()
        return True
    except (ImportError, curses.error):
        return False
    except Exception:
        return False
