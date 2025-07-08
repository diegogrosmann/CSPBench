"""
Factory para criação de consoles de execução.

Permite criar diferentes tipos de consoles baseados no ambiente
e preferências do usuário.
"""

import logging
from typing import Any, Optional

from ..execution_console import IConsole
from .execution_curses_console import CursesConsole
from .execution_simple_console import SimpleConsole

logger = logging.getLogger(__name__)


def create_console(console_type: str = "auto", **kwargs) -> IConsole:
    """
    Cria um console baseado no tipo especificado.

    Args:
        console_type: Tipo do console ("simple", "curses", "auto")
        **kwargs: Argumentos adicionais para o console

    Returns:
        IConsole: Instância do console criado

    Raises:
        ValueError: Se o tipo do console não for suportado
    """
    if console_type == "simple":
        return create_simple_console(**kwargs)
    elif console_type == "curses":
        return create_curses_console(**kwargs)
    elif console_type == "auto":
        return create_auto_console(**kwargs)
    else:
        raise ValueError(f"Tipo de console não suportado: {console_type}")


def create_simple_console(verbose: bool = True) -> SimpleConsole:
    """
    Cria um console simples baseado em prints.

    Args:
        verbose: Se True, mostra mensagens detalhadas

    Returns:
        SimpleConsole: Instância do console simples
    """
    logger.info("Criando SimpleConsole")
    return SimpleConsole(verbose=verbose)


def create_curses_console(stdscr: Optional[Any] = None) -> CursesConsole:
    """
    Cria um console curses.

    Args:
        stdscr: Tela curses (opcional)

    Returns:
        CursesConsole: Instância do console curses
    """
    logger.info("Criando CursesConsole")
    return CursesConsole(stdscr=stdscr)


def create_auto_console(**kwargs) -> IConsole:
    """
    Cria um console automaticamente baseado no ambiente.

    Args:
        **kwargs: Argumentos adicionais

    Returns:
        IConsole: Instância do console criado
    """
    # Verificar se curses está disponível
    try:
        import curses

        # Verificar se está em um terminal adequado
        if hasattr(curses, "initscr"):
            logger.info("Ambiente suporta curses, criando CursesConsole")
            return create_curses_console(**kwargs)
    except ImportError:
        logger.info("Curses não disponível")
    except Exception as e:
        logger.info("Erro ao verificar curses: %s", e)

    # Fallback para console simples
    logger.info("Usando SimpleConsole como fallback")
    return create_simple_console(**kwargs)
