#!/usr/bin/env python3
"""
CSPBench Main Entry Point - Hexagonal Architecture

Features:
  - Load configuration (config/settings.yaml + ENV)
  - Dependency injection bootstrap
  - Modern CLI (Typer)
  - Interactive menu when no args

Key CLI (usar .venv/bin/python):
  .venv/bin/python main.py
  .venv/bin/python main.py run <Algorithm> <dataset>
  .venv/bin/python main.py batch <arquivo.yaml>
  .venv/bin/python main.py work-item-create <name> --config conf.yaml
  .venv/bin/python main.py work-item-list

Variáveis de ambiente comuns:
  DATASET_PATH, LOG_LEVEL, OUTPUT_BASE_DIRECTORY, NCBI_EMAIL, NCBI_API_KEY
"""

import sys
from pathlib import Path
from typing import Optional
import typer

# Early startup notice as soon as possible when running as script
if __name__ == "__main__":
    try:
        print(
            "🔧 Inicializando CSPBench... aguarde, carregando módulos e configurações.",
            flush=True,
        )
    except Exception:
        pass

# Load environment variables from .env file
try:
    from dotenv import load_dotenv

    load_dotenv(override=True)
except ImportError:
    # If python-dotenv is not installed, continue without it
    pass

# Add the root directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

# Initialize logging system early using environment variables
try:
    from src.infrastructure.logging_config import LoggerConfig

    LoggerConfig.initialize()
    logger = LoggerConfig.get_logger("CSPBench.Main")
    logger.info("Sistema de logging inicializado com sucesso")
    logger.info("Iniciando CSPBench - modo principal")
except Exception as e:
    print(f"⚠️  Aviso: falha ao inicializar logging: {e}")
    # Continue execution even if logging fails
    logger = None


# IMPORTANT: Import algorithms first to load the global_registry via auto-discovery
try:  # noqa: WPS501
    if logger:
        logger.debug("Carregando pacote de algoritmos via auto-discovery")
    import algorithms  # noqa: F401

    if logger:
        logger.info("Pacote 'algorithms' carregado com sucesso")
except Exception as _e:  # noqa: BLE001
    # Falha ao carregar algoritmos não deve impedir CLI básica
    error_msg = f"⚠️  Aviso: não foi possível carregar pacote 'algorithms': {_e}"
    print(error_msg)
    if logger:
        logger.warning(f"Falha ao carregar algoritmos: {_e}", exc_info=True)

# Initialize WorkService and pause orphaned running work
try:
    if logger:
        logger.debug("Inicializando WorkService global")
    from src.application.services.work_service import initialize_work_service, pause_orphaned_running_work
    
    work_service = initialize_work_service()
    pause_orphaned_running_work(work_service)
    
    if logger:
        logger.info("WorkService inicializado e trabalhos órfãos pausados")
except Exception as _e:  # noqa: BLE001
    # Falha ao inicializar WorkService não deve impedir CLI básica
    error_msg = f"⚠️  Aviso: não foi possível inicializar WorkService: {_e}"
    print(error_msg)
    if logger:
        logger.warning(f"Falha ao inicializar WorkService: {_e}", exc_info=True)

from src.presentation.cli.commands import register_commands
from src.presentation.cli.interactive_menu import show_interactive_menu

app = typer.Typer(
    name="cspbench",
    help="CSPBench - Framework for Closest String Problem",
    add_completion=False,
    no_args_is_help=False,
    rich_markup_mode=None,
)

register_commands(app)


def _dispatch_args(args: list[str]) -> None:
    """Dispatch parsed CLI arguments (expects no interactive concern)."""
    import sys as _sys

    if logger:
        logger.debug(f"Despachando argumentos CLI: {args}")

    # Normalize invocation with only a YAML file to: batch <file>
    # (anteriormente usava '--batch' que não existe como opção Typer)
    if len(args) == 1 and args[0].endswith((".yaml", ".yml")):
        if logger:
            logger.info(f"Normalizando arquivo YAML para comando batch: {args[0]}")
        args = ["batch", args[0]]

    original = _sys.argv[:]
    try:
        if logger:
            logger.debug(f"Executando comando Typer com argumentos: {args}")
        _sys.argv = ["main.py"] + args
        app()
    finally:
        _sys.argv = original
        if logger:
            logger.debug("Argumentos sys.argv restaurados")


def main(args: Optional[list] = None):
    """Programmatic entrypoint (does NOT show interactive menu on empty args)."""
    import sys

    if args is None:
        args = sys.argv[1:]

    if logger:
        logger.info(f"Ponto de entrada programático iniciado com argumentos: {args}")

    try:
        _dispatch_args(args)
        if logger:
            logger.info("Execução programática concluída com sucesso")
    except KeyboardInterrupt:
        error_msg = "\n🚫 Operation cancelled by user (Ctrl+C)"
        print(error_msg)
        if logger:
            logger.warning("Operação cancelada pelo usuário (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        if logger:
            logger.debug("SystemExit capturado - propagando")
        raise
    except Exception as e:
        error_msg = f"❌ Unexpected error: {e}"
        print(error_msg)
        if logger:
            logger.error(f"Erro inesperado na execução: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    import sys

    try:
        if len(sys.argv) == 1:
            if logger:
                logger.info("Nenhum argumento fornecido - iniciando menu interativo")

            def _invoke(cmd_args: list[str]):  # local to avoid exporting
                if logger:
                    logger.debug(f"Menu interativo invocando comando: {cmd_args}")
                _dispatch_args(cmd_args)

            show_interactive_menu(_invoke)
        else:
            if logger:
                logger.info(f"Argumentos CLI fornecidos: {sys.argv[1:]} - modo direto")
            _dispatch_args(sys.argv[1:])
    except KeyboardInterrupt:
        error_msg = "\n🚫 Operation cancelled by user (Ctrl+C)"
        print(error_msg)
        if logger:
            logger.warning("Operação cancelada pelo usuário (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        if logger:
            logger.debug("SystemExit capturado no main - propagando")
        raise
    except Exception as e:
        error_msg = f"❌ Unexpected error: {e}"
        print(error_msg)
        if logger:
            logger.error(f"Erro inesperado no main: {e}", exc_info=True)
        sys.exit(1)
