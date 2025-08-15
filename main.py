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


# IMPORTANT: Import algorithms first to load the global_registry via auto-discovery
try:  # noqa: WPS501
    import algorithms  # noqa: F401
except Exception as _e:  # noqa: BLE001
    # Falha ao carregar algoritmos não deve impedir CLI básica
    print(f"⚠️  Aviso: não foi possível carregar pacote 'algorithms': {_e}")
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

    # Normalize invocation with only a YAML file to: batch <file>
    # (anteriormente usava '--batch' que não existe como opção Typer)
    if len(args) == 1 and args[0].endswith((".yaml", ".yml")):
        args = ["batch", args[0]]

    original = _sys.argv[:]
    try:
        _sys.argv = ["main.py"] + args
        app()
    finally:
        _sys.argv = original


def main(args: Optional[list] = None):
    """Programmatic entrypoint (does NOT show interactive menu on empty args)."""
    import sys

    if args is None:
        args = sys.argv[1:]
    try:
        _dispatch_args(args)
    except KeyboardInterrupt:
        print("\n🚫 Operation cancelled by user (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        raise
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    import sys

    try:
        if len(sys.argv) == 1:

            def _invoke(cmd_args: list[str]):  # local to avoid exporting
                _dispatch_args(cmd_args)

            show_interactive_menu(_invoke)
        else:
            _dispatch_args(sys.argv[1:])
    except KeyboardInterrupt:
        print("\n🚫 Operation cancelled by user (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        raise
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)
