#!/usr/bin/env python3
"""CSPBench Main Entry Point - Hexagonal Architecture.

This module serves as the main entry point for CSPBench, implementing a hexagonal
architecture pattern with dependency injection and modern CLI capabilities.

Features:
  - Load configuration (config/settings.yaml + ENV)
  - Dependency injection bootstrap
  - Modern CLI (Typer)
  - Interactive menu when no args

Key CLI usage (use .venv/bin/python):
  .venv/bin/python main.py
  .venv/bin/python main.py run <Algorithm> <dataset>
  .venv/bin/python main.py batch <file.yaml>
  .venv/bin/python main.py work-item-create <name> --config conf.yaml
  .venv/bin/python main.py work-item-list

Common environment variables:
  DATASET_DIRECTORY, LOG_LEVEL, OUTPUT_BASE_DIRECTORY, NCBI_EMAIL, NCBI_API_KEY
"""

import os
import sys
from pathlib import Path
from typing import Optional

import typer

# Early startup notice as soon as possible when running as script
if __name__ == "__main__":
    try:
        print(
            "🔧 Initializing CSPBench... please wait, loading modules and configurations.",
            flush=True,
        )
    except Exception:
        pass

# Load environment variables from .env file
try:
    from dotenv import load_dotenv

    load_dotenv(override=False)
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
    logger.info("Logging system initialized successfully")
    logger.info("Starting CSPBench - main mode")
except Exception as e:
    print(f"⚠️  Warning: failed to initialize logging: {e}")
    # Continue execution even if logging fails
    logger = None


def _ensure_runtime_directories():
    """Create (if needed) and report status for key runtime directories.

    This function ensures that essential runtime directories exist by checking
    environment variables and creating directories as needed. It handles three
    key directories: DATASET_DIRECTORY, BATCH_DIRECTORY, and OUTPUT_BASE_DIRECTORY.

    The function will:
    - Check if each environment variable is defined
    - Create the directory if it doesn't exist
    - Log the status of each directory operation
    - Continue execution even if some directories fail to be created

    Environment Variables:
        DATASET_DIRECTORY: Directory for storing datasets
        BATCH_DIRECTORY: Directory for batch processing files
        OUTPUT_BASE_DIRECTORY: Base directory for output files

    Note:
        This function is designed to be fault-tolerant and will not raise
        exceptions that could prevent the application from starting.
    """
    dir_env_vars = [
        "DATASET_DIRECTORY",
        "BATCH_DIRECTORY",
        "OUTPUT_BASE_DIRECTORY",
    ]

    results: list[str] = []
    for var in dir_env_vars:
        raw = os.getenv(var)
        if not raw:
            msg = f"⚠️  {var} not defined in environment (.env)."
            print(msg)
            if logger:
                logger.warning(msg)
            results.append(f"{var}: MISSING")
            continue

        try:
            p = Path(raw).expanduser()
            # Only create if it doesn't exist
            if not p.exists():
                p.mkdir(parents=True, exist_ok=True)
                msg = f"✅ {var} created: {p}"
                if logger:
                    logger.info(msg)
                results.append(f"{var}: CREATED")
            else:
                msg = f"ℹ️  {var} already exists: {p}"
                if logger:
                    logger.debug(msg)
                results.append(f"{var}: OK")
        except Exception as exc:  # noqa: BLE001
            msg = f"❌ Failed to ensure directory for {var}: {exc} (value='{raw}')"
            print(msg)
            if logger:
                logger.error(msg, exc_info=True)
            results.append(f"{var}: ERROR")

    # Compact final summary (useful in aggregated logs)
    if logger:
        logger.info("Base directories status: " + ", ".join(results))


# Execute verification/creation immediately after loading logging & .env
try:
    _ensure_runtime_directories()
except Exception as _e:  # noqa: BLE001
    # Never prevent initialization because of this
    print(f"⚠️  Warning: unexpected error ensuring base directories: {_e}")
    if logger:
        logger.warning(
            f"Unexpected error ensuring base directories: {_e}", exc_info=True
        )


# IMPORTANT: Import algorithms first to load the global_registry via auto-discovery
try:  # noqa: WPS501
    if logger:
        logger.debug("Loading algorithms package via auto-discovery")
    import algorithms  # noqa: F401

    if logger:
        logger.info("'algorithms' package loaded successfully")
except Exception as _e:  # noqa: BLE001
    # Failure to load algorithms should not prevent basic CLI
    error_msg = f"⚠️  Warning: could not load 'algorithms' package: {_e}"
    print(error_msg)
    if logger:
        logger.warning(f"Failed to load algorithms: {_e}", exc_info=True)

# Initialize WorkService and pause orphaned running work
try:
    if logger:
        logger.debug("Initializing global WorkService")
    from src.application.services.work_service import (
        initialize_work_service,
    )

    work_service = initialize_work_service()

    if logger:
        logger.info("WorkService initialized successfully")
except Exception as _e:  # noqa: BLE001
    # WorkService initialization failure should not prevent basic CLI
    error_msg = f"⚠️  Warning: could not initialize WorkService: {_e}"
    print(error_msg)
    if logger:
        logger.warning(f"Failed to initialize WorkService: {_e}", exc_info=True)

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
    """Dispatch parsed CLI arguments without interactive concerns.

    This function handles the execution of CLI commands by normalizing arguments
    and invoking the appropriate Typer command handlers.

    Args:
        args: List of command-line arguments to process

    Note:
        This function performs argument normalization, specifically converting
        standalone YAML files to 'batch <file>' commands for backward compatibility.

        The function temporarily modifies sys.argv to ensure proper Typer execution
        and restores it afterward.
    """
    import sys as _sys

    if logger:
        logger.debug(f"Dispatching CLI arguments: {args}")

    # Normalize invocation with only a YAML file to: batch <file>
    # (previously used '--batch' which doesn't exist as a Typer option)
    if len(args) == 1 and args[0].endswith((".yaml", ".yml")):
        if logger:
            logger.info(f"Normalizing YAML file to batch command: {args[0]}")
        args = ["batch", args[0]]

    original = _sys.argv[:]
    try:
        if logger:
            logger.debug(f"Executing Typer command with arguments: {args}")
        _sys.argv = ["main.py"] + args
        app()
    finally:
        _sys.argv = original
        if logger:
            logger.debug("sys.argv arguments restored")


def main(args: Optional[list] = None):
    """Programmatic entrypoint that does NOT show interactive menu on empty args.

    This function serves as the main programmatic interface for CSPBench,
    allowing the application to be invoked from other Python code.

    Args:
        args: Optional list of command arguments. If None, uses sys.argv[1:]

    Raises:
        SystemExit: When the application needs to exit with a specific code
        KeyboardInterrupt: When user interrupts execution (handled gracefully)

    Note:
        Unlike the __main__ execution path, this function will not show the
        interactive menu when no arguments are provided.
    """
    import sys

    if args is None:
        args = sys.argv[1:]

    if logger:
        logger.info(f"Programmatic entry point started with arguments: {args}")

    try:
        _dispatch_args(args)
        if logger:
            logger.info("Programmatic execution completed successfully")
    except KeyboardInterrupt:
        error_msg = "\n🚫 Operation cancelled by user (Ctrl+C)"
        print(error_msg)
        if logger:
            logger.warning("Operation cancelled by user (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        if logger:
            logger.debug("SystemExit caught - propagating")
        raise
    except Exception as e:
        error_msg = f"❌ Unexpected error: {e}"
        print(error_msg)
        if logger:
            logger.error(f"Unexpected error in execution: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    import sys

    try:
        if len(sys.argv) == 1:
            if logger:
                logger.info("No arguments provided - starting interactive menu")

            def _invoke(cmd_args: list[str]):  # local to avoid exporting
                """Local function to invoke commands from interactive menu.

                Args:
                    cmd_args: List of command arguments from interactive menu
                """
                if logger:
                    logger.debug(f"Interactive menu invoking command: {cmd_args}")
                _dispatch_args(cmd_args)

            show_interactive_menu(_invoke)
        else:
            if logger:
                logger.info(f"CLI arguments provided: {sys.argv[1:]} - direct mode")
            _dispatch_args(sys.argv[1:])
    except KeyboardInterrupt:
        error_msg = "\n🚫 Operation cancelled by user (Ctrl+C)"
        print(error_msg)
        if logger:
            logger.warning("Operation cancelled by user (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        if logger:
            logger.debug("SystemExit caught in main - propagating")
        raise
    except Exception as e:
        error_msg = f"❌ Unexpected error: {e}"
        print(error_msg)
        if logger:
            logger.error(f"Unexpected error in main: {e}", exc_info=True)
        sys.exit(1)
