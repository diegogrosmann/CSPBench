"""Console interactive menu for selecting and executing batch files.

This module provides an interactive terminal interface for selecting and executing
CSPBench batch configuration files. It's separated from main.py to keep the
entrypoint slim and focused only on argument routing and bootstrapping.

The module provides functionality to:
- List available batch files from the configured directory
- Extract descriptions from batch file comments
- Present manual command options
- Handle user selection and execute the chosen batch file
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, List, Optional

from src.infrastructure.utils.path_utils import get_batch_directory

# Receives injected callable representing Typer app; avoids circular import
AppInvoker = Callable[[list[str]], None]


def show_interactive_menu(app_invoke: AppInvoker) -> None:
    """Show interactive interface for batch file selection.

    This function displays an interactive menu that allows users to:
    1. Browse available batch files in the configured directory
    2. View file descriptions extracted from comments
    3. Select and execute a batch file
    4. Access manual command options

    Args:
        app_invoke: Callable that can invoke Typer app commands with arguments
    """
    print("CSPBench v0.1.0 - Framework for Closest String Problem")
    print("=" * 60)
    print()

    batches_dir = get_batch_directory()
    if not batches_dir.exists():
        _show_commands_help(f"📁 Directory '{batches_dir}' not found")
        return

    batch_files: List[Path] = list(batches_dir.glob("*.yaml")) + list(
        batches_dir.glob("*.yml")
    )

    if not batch_files:
        _show_commands_help(f"📋 No batch files found in '{batches_dir}'")
        return

    _display_batch_files(batch_files)
    _display_manual_commands()

    selected_file = _get_user_selection(batch_files)
    if selected_file:
        # Execute 'batch' command directly (previously used non-existent '--batch')
        app_invoke(["batch", str(selected_file)])


def _show_commands_help(message: str) -> None:
    """Display help message with available commands when no batch files are found.

    Args:
        message: Message to display before the commands list
    """
    print(message)
    print()
    print("Available commands:")
    print("  test         - Basic system test")
    print("  algorithms   - List available algorithms")
    print("  config-info  - Show configuration")
    print("  run          - Execute single algorithm")
    print()
    print("Usage: python main.py <command> [arguments]")
    print("Example: python main.py test")


def _display_batch_files(batch_files: List[Path]) -> None:
    """Display numbered list of available batch files with descriptions.

    Args:
        batch_files: List of batch file paths to display
    """
    print("📋 Available batch files:")
    print()
    for i, batch_file in enumerate(batch_files, 1):
        description = _extract_file_description(batch_file)
        print(f"  {i}. {batch_file.name}")
        print(f"     {description}")
        print()


def _extract_file_description(batch_file: Path) -> str:
    """Extract description from batch file comments.

    Searches for the first meaningful comment line in the file that starts with '#'
    and doesn't contain separator characters like '='.

    Args:
        batch_file: Path to the batch file to analyze

    Returns:
        Description string extracted from comments, or default message
    """
    try:
        with open(batch_file, encoding="utf-8") as f:
            first_chars = f.read(400)
        for line in first_chars.splitlines():
            line = line.strip()
            if line.startswith("#") and len(line) > 1:
                text = line[1:].strip()
                if text and not text.startswith("="):
                    return text
        return "Batch configuration file"
    except Exception:
        return "Batch configuration file"


def _display_manual_commands() -> None:
    """Display available manual commands that can be executed."""
    print("📋 Available manual commands:")
    print("  test         - Basic system test")
    print("  algorithms   - List available algorithms")
    print("  config-info  - Show configuration")
    print("  run          - Execute single algorithm")
    print()


def _get_user_selection(batch_files: List[Path]) -> Optional[Path]:
    """Get user selection from available batch files.

    Prompts the user to select a batch file by number and validates the input.

    Args:
        batch_files: List of available batch files

    Returns:
        Selected batch file path, or None if user cancels or provides invalid input
    """
    try:
        choice = input("💡 Select a file (number) or press Enter to exit: ").strip()
        if choice == "":
            print("👋 Goodbye!")
            return None
        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(batch_files):
                return batch_files[idx - 1]
            print("❌ Invalid number!")
            return None
        print("❌ Invalid input!")
        return None
    except (KeyboardInterrupt, EOFError):
        print("\n👋 Goodbye!")
        return None
