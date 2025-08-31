"""Console interactive menu for selecting and executing batch files.

Separated from main.py to keep the entrypoint slim and focused only on
argument routing and bootstrapping.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional

from typing import Callable
from src.infrastructure.utils.path_utils import get_batch_directory

# Recebe injetado um callable que representa o app Typer; evitamos import circular.
AppInvoker = Callable[[list[str]], None]


def show_interactive_menu(app_invoke: AppInvoker) -> None:
    """Show interactive interface for batch file selection."""
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
        # Aciona comando 'batch' diretamente (antes usava '--batch' inexistente)
        app_invoke(["batch", str(selected_file)])


def _show_commands_help(message: str) -> None:
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
    print("📋 Available batch files:")
    print()
    for i, batch_file in enumerate(batch_files, 1):
        description = _extract_file_description(batch_file)
        print(f"  {i}. {batch_file.name}")
        print(f"     {description}")
        print()


def _extract_file_description(batch_file: Path) -> str:
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
    print("📋 Available manual commands:")
    print("  test         - Basic system test")
    print("  algorithms   - List available algorithms")
    print("  config-info  - Show configuration")
    print("  run          - Execute single algorithm")
    print()


def _get_user_selection(batch_files: List[Path]) -> Optional[Path]:
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
