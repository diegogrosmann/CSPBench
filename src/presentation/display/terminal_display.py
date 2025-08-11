"""Terminal display adapter for progress events."""

import sys
from typing import Optional, TYPE_CHECKING, Any
from datetime import datetime
import logging

if TYPE_CHECKING:
    from src.application.monitoring.progress_events import (
        ProgressEvent,
        AlgorithmProgressEvent,
        ExecutionProgressEvent,
    )


class TerminalDisplay:
    """
    Display adapter for terminal output.

    Listens to progress events and displays them in the terminal
    using simple text output with progress indicators.
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize terminal display.

        Args:
            verbose: Whether to show detailed progress
        """
        self._verbose = verbose
        self._logger = logging.getLogger(__name__)
        self._current_line = ""
        self._last_update = datetime.now()

    def handle_event(self, event: Any) -> None:
        """
        Handle a progress event.

        Args:
            event: Progress event to handle
        """
        try:
            event_type = type(event).__name__

            if hasattr(self, f'_handle_{event_type.lower().replace("event", "")}'):
                handler = getattr(
                    self, f'_handle_{event_type.lower().replace("event", "")}'
                )
                handler(event)
            else:
                self._handle_generic(event)

        except Exception as e:
            self._logger.error(f"Error handling event {type(event).__name__}: {e}")

    def _handle_taskstarted(self, event) -> None:
        """Handle TaskStartedEvent."""
        self._print_line(f"\n🚀 Starting {event.task_type.value}: {event.task_name}")
        if self._verbose and event.metadata:
            for key, value in event.metadata.items():
                self._print_line(f"   {key}: {value}")

    def _handle_taskfinished(self, event) -> None:
        """Handle TaskFinishedEvent."""
        if event.success:
            self._print_line(f"✅ Task completed: {event.task_name}")
        else:
            self._print_line(f"❌ Task failed: {event.task_name}")
            if event.error_message:
                self._print_line(f"   Error: {event.error_message}")

    def _handle_executionstarted(self, event) -> None:
        """Handle ExecutionStartedEvent."""
        self._print_line(
            f"\n📊 Task: {event.execution_name} ({event.total_items} items)"
        )

    def _handle_executionprogress(self, event) -> None:
        """Handle ExecutionProgressEvent."""
        # Show progress bar for execution
        progress_bar = self._create_progress_bar(event.progress_percent)
        status = f"[{progress_bar}] {event.progress_percent:5.1f}% - {event.item_name}"
        if event.message:
            status += f" - {event.message}"

        self._update_line(status)

    def _handle_executionfinished(self, event) -> None:
        """Handle ExecutionFinishedEvent."""
        if event.success:
            self._print_line(
                f"\n✅ Execution completed: {event.execution_name} ({event.total_processed} items)"
            )
        else:
            self._print_line(f"\n❌ Execution failed: {event.execution_name}")
            if event.error_message:
                self._print_line(f"   Error: {event.error_message}")

    def _handle_algorithmprogress(self, event) -> None:
        """Handle AlgorithmProgressEvent."""
        if self._verbose:
            progress_bar = self._create_progress_bar(event.progress_percent, width=20)
            status = f"  {event.algorithm_name}: [{progress_bar}] {event.progress_percent:5.1f}%"
            if event.message:
                status += f" - {event.message}"
            self._update_line(status)

    def _handle_error(self, event) -> None:
        """Handle ErrorEvent."""
        self._print_line(f"\n❌ Error ({event.error_type}): {event.error_message}")

    def _handle_generic(self, event) -> None:
        """Handle generic events."""
        if self._verbose:
            self._print_line(
                f"📝 {type(event).__name__}: {getattr(event, 'message', 'No message')}"
            )

    def _create_progress_bar(self, percentage: float, width: int = 30) -> str:
        """
        Create a text progress bar.

        Args:
            percentage: Progress percentage (0-100)
            width: Width of the progress bar

        Returns:
            Progress bar string
        """
        filled = int(width * percentage / 100)
        bar = "█" * filled + "░" * (width - filled)
        return bar

    def _print_line(self, text: str) -> None:
        """
        Print a line to the terminal.

        Args:
            text: Text to print
        """
        # Clear current line if needed
        if self._current_line:
            sys.stdout.write("\r" + " " * len(self._current_line) + "\r")

        print(text)
        self._current_line = ""
        sys.stdout.flush()

    def _update_line(self, text: str) -> None:
        """
        Update the current line (for progress indicators).

        Args:
            text: Text to display
        """
        # Limit update frequency to avoid terminal flooding
        now = datetime.now()
        if (now - self._last_update).total_seconds() < 0.1:  # 100ms minimum
            return

        # Clear previous line
        if self._current_line:
            sys.stdout.write("\r" + " " * len(self._current_line) + "\r")

        # Write new line
        sys.stdout.write(text)
        sys.stdout.flush()
        self._current_line = text
        self._last_update = now

    def finish(self) -> None:
        """Finish display and clean up."""
        if self._current_line:
            sys.stdout.write("\n")
            sys.stdout.flush()
            self._current_line = ""
