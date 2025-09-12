"""
CLI Progress Monitor using curses for real-time algorithm execution monitoring.

This module provides a terminal-based progress monitor that displays real-time
execution status of CSPBench work items in an htop-style interface.
"""

import argparse
import curses
import signal
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional, Dict, Any

from src.infrastructure.persistence.work_state.core import (
    WorkPersistence,
    ExecutionDetail,
    ProgressSummary,
    ErrorSummary,
)
from src.infrastructure.logging_config import get_logger

pm_logger = get_logger("CSPBench.CLI.ProgressMonitor")


class ProgressMonitor:
    """Real-time progress monitor using curses, htop-style interface.
    
    This class provides a terminal-based interface for monitoring the progress
    of CSPBench work executions. It displays real-time information about:
    - Overall work progress
    - Individual algorithm executions
    - System status and logs
    - Resource usage
    
    The interface supports keyboard navigation and can pause/resume work execution.
    
    Attributes:
        work_id: The identifier of the work item being monitored
        running: Flag indicating if the monitor is active
        refresh_interval: Screen refresh interval in seconds
        start_time: Monitor start timestamp
        scroll_pos: Current scroll position in execution list
        executions: List of execution details
        logs: List of recent log entries
        previous_status: Previous work status for change detection
        persistence: Work persistence interface
    """

    def __init__(self, work_id: str):
        """Initialize the progress monitor.
        
        Args:
            work_id: The identifier of the work item to monitor
            
        Raises:
            ValueError: If the work item is not found or fails to initialize
        """
        print(f"🚀 Initializing monitor for work: {work_id}")

        self.work_id = work_id
        self.running = True
        self.refresh_interval = 0.1  # seconds
        self.start_time = time.time()
        self.scroll_pos = 0
        self.executions: List[ExecutionDetail] = []
        self.logs: List[dict] = []
        self.previous_status = None

        # Create WorkPersistence with default configuration
        print(f"⚙️ Connecting to persistence system...", end="\r")
        self.persistence = WorkPersistence()

        # Wait for work to exist and be in execution
        print(f"⏳ Waiting for work to be available and running...", end="\r")
        if not self._wait_for_work_ready(timeout=60):
            raise ValueError(f"Timeout: Work '{work_id}' not found or failed to start execution")

        print(" " * 60, end="\r")
        print("✅ Monitor initialized successfully!")
        print("🎯 Starting interface...")

    def _wait_for_work_ready(self, timeout: int = 60) -> bool:
        """Wait for work to exist and be available for monitoring.

        Args:
            timeout: Timeout in seconds to wait for work readiness

        Returns:
            True if work is ready for monitoring, False if timeout occurred
        """
        start_time = time.time()
        last_log = 0.0

        while time.time() - start_time < timeout:
            current_time = time.time()
            
            # Periodic logging
            if current_time - last_log >= 5:
                elapsed = int(current_time - start_time)
                print(f"   Waiting for work ({elapsed}s)...", end="\r")
                last_log = current_time
            
            try:
                # Check if work exists
                work = self.persistence.work_get(self.work_id)
                if not work:
                    time.sleep(1)
                    continue
                
                status = work.get('status')
                print(f"   Work found with status: {status}", end="\r")
                
                # For final statuses, allow monitoring
                if status in ["completed", "failed", "canceled", "error"]:
                    print(" " * 60, end="\r")
                    print(f"   Work in final state: {status}")
                    self.previous_status = status
                    return True
                
                # For execution statuses, check if there's progress data
                if status in ["running", "queued"]:
                    try:
                        progress_data = self.persistence.get_work_progress_summary(self.work_id)
                        if progress_data and progress_data.global_execution.get("Total", 0) > 0:
                            print(f"   Work ready: status={status}, total_sequences={progress_data.global_execution.get('Total', 0)}", end="\r")
                            self.previous_status = status
                            return True
                    except Exception as e:
                        pm_logger.error("[PM] Error getting work progress: %s", e)
                        pass  # Still no progress data
                        
            except Exception as e:
                pm_logger.debug("[PM] Error checking work: %s", e)
            
            time.sleep(1)
        
        return False

    def signal_handler(self, signum, frame):
        """Handle Ctrl+C gracefully by pausing the work.
        
        Args:
            signum: Signal number
            frame: Current stack frame
        """
        try:
            # Import and use WorkService to pause the work
            from src.application.services.work_service import get_work_service
            
            work_service = get_work_service()
            success = work_service.pause(self.work_id)
            
            if success:
                pm_logger.info(f"Work {self.work_id} paused via Ctrl+C")
                print(f"\n✅ Work {self.work_id} paused successfully!")
                print("🔄 Use 'cspbench monitor {work_id}' to continue monitoring")
            else:
                pm_logger.warning(f"Failed to pause work {self.work_id} via Ctrl+C")
                print(f"\n⚠️  Could not pause work {self.work_id}")
                
        except Exception as e:
            pm_logger.error(f"Error pausing work {self.work_id} via Ctrl+C: {e}")
            print(f"\n❌ Error pausing work: {e}")
        
        self.running = False

    def format_progress_bar(
        self, progress: float, width: int, show_percent: bool = True
    ) -> str:
        """Create a text-based progress bar.
        
        Args:
            progress: Progress value between 0.0 and 1.0
            width: Width of the progress bar in characters
            show_percent: Whether to show percentage alongside the bar
            
        Returns:
            Formatted progress bar string
        """
        filled_len = int(width * progress)
        bar = "█" * filled_len + "─" * (width - filled_len)
        if show_percent:
            return f"[{bar}] {progress:.1%}"
        else:
            return bar

    def format_execution_time(self, elapsed_seconds: float) -> str:
        """Format execution time as mm:ss.mmm.
        
        Args:
            elapsed_seconds: Elapsed time in seconds
            
        Returns:
            Formatted time string
        """
        total_ms = int(elapsed_seconds * 1000)
        minutes = total_ms // 60000
        seconds = (total_ms % 60000) // 1000
        milliseconds = total_ms % 1000

        if minutes > 0:
            return f"{minutes:02d}:{seconds:02d}.{milliseconds:03d}"
        else:
            return f"{seconds:02d}.{milliseconds:03d}s"

    def draw_header(self, stdscr, progress: ProgressSummary, work_status: str):
        """Draw the header with summary information.
        
        Args:
            stdscr: Curses screen object
            progress: Progress summary data
            work_status: Current work status
            
        Returns:
            Number of lines used by the header
        """
        h, w = stdscr.getmaxyx()
        if h < 5:
            return 0

        # Line 1: Work Info
        uptime = timedelta(seconds=int(time.time() - self.start_time))
        work_info = (
            f"Work: {self.work_id}  Status: {work_status.upper()}  Uptime: {uptime}"
        )
        stdscr.addstr(
            0, 1, work_info.ljust(w - 2), curses.color_pair(4) | curses.A_BOLD
        )

        # Line 2: Global Progress
        glob_prog = progress.global_progress
        glob_label = f"Overall Progress: {progress.global_execution['Finished']}/{progress.global_execution['Total']}"

        # Calculate available width for progress bar
        available_width = w - len(glob_label) - 10  # Reserve space for label + padding
        bar_width = max(
            10, min(40, available_width)
        )  # Between 10-40 chars, but fit screen

        glob_bar = self.format_progress_bar(glob_prog, bar_width)

        stdscr.addstr(1, 1, glob_label, curses.A_BOLD)
        if len(glob_label) + len(glob_bar) + 4 <= w - 2:
            stdscr.addstr(1, len(glob_label) + 2, glob_bar, curses.color_pair(3))
        else:
            # If still too long, show on next line
            stdscr.addstr(2, 3, glob_bar, curses.color_pair(3))
            # Adjust return value to account for extra line
            return 7

        # Line 3-4: Hierarchical Progress in aligned columns
        tasks = progress.tasks
        dsets = progress.datasets
        cfgs = progress.configs
        algs = progress.algorithms

        task_running_count = 1 if tasks.get("Running") else 0
        dset_running_count = 1 if dsets.get("Running") else 0
        cfg_running_count = 1 if cfgs.get("Running") else 0
        alg_running_count = 1 if algs.get("Running") else 0

        # Calculate column widths for alignment
        col1_width = w // 2 - 2  # Half screen minus padding
        col2_width = w - col1_width - 4  # Remaining width

        # Line 3: Task and Dataset
        task_str = f"Task: {tasks.get('Running', '-') or '-'} ({tasks['Finished'] + task_running_count}/{tasks['Total']})"
        dset_str = f"Dataset: {dsets.get('Running', '-') or '-'} ({dsets['Finished'] + dset_running_count}/{dsets['Total']})"

        stdscr.addstr(2, 1, task_str[:col1_width].ljust(col1_width))
        stdscr.addstr(2, col1_width + 2, dset_str[:col2_width])

        # Line 4: Config and Algorithm
        cfg_str = f"Config: {cfgs.get('Running', '-') or '-'} ({cfgs['Finished'] + cfg_running_count}/{cfgs['Total']})"
        alg_str = f"Algorithm: {algs.get('Running', '-') or '-'} ({algs['Finished'] + alg_running_count}/{algs['Total']})"

        stdscr.addstr(3, 1, cfg_str[:col1_width].ljust(col1_width))
        stdscr.addstr(3, col1_width + 2, alg_str[:col2_width])

        # Line 5: Current Combination Description
        combo = progress.current_combination_details
        if combo:
            exec_info = progress.execution
            current_seq = exec_info["Finished"] + exec_info["Running"]
            total_seq = exec_info["Total"]
            combo_desc = f"Running combination: {combo['task_id']} > {combo['dataset_id']} > {combo['preset_id']}. Sequences: {current_seq}/{total_seq}"
            stdscr.addstr(4, 1, combo_desc.ljust(w - 2)[: w - 2], curses.A_BOLD)

        # Line 6: Execution Table Header
        header = (
            f"{'ST SEQ':<6} {'ALGORITHM':<19} {'PROGRESS':<17} {'TIME':<10} {'DETAILS'}"
        )
        stdscr.addstr(5, 0, header.ljust(w), curses.color_pair(5) | curses.A_REVERSE)

        return 6

    def draw_executions_list(self, stdscr, y_pos: int):
        """Draw the list of executions with visual status indicators.
        
        Args:
            stdscr: Curses screen object
            y_pos: Starting Y position for the execution list
        """
        h, w = stdscr.getmaxyx()
        max_rows = h - y_pos - 1

        for i in range(max_rows):
            exec_idx = self.scroll_pos + i
            if exec_idx >= len(self.executions):
                break

            exec_detail = self.executions[exec_idx]

            seq = exec_detail.sequencia
            algo = exec_detail.algorithm_id[:19]

            # Visual status indicator
            if exec_detail.status == "completed":
                status_icon = "✓"
                color = curses.color_pair(4)  # Cyan for completed
                progress_bar = f"[{'█' * 15}]"  # Full bar for completed
            elif exec_detail.status == "running":
                status_icon = "▶"
                color = (
                    curses.color_pair(3)
                    if exec_detail.progress > 0
                    else curses.color_pair(2)
                )  # Green/Yellow
                progress_bar = f"[{self.format_progress_bar(exec_detail.progress, 15, show_percent=False)}]"
            elif exec_detail.status in ("failed", "error"):
                status_icon = "✗"
                color = curses.color_pair(1)  # Red for failed
                progress_bar = f"[{'▓' * 15}]"  # Error pattern
            elif exec_detail.status == "queued":
                status_icon = "⋯"
                color = (
                    curses.color_pair(6)
                    if hasattr(curses, "COLOR_MAGENTA")
                    else curses.color_pair(2)
                )  # Magenta or fallback
                progress_bar = f"[{'─' * 15}]"  # Empty bar for queued
            else:
                status_icon = "?"
                color = curses.color_pair(2)  # Yellow for unknown
                progress_bar = f"[{self.format_progress_bar(exec_detail.progress, 15, show_percent=False)}]"

            # Calculate elapsed time only if execution has started
            if exec_detail.finished_at and exec_detail.started_at:
                elapsed_time = exec_detail.finished_at - exec_detail.started_at
                time_str = self.format_execution_time(elapsed_time)
            elif exec_detail.started_at:
                elapsed_time = time.time() - exec_detail.started_at
                time_str = self.format_execution_time(elapsed_time)
            else:
                # No time display for queued executions
                time_str = "-"

            details = (exec_detail.progress_message or "")[: w - 70]

            line = f"{status_icon} {seq:<4} {algo:<19} {progress_bar:<17} {time_str:<10} {details}"

            stdscr.addstr(y_pos + i, 0, line.ljust(w)[:w], color)

    def draw_footer(self, stdscr):
        """Draw the footer with controls, legend and status.
        
        Args:
            stdscr: Curses screen object
        """
        h, w = stdscr.getmaxyx()
        if h <= 1:
            return

        controls = "Q:Quit P:Pause Ctrl+C:Pause ↑↓:Scroll"
        legend = "✓:Done ▶:Running ✗:Failed ⋯:Queued"
        timestamp = datetime.now().strftime("%H:%M:%S")

        # If screen is wide enough, show legend and controls on same line
        if len(controls) + len(legend) + len(timestamp) + 6 <= w:
            footer_text = (
                f"{controls}  {legend}".ljust(w - len(timestamp) - 1) + timestamp
            )
        else:
            # Otherwise just show controls and timestamp
            footer_text = f"{controls}".ljust(w - len(timestamp) - 1) + timestamp

        stdscr.addstr(h - 1, 0, footer_text, curses.color_pair(5) | curses.A_REVERSE)

    def draw_log_panel(self, stdscr, y_pos: int, max_height: int):
        """Draw a small panel with recent error and warning logs.
        
        Args:
            stdscr: Curses screen object
            y_pos: Starting Y position for the log panel
            max_height: Maximum height available for the log panel
        """
        h, w = stdscr.getmaxyx()
        if max_height <= 1 or not self.logs:
            return

        # Small panel header
        title = " Recent Events "
        stdscr.addstr(
            y_pos,
            0,
            "┌" + "─" * (len(title)) + "┐" + "─" * (w - len(title) - 3),
            curses.color_pair(5),
        )
        stdscr.addstr(y_pos, 1, title, curses.color_pair(5) | curses.A_REVERSE)
        y_pos += 1

        # Show only the most recent 2 entries
        for i, log_entry in enumerate(self.logs[: min(max_height - 1, 2)]):
            if y_pos >= h - 1:
                break

            log_time = datetime.fromtimestamp(log_entry["timestamp"]).strftime(
                "%H:%M:%S"
            )

            # Choose color and icon based on event type
            if log_entry["type"] == "error":
                color = curses.color_pair(1)  # Red
                icon = "✗"
            elif log_entry["type"] == "warning":
                color = curses.color_pair(2)  # Yellow
                icon = "⚠"
            else:  # info, progress, or other
                color = curses.color_pair(4)  # Cyan
                icon = "ℹ"

            # Truncate message to fit width
            max_msg_len = w - 15  # Reserve space for time and icon
            message = log_entry["message"]
            if len(message) > max_msg_len:
                message = message[: max_msg_len - 3] + "..."

            line = f"│ {icon} [{log_time}] {message}"
            stdscr.addstr(y_pos + i, 0, line.ljust(w)[:w], color)

    def run(self, stdscr):
        """Main monitoring loop.
        
        Args:
            stdscr: Curses screen object
        """
        curses.curs_set(0)
        stdscr.nodelay(1)
        stdscr.timeout(50)  # 50ms timeout for getch()

        curses.start_color()
        curses.init_pair(1, curses.COLOR_RED, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)
        curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)
        curses.init_pair(4, curses.COLOR_CYAN, curses.COLOR_BLACK)
        curses.init_pair(5, curses.COLOR_WHITE, curses.COLOR_BLUE)
        curses.init_pair(6, curses.COLOR_MAGENTA, curses.COLOR_BLACK)

        last_refresh = 0

        while self.running:
            h, w = stdscr.getmaxyx()
            current_time = time.time()

            key = stdscr.getch()
            if key in [ord("q"), ord("Q"), 27]:
                break
            elif key in [ord("p"), ord("P")]:
                self.handle_pause(stdscr)
            elif key == curses.KEY_UP:
                self.scroll_pos = max(0, self.scroll_pos - 1)
            elif key == curses.KEY_DOWN:
                # Adjust scroll calculation for new layout
                log_panel_height = 3 if self.logs else 0  # Fixed small size
                max_scroll = max(
                    0, len(self.executions) - (h - 6 - log_panel_height - 1)
                )
                self.scroll_pos = min(max_scroll, self.scroll_pos + 1)

            # Check if was paused or canceled and should exit
            if not self.running:
                break

            if current_time - last_refresh >= self.refresh_interval:
                try:
                    stdscr.clear()

                    progress = self.persistence.get_work_progress_summary(self.work_id)
                    if not progress:
                        stdscr.addstr(
                            0,
                            0,
                            f"Work '{self.work_id}' not found!",
                            curses.color_pair(1),
                        )
                        stdscr.refresh()
                        time.sleep(1)
                        continue

                    work = self.persistence.work_get(self.work_id)
                    work_status = work.get('status') if work else "unknown"

                    # Detect status change from 'queued' to something other than 'running'
                    if (
                        self.previous_status == "queued"
                        and work_status not in ["queued", "running"]
                        and work_status != "unknown"
                    ):
                        self.show_status_change_screen(
                            stdscr, self.previous_status, work_status, progress
                        )
                        break

                    # Update previous status
                    if work_status != "unknown":
                        self.previous_status = work_status

                    if work_status == "paused":
                        self.show_paused_screen(stdscr, work_status, progress)
                        # Exit monitor when paused
                        break
                    elif work_status not in ["running", "queued"]:
                        self.show_final_screen(stdscr, work_status, progress)
                        break

                    # Fetch executions for the current combination
                    current_combo_id = (
                        progress.current_combination_details.get("combination_id")
                        if progress.current_combination_details
                        else None
                    )
                    if current_combo_id:
                        old_executions_count = len(self.executions)
                        # Get ALL executions for the combination, not just running ones
                        self.executions = (
                            self.persistence.get_combination_executions_detail(
                                current_combo_id
                            )
                        )
                        # Sort by sequence number for better visualization
                        self.executions.sort(key=lambda x: x.sequencia)
                        new_executions_count = len(self.executions)
                    else:
                        old_executions_count = 0
                        new_executions_count = 0
                        self.executions = []

                    # Fetch only events from events table (not execution errors)
                    # Get recent events of relevant types (warnings, errors, info, progress)
                    recent_events = self.persistence.get_events(
                        self.work_id, limit=4, event_types=['warning', 'error', 'info', 'progress']
                    )

                    logs = []
                    # Only show actual events from events table, not algorithm execution errors
                    for event in recent_events:
                        event_type = event.get('event_type', 'unknown')
                        message = event.get('message', 'No message')
                        unit_id = event.get('entity_data', {}).get('unit_id', 'N/A')
                        
                        # Use appropriate icon based on event type
                        if event_type == 'error':
                            log_type = 'error'
                        elif event_type == 'warning':
                            log_type = 'warning'
                        else:
                            log_type = 'info'  # For progress, info, and other types
                        
                        logs.append(
                            {
                                "type": log_type,
                                "timestamp": event.get('timestamp', 0),
                                "message": f"Unit {unit_id[:8]}: {message}",
                            }
                        )
                    
                    self.logs = sorted(
                        logs, key=lambda x: x.get("timestamp", 0), reverse=True
                    )[
                        :2
                    ]  # Keep only 2 most recent

                    # Calculate layout - fixed small log panel
                    log_panel_height = (
                        3 if self.logs else 0
                    )  # 1 header + max 2 log lines
                    exec_list_height = (
                        h - 6 - log_panel_height - 1
                    )  # header - logs - footer

                    # Smart auto-scroll logic - focus on running/next executions
                    if exec_list_height > 0 and new_executions_count > 0:
                        max_scroll = max(0, new_executions_count - exec_list_height)

                        # Find first running or queued execution
                        first_active_idx = None
                        for i, exec_detail in enumerate(self.executions):
                            if exec_detail.status in ("running", "queued"):
                                first_active_idx = i
                                break

                        # Check if we were showing the last line before update
                        was_at_bottom = (old_executions_count <= exec_list_height) or (
                            self.scroll_pos
                            >= max(0, old_executions_count - exec_list_height)
                        )

                        # If new executions were added and we were at bottom, or if no manual scrolling occurred
                        if (
                            new_executions_count > old_executions_count
                            and was_at_bottom
                        ):
                            # Try to center on first active execution, but don't go past bottom
                            if first_active_idx is not None:
                                target_scroll = max(
                                    0,
                                    min(
                                        max_scroll,
                                        first_active_idx - exec_list_height // 2,
                                    ),
                                )
                                self.scroll_pos = target_scroll
                            else:
                                self.scroll_pos = max_scroll
                        elif self.scroll_pos > max_scroll:
                            # Ensure scroll position is valid
                            self.scroll_pos = max_scroll

                    # Draw interface
                    y_pos = self.draw_header(stdscr, progress, work_status)

                    if exec_list_height > 0:
                        self.draw_executions_list(stdscr, y_pos)

                    # Draw small log panel just above footer
                    if self.logs:
                        log_y_pos = h - 1 - log_panel_height
                        self.draw_log_panel(stdscr, log_y_pos, log_panel_height)

                    self.draw_footer(stdscr)

                    stdscr.refresh()
                    last_refresh = current_time

                except Exception as e:
                    stdscr.clear()
                    stdscr.addstr(0, 0, f"Error: {str(e)}", curses.color_pair(1))
                    stdscr.addstr(1, 0, "Press 'q' to quit")
                    stdscr.refresh()
                    time.sleep(1)

            time.sleep(0.02)

    def handle_pause(self, stdscr):
        """Show a confirmation dialog for pausing the work.
        
        Args:
            stdscr: Curses screen object
        """
        h, w = stdscr.getmaxyx()

        # Create dialog box with background
        dialog_width = 50
        dialog_height = 7
        start_y = (h - dialog_height) // 2
        start_x = (w - dialog_width) // 2

        # Create dialog window
        dialog_win = curses.newwin(dialog_height, dialog_width, start_y, start_x)
        dialog_win.bkgd(" ", curses.color_pair(5))  # Blue background
        dialog_win.box()

        # Title
        title = "⏸️ PAUSE WORK"
        dialog_win.addstr(1, (dialog_width - len(title)) // 2, title, curses.A_BOLD)

        # Message
        msg1 = f"Pause work '{self.work_id[:20]}...'?"
        msg2 = "Press Y to confirm or N to cancel"
        dialog_win.addstr(3, (dialog_width - len(msg1)) // 2, msg1)
        dialog_win.addstr(4, (dialog_width - len(msg2)) // 2, msg2)

        # Options
        options = "[Y] Yes    [N] No"
        dialog_win.addstr(5, (dialog_width - len(options)) // 2, options, curses.A_BOLD)

        dialog_win.refresh()

        confirm_key = -1
        while confirm_key not in [
            ord("y"),
            ord("Y"),
            ord("n"),
            ord("N"),
            27,
        ]:  # 27 = ESC
            confirm_key = stdscr.getch()
            time.sleep(0.1)

        if confirm_key in [ord("y"), ord("Y")]:
            try:
                # Import and use WorkService to actually pause the work
                from src.application.services.work_service import get_work_service
                
                work_service = get_work_service()
                success = work_service.pause(self.work_id)
                
                if success:
                    self.show_paused_confirmation_screen(stdscr)
                    self.running = False
                else:
                    # Show error if couldn't pause
                    self.show_action_error_screen(stdscr, "pause")
                    
            except Exception as e:
                pm_logger.error(f"Error pausing work {self.work_id}: {e}")
                self.show_action_error_screen(stdscr, "pause")
                pass  # Continue monitoring in case of error

    def show_paused_confirmation_screen(self, stdscr):
        """Display a screen when the work is successfully paused.
        
        Args:
            stdscr: Curses screen object
        """
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = "⏸️ WORK PAUSED"
        msg1 = f"Work '{self.work_id}' was paused successfully."
        msg2 = "Monitor will exit in 3 seconds..."

        stdscr.addstr(
            h // 2 - 1,
            (w - len(title)) // 2,
            title,
            curses.color_pair(2) | curses.A_BOLD,
        )
        stdscr.addstr(h // 2, (w - len(msg1)) // 2, msg1)
        stdscr.addstr(h // 2 + 1, (w - len(msg2)) // 2, msg2)
        stdscr.refresh()

        for _ in range(30):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def show_paused_screen(self, stdscr, status: str, progress):
        """Display a paused screen when the work is paused.
        
        Args:
            stdscr: Curses screen object
            status: Current work status
            progress: Progress summary data
        """
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = f"⏸️ Work Paused: {status.upper()}"
        finished = progress.global_execution["Finished"]
        total = progress.global_execution["Total"]
        pause_progress = (
            f"Progress at pause: {finished}/{total} executions completed."
        )

        stdscr.addstr(
            h // 2 - 1,
            (w - len(title)) // 2,
            title,
            curses.color_pair(2) | curses.A_BOLD,
        )
        stdscr.addstr(h // 2, (w - len(pause_progress)) // 2, pause_progress)
        stdscr.addstr(h // 2 + 1, (w - 60) // 2, "Work was paused successfully.")
        stdscr.addstr(
            h // 2 + 2, (w - 40) // 2, "Monitor will exit in 3 seconds..."
        )
        stdscr.refresh()

        for _ in range(30):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def show_status_change_screen(
        self, stdscr, old_status: str, new_status: str, progress: ProgressSummary
    ):
        """Display a screen when work status changes from queued to non-running status.
        
        Args:
            stdscr: Curses screen object
            old_status: Previous work status
            new_status: New work status
            progress: Progress summary data
        """
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        # Determine color and message based on new status
        if new_status in ["failed", "error"]:
            color = curses.color_pair(1)  # Red
            title = f"Work Failed: {new_status.upper()}"
            message = "Work failed before starting execution."
        elif new_status == "canceled":
            color = curses.color_pair(2)  # Yellow
            title = f"Work Canceled: {new_status.upper()}"
            message = "Work was canceled before starting execution."
        elif new_status == "completed":
            color = curses.color_pair(3)  # Green
            title = f"Work Completed: {new_status.upper()}"
            message = "Work completed quickly."
        else:
            color = curses.color_pair(2)  # Yellow
            title = f"Status Changed: {old_status.upper()} → {new_status.upper()}"
            message = (
                f"Work changed from '{old_status}' to '{new_status}' without executing."
            )

        finished = progress.global_execution["Finished"]
        total = progress.global_execution["Total"]
        progress_info = f"Progress: {finished}/{total} executions completed."

        # Calculate centered positions
        title_y = h // 2 - 2
        message_y = h // 2 - 1
        progress_y = h // 2
        exit_y = h // 2 + 2

        stdscr.addstr(title_y, (w - len(title)) // 2, title, color | curses.A_BOLD)
        stdscr.addstr(message_y, (w - len(message)) // 2, message)
        stdscr.addstr(progress_y, (w - len(progress_info)) // 2, progress_info)
        stdscr.addstr(exit_y, (w - 40) // 2, "Monitor will exit in 5 seconds...")
        stdscr.refresh()

        # Wait 5 seconds or key press
        for _ in range(50):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def show_final_screen(self, stdscr, status: str, progress: ProgressSummary):
        """Display a final summary screen when the work is done.
        
        Args:
            stdscr: Curses screen object
            status: Final work status
            progress: Progress summary data
        """
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        # Choose emoji based on status
        if status == "completed":
            emoji = "🎉"
        elif status in ["failed", "error"]:
            emoji = "❌"
        elif status == "canceled":
            emoji = "🛑"
        else:
            emoji = "ℹ️"

        title = f"{emoji} Work Finished with status: {status.upper()}"
        finished = progress.global_execution["Finished"]
        total = progress.global_execution["Total"]
        final_progress = f"Final Progress: {finished}/{total} executions completed."

        stdscr.addstr(
            h // 2 - 1,
            (w - len(title)) // 2,
            title,
            curses.color_pair(3) | curses.A_BOLD,
        )
        stdscr.addstr(h // 2, (w - len(final_progress)) // 2, final_progress)
        stdscr.addstr(h // 2 + 2, (w - 30) // 2, "Monitor will exit in 3 seconds...")
        stdscr.refresh()

        for _ in range(30):
            if stdscr.getch() != -1:
                break
            time.sleep(0.1)

    def show_action_error_screen(self, stdscr, action: str):
        """Display an error screen when an action fails.
        
        Args:
            stdscr: Curses screen object
            action: The action that failed (e.g., "pause")
        """
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = f"❌ ERROR {action.upper()}"
        msg1 = f"Could not {action} work '{self.work_id}'"
        msg2 = "Check logs or try again"
        msg3 = "Press any key to continue..."

        stdscr.addstr(
            h // 2 - 2,
            (w - len(title)) // 2,
            title,
            curses.color_pair(1) | curses.A_BOLD,
        )
        stdscr.addstr(h // 2 - 1, (w - len(msg1)) // 2, msg1)
        stdscr.addstr(h // 2, (w - len(msg2)) // 2, msg2)
        stdscr.addstr(h // 2 + 2, (w - len(msg3)) // 2, msg3, curses.A_BOLD)
        stdscr.refresh()

        # Wait for key press
        stdscr.getch()

    def start(self):
        """Start the progress monitor.
        
        Returns:
            Exit code (0 for success, 1 for error)
        """
        signal.signal(signal.SIGINT, self.signal_handler)
        try:
            if not self.persistence.work_get(self.work_id):
                print(f"Error: Work '{self.work_id}' not found!")
                return 1
            curses.wrapper(self.run)
            return 0
        except KeyboardInterrupt:
            return 0
        except Exception as e:
            print(f"Error starting monitor: {e}")
            return 1
        finally:
            if hasattr(self, 'persistence'):
                self.persistence.close()


def main():
    """CLI entry point for the progress monitor."""

    parser = argparse.ArgumentParser(
        description="CSPBench Progress Monitor (htop-style)"
    )
    parser.add_argument("work_id", help="Work ID to monitor")

    args = parser.parse_args()

    try:
        monitor = ProgressMonitor(args.work_id)
        return monitor.start()
    except RuntimeError as e:
        if "Timeout" in str(e):
            if "WorkService" in str(e):
                print(
                    "\n❌ Error: WorkService system was not initialized in time.",
                    file=sys.stderr,
                )
            else:
                print(f"\n❌ Runtime error: {e}", file=sys.stderr)
        else:
            print(f"\n❌ Runtime error: {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        if "not found" in str(e):
            print(f"\n🔍 Error: Work '{args.work_id}' not found.", file=sys.stderr)
        elif "Timeout" in str(e):
            print(
                f"\n⏰ Error: Timeout waiting for work '{args.work_id}' initialization.",
                file=sys.stderr,
            )
        else:
            print(f"\n⚠️ Configuration error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\n💥 Unexpected error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
