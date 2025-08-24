"""
CLI Progress Monitor using curses for real-time algorithm execution monitoring.
"""

import curses
import time
import signal
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional

from src.infrastructure.persistence.work_state.queries import (
    WorkStateQueries, 
    ProgressSummary, 
    ExecutionDetail, 
    ErrorSummary
)


class ProgressMonitor:
    """Real-time progress monitor using curses, htop-style."""
    
    def __init__(self, work_id: str):
        """Initialize monitor."""
        from src.application.services.work_service import get_work_service
        
        self.work_id = work_id
        self.work_service = get_work_service()
        self.running = True
        self.refresh_interval = 0.1  # seconds
        
        work_details = self.work_service.get(work_id)
        if not work_details or not work_details.get("output_path"):
            raise ValueError(f"Work {work_id} not found or has no output path")
            
        state_db_path = Path(work_details["output_path"]) / "state.db"
        if not state_db_path.exists():
            raise ValueError(f"State database not found at {state_db_path}")
            
        self.queries = WorkStateQueries(str(state_db_path))
        self.start_time = time.time()
        self.scroll_pos = 0
        self.executions: List[ExecutionDetail] = []
        self.logs: List[dict] = []

    def signal_handler(self, signum, frame):
        """Handle Ctrl+C gracefully."""
        self.running = False

    def format_progress_bar(self, progress: float, width: int) -> str:
        """Creates a text-based progress bar."""
        filled_len = int(width * progress)
        bar = '█' * filled_len + '─' * (width - filled_len)
        return f"[{bar}] {progress:.1%}"

    def format_execution_time(self, elapsed_seconds: float) -> str:
        """Format execution time as mm:ss.mmm"""
        total_ms = int(elapsed_seconds * 1000)
        minutes = total_ms // 60000
        seconds = (total_ms % 60000) // 1000
        milliseconds = total_ms % 1000
        
        if minutes > 0:
            return f"{minutes:02d}:{seconds:02d}.{milliseconds:03d}"
        else:
            return f"{seconds:02d}.{milliseconds:03d}s"

    def draw_header(self, stdscr, progress: ProgressSummary, work_status: str):
        """Draws the header with summary information."""
        h, w = stdscr.getmaxyx()
        if h < 5: return 0

        # Line 1: Work Info
        uptime = timedelta(seconds=int(time.time() - self.start_time))
        work_info = f"Work: {self.work_id}  Status: {work_status.upper()}  Uptime: {uptime}"
        stdscr.addstr(0, 1, work_info.ljust(w - 2), curses.color_pair(4) | curses.A_BOLD)

        # Line 2: Global Progress
        glob_prog = progress.global_progress
        glob_label = f"Overall Progress: {progress.global_execution['Finished']}/{progress.global_execution['Total']}"
        
        # Calculate available width for progress bar
        available_width = w - len(glob_label) - 10  # Reserve space for label + padding
        bar_width = max(10, min(40, available_width))  # Between 10-40 chars, but fit screen
        
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
        
        task_running_count = 1 if tasks.get('Running') else 0
        dset_running_count = 1 if dsets.get('Running') else 0
        cfg_running_count = 1 if cfgs.get('Running') else 0
        alg_running_count = 1 if algs.get('Running') else 0

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
            current_seq = exec_info['Finished'] + exec_info['Running']
            total_seq = exec_info['Total']
            combo_desc = f"Running combination: {combo['task_id']} > {combo['dataset_id']} > {combo['preset_id']}. Sequences: {current_seq}/{total_seq}"
            stdscr.addstr(4, 1, combo_desc.ljust(w-2)[:w-2], curses.A_BOLD)

        # Line 6: Execution Table Header
        header = f"{'SEQ':<5} {'ALGORITHM':<20} {'PROGRESS':<25} {'TIME':<10} {'DETAILS'}"
        stdscr.addstr(5, 0, header.ljust(w), curses.color_pair(5) | curses.A_REVERSE)
        
        return 6

    def draw_executions_list(self, stdscr, y_pos: int):
        """Draws the list of running executions."""
        h, w = stdscr.getmaxyx()
        max_rows = h - y_pos - 1

        for i in range(max_rows):
            exec_idx = self.scroll_pos + i
            if exec_idx >= len(self.executions):
                break
            
            exec_detail = self.executions[exec_idx]
            
            seq = exec_detail.sequencia
            algo = exec_detail.algorithm_id[:19]
            
            progress_bar = self.format_progress_bar(exec_detail.progress, 15)
            
            if exec_detail.finished_at and exec_detail.started_at:
                elapsed_time = exec_detail.finished_at - exec_detail.started_at
            elif exec_detail.started_at:
                elapsed_time = time.time() - exec_detail.started_at
            else:
                elapsed_time = 0
            time_str = self.format_execution_time(elapsed_time)
            
            details = (exec_detail.progress_message or '')[:w - 65]

            line = f"{seq:<5} {algo:<20} {progress_bar:<25} {time_str:<10} {details}"
            
            color = curses.color_pair(3)
            if exec_detail.status == 'running':
                color = curses.color_pair(3) if exec_detail.progress > 0 else curses.color_pair(2)
            elif exec_detail.status in ('failed', 'error'):
                color = curses.color_pair(1)
            elif exec_detail.status == 'completed':
                color = curses.color_pair(4) # Cyan for completed

            stdscr.addstr(y_pos + i, 0, line.ljust(w)[:w], color)

    def draw_footer(self, stdscr):
        """Draws the footer with controls and status."""
        h, w = stdscr.getmaxyx()
        if h <= 1: return
        
        controls = "Q:Quit C:Cancel"
        timestamp = datetime.now().strftime("%H:%M:%S")
        
        footer_text = f"{controls}".ljust(w - len(timestamp) - 1) + timestamp
        stdscr.addstr(h - 1, 0, footer_text, curses.color_pair(5) | curses.A_REVERSE)

    def draw_log_panel(self, stdscr, y_pos: int, max_height: int):
        """Draws a small panel with recent error and warning logs."""
        h, w = stdscr.getmaxyx()
        if max_height <= 1 or not self.logs: 
            return

        # Small panel header
        title = " Recent Events "
        stdscr.addstr(y_pos, 0, "┌" + "─" * (len(title)) + "┐" + "─" * (w - len(title) - 3), curses.color_pair(5))
        stdscr.addstr(y_pos, 1, title, curses.color_pair(5) | curses.A_REVERSE)
        y_pos += 1

        # Show only the most recent 2 entries
        for i, log_entry in enumerate(self.logs[:min(max_height - 1, 2)]):
            if y_pos >= h - 1:
                break
                
            log_time = datetime.fromtimestamp(log_entry['timestamp']).strftime('%H:%M:%S')
            
            if log_entry['type'] == 'error':
                color = curses.color_pair(1)
                icon = "✗"
            else:
                color = curses.color_pair(2)
                icon = "!"
            
            # Truncate message to fit width
            max_msg_len = w - 15  # Reserve space for time and icon
            message = log_entry['message']
            if len(message) > max_msg_len:
                message = message[:max_msg_len-3] + "..."
            
            line = f"│ {icon} [{log_time}] {message}"
            stdscr.addstr(y_pos + i, 0, line.ljust(w)[:w], color)

    def run(self, stdscr):
        """Main monitoring loop."""
        curses.curs_set(0)
        stdscr.nodelay(1)
        stdscr.timeout(50) # 50ms timeout for getch()
        
        curses.start_color()
        curses.init_pair(1, curses.COLOR_RED, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)
        curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)
        curses.init_pair(4, curses.COLOR_CYAN, curses.COLOR_BLACK)
        curses.init_pair(5, curses.COLOR_WHITE, curses.COLOR_BLUE)
        
        last_refresh = 0
        
        while self.running:
            h, w = stdscr.getmaxyx()
            current_time = time.time()
            
            key = stdscr.getch()
            if key in [ord('q'), ord('Q'), 27]:
                break
            elif key in [ord('c'), ord('C')]:
                self.handle_cancel(stdscr)
            elif key == curses.KEY_UP:
                self.scroll_pos = max(0, self.scroll_pos - 1)
            elif key == curses.KEY_DOWN:
                # Adjust scroll calculation for new layout
                log_panel_height = 3 if self.logs else 0  # Fixed small size
                max_scroll = max(0, len(self.executions) - (h - 6 - log_panel_height - 1))
                self.scroll_pos = min(max_scroll, self.scroll_pos + 1)

            if current_time - last_refresh >= self.refresh_interval:
                try:
                    stdscr.clear()
                    
                    progress = self.queries.get_work_progress_summary(self.work_id)
                    if not progress:
                        stdscr.addstr(0, 0, f"Work '{self.work_id}' not found!", curses.color_pair(1))
                        stdscr.refresh()
                        time.sleep(1)
                        continue
                    
                    work_status = self.queries.get_work_status(self.work_id) or 'unknown'
                    if work_status not in ['running', 'queued']:
                        self.show_final_screen(stdscr, work_status, progress)
                        break
                    
                    # Fetch executions for the current combination
                    current_combo_id = progress.current_combination_details.get('combination_id') if progress.current_combination_details else None
                    if current_combo_id:
                        old_executions_count = len(self.executions)
                        self.executions = self.queries.get_combination_executions_detail(current_combo_id)
                        new_executions_count = len(self.executions)
                    else:
                        old_executions_count = 0
                        new_executions_count = 0
                        self.executions = []

                    # Fetch and combine logs
                    errors = self.queries.get_error_summary(self.work_id, limit=2)
                    warnings = self.queries.get_execution_warnings(self.work_id, limit=2)
                    
                    logs = []
                    for e in errors:
                        logs.append({'type': 'error', 'timestamp': e.timestamp, 'message': f"Unit {e.unit_id[:8]}: {e.error_message}"})
                    for w in warnings:
                        logs.append({'type': 'warning', 'timestamp': w['timestamp'], 'message': f"Unit {(w.get('unit_id') or 'N/A')[:8]}: {w.get('message', 'Unknown warning')}"})
                    self.logs = sorted(logs, key=lambda x: x.get('timestamp', 0), reverse=True)[:2]  # Keep only 2 most recent

                    # Calculate layout - fixed small log panel
                    log_panel_height = 3 if self.logs else 0 # 1 header + max 2 log lines
                    exec_list_height = h - 6 - log_panel_height - 1  # header - logs - footer
                    
                    # Smart auto-scroll logic
                    if exec_list_height > 0 and new_executions_count > 0:
                        max_scroll = max(0, new_executions_count - exec_list_height)
                        
                        # Check if we were showing the last line before update
                        was_at_bottom = (old_executions_count <= exec_list_height) or \
                                       (self.scroll_pos >= max(0, old_executions_count - exec_list_height))
                        
                        # If new executions were added and we were at bottom, auto-scroll
                        if new_executions_count > old_executions_count and was_at_bottom:
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

    def handle_cancel(self, stdscr):
        """Shows a confirmation dialog for cancelling the work."""
        h, w = stdscr.getmaxyx()
        msg = f"Cancel work '{self.work_id}'? [Y/N]"
        stdscr.addstr(h // 2, (w - len(msg)) // 2, msg, curses.color_pair(1) | curses.A_BOLD)
        stdscr.refresh()
        
        confirm_key = -1
        while confirm_key not in [ord('y'), ord('Y'), ord('n'), ord('N')]:
            confirm_key = stdscr.getch()
            time.sleep(0.1)

        if confirm_key in [ord('y'), ord('Y')]:
            try:
                self.work_service.cancel(self.work_id)
            except Exception as e:
                pass # Ignore errors, will be reflected in status
    
    def show_final_screen(self, stdscr, status: str, progress: ProgressSummary):
        """Displays a final summary screen when the work is done."""
        stdscr.clear()
        h, w = stdscr.getmaxyx()
        
        title = f"Work Finished with status: {status.upper()}"
        finished = progress.global_execution["Finished"]
        total = progress.global_execution["Total"]
        final_progress = f"Final Progress: {finished}/{total} executions completed."
        
        stdscr.addstr(h//2 - 1, (w - len(title)) // 2, title, curses.color_pair(3) | curses.A_BOLD)
        stdscr.addstr(h//2, (w - len(final_progress)) // 2, final_progress)
        stdscr.addstr(h//2 + 2, (w - 30) // 2, "Monitor will exit in 3 seconds...")
        stdscr.refresh()
        
        for _ in range(30):
            if stdscr.getch() != -1: break
            time.sleep(0.1)

    def start(self):
        """Start the monitor."""
        signal.signal(signal.SIGINT, self.signal_handler)
        try:
            if not self.queries.work_exists(self.work_id):
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
            self.queries.close()


def main():
    """CLI entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description="CSPBench Progress Monitor (htop-style)")
    parser.add_argument("work_id", help="Work ID to monitor")
    
    args = parser.parse_args()
    
    try:
        monitor = ProgressMonitor(args.work_id)
        return monitor.start()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())