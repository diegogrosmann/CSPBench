"""
CLI Progress Monitor using curses for real-time algorithm execution monitoring.
"""

import curses
import time
import signal
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from src.infrastructure.persistence.work_state.queries import (
    WorkStateQueries, 
    ProgressSummary, 
    ExecutionDetail, 
    ErrorSummary
)


class ProgressMonitor:
    """Real-time progress monitor using curses."""
    
    def __init__(self, work_id: str):
        """Initialize monitor.
        
        Args:
            work_id: Work ID to monitor
        """
        from src.application.services.work_service import get_work_service
        
        self.work_id = work_id
        self.work_service = get_work_service()
        self.running = True
        self.refresh_interval = 2  # seconds
        
        # Find state database path
        work_details = self.work_service.get(work_id)
        if not work_details or not work_details.get("output_path"):
            raise ValueError(f"Work {work_id} not found or has no output path")
            
        state_db_path = Path(work_details["output_path"]) / "state.db"
        if not state_db_path.exists():
            raise ValueError(f"State database not found at {state_db_path}")
            
        self.queries = WorkStateQueries(str(state_db_path))
        
    def signal_handler(self, signum, frame):
        """Handle Ctrl+C gracefully."""
        self.running = False
    
    def format_progress_bar(self, current: int, total: int, width: int = 25, show_percentage: bool = True) -> str:
        """Create a beautiful text progress bar with colors."""
        if total == 0:
            return "█" * width + " 0/0  0%"
        
        progress = current / total
        filled = int(width * progress)
        
        # Use different characters for better visual appeal
        complete_char = "█"
        partial_char = "▓"
        empty_char = "░"
        
        # Create smooth progress bar
        bar = complete_char * filled + empty_char * (width - filled)
        percentage = int(progress * 100)
        
        if show_percentage:
            return f"{bar} {current:,}/{total:,} ({percentage:3d}%)"
        else:
            return f"{bar} {current:,}/{total:,}"
    
    def draw_header(self, stdscr, progress: ProgressSummary, y_pos: int) -> int:
        """Draw an enhanced header with overall progress."""
        height, width = stdscr.getmaxyx()
        
        # Title with decorative border
        title = f"🚀 CSPBench Progress Monitor"
        work_info = f"Work ID: {self.work_id[:12]}{'...' if len(self.work_id) > 12 else ''}"
        
        stdscr.addstr(y_pos, (width - len(title)) // 2, title, curses.A_BOLD | curses.color_pair(4))
        y_pos += 1
        stdscr.addstr(y_pos, (width - len(work_info)) // 2, work_info, curses.color_pair(5))
        y_pos += 2
        
        # Overall progress bar
        overall_bar = self.format_progress_bar(
            progress.global_execution["Finished"],
            progress.global_execution["Total"],
            width - 20,
            True
        )
        overall_label = "Overall Progress:"
        stdscr.addstr(y_pos, 2, overall_label, curses.A_BOLD)
        stdscr.addstr(y_pos, len(overall_label) + 3, overall_bar, curses.color_pair(3))
        y_pos += 2
        
        # Hierarchical progress in a nice table format
        stdscr.addstr(y_pos, 0, "┌" + "─" * (width - 2) + "┐")
        y_pos += 1
        
        # Headers
        hierarchy_header = "│ Component      │ Current Status                                        │"
        stdscr.addstr(y_pos, 0, hierarchy_header[:width], curses.A_BOLD)
        y_pos += 1
        
        stdscr.addstr(y_pos, 0, "├" + "─" * 15 + "┼" + "─" * (width - 18) + "┤")
        y_pos += 1
        
        # Task row
        current_task = progress.tasks.get("Running", "") or ""
        task_line = f"│ 📋 Task        │ {current_task:<{width-20}} │"
        stdscr.addstr(y_pos, 0, task_line[:width])
        y_pos += 1
        
        # Dataset row  
        current_dataset = progress.datasets.get("Running", "") or ""
        dataset_line = f"│ 📊 Dataset     │ {current_dataset:<{width-20}} │"
        stdscr.addstr(y_pos, 0, dataset_line[:width])
        y_pos += 1
        
        # Config row
        current_config = progress.configs.get("Running", "") or ""
        config_line = f"│ ⚙️  Config      │ {current_config:<{width-20}} │"
        stdscr.addstr(y_pos, 0, config_line[:width])
        y_pos += 1
        
        # Algorithm row
        current_algorithm = progress.algorithms.get("Running", "") or ""
        algorithm_line = f"│ 🔬 Algorithm   │ {current_algorithm:<{width-20}} │"
        stdscr.addstr(y_pos, 0, algorithm_line[:width])
        y_pos += 1
        
        # Sequence row
        current_sequence = (progress.execution.get("Finished", 0)
                            + progress.execution.get("Running", 0))
        total_sequences = progress.execution.get("Total", 0)
        seq_progress = f"{current_sequence}/{total_sequences}"
        sequence_line = f"│ 🔄 Sequences   │ {seq_progress:<{width-20}} │"
        stdscr.addstr(y_pos, 0, sequence_line[:width])
        y_pos += 1
        
        stdscr.addstr(y_pos, 0, "└" + "─" * (width - 2) + "┘")
        y_pos += 2
        
        return y_pos
    
    def draw_execution_table(self, stdscr, executions: List[ExecutionDetail], warnings: List[dict], y_pos: int) -> int:
        """Draw an enhanced execution table."""
        height, width = stdscr.getmaxyx()
        
        if not executions:
            no_exec_msg = "📝 No active executions at the moment"
            stdscr.addstr(y_pos, (width - len(no_exec_msg)) // 2, no_exec_msg, curses.color_pair(5))
            return y_pos + 2
        
        # Table title
        table_title = "🔄 Active Executions"
        stdscr.addstr(y_pos, 2, table_title, curses.A_BOLD | curses.color_pair(4))
        y_pos += 2
        
        # Table header with enhanced design
        table_width = min(width - 2, 90)
        
        stdscr.addstr(y_pos, 0, "┌" + "─" * (table_width - 2) + "┐")
        y_pos += 1
        
        header_line = f"│ SEQ │ Unit ID      │ Progress                    │ Details                 │"
        stdscr.addstr(y_pos, 0, header_line[:table_width], curses.A_BOLD)
        y_pos += 1
        
        stdscr.addstr(y_pos, 0, "├" + "─" * (table_width - 2) + "┤")
        y_pos += 1
        
        # Execution rows with enhanced formatting
        max_rows = min(len(executions), (height - y_pos - 8))
        
        for i, exec_detail in enumerate(executions[:max_rows]):
            if y_pos >= height - 6:
                break
                
            # Main execution row
            progress_value = int(exec_detail.progress * 100) if exec_detail.progress else 0
            progress_bar = self.format_progress_bar(progress_value, 100, 15, False)
            
            seq_str = f"{exec_detail.sequencia:3d}"
            unit_id_short = exec_detail.unit_id[:12] if exec_detail.unit_id else "Unknown"
            
            # Status emoji
            status_emoji = "🟢" if exec_detail.status == "running" else "⏸️"
            
            # Combination info  
            combo_info = f"{exec_detail.algorithm_id[:8]}..."
            
            main_line = f"│ {seq_str} │ {unit_id_short:<12} │ {progress_bar:<27} │ {status_emoji} {combo_info:<19} │"
            stdscr.addstr(y_pos, 0, main_line[:table_width])
            y_pos += 1
            
            # Progress message line (if exists)
            if exec_detail.progress_message and y_pos < height - 5:
                msg = exec_detail.progress_message[:table_width-20] if len(exec_detail.progress_message) > table_width-20 else exec_detail.progress_message
                progress_msg_line = f"│     │ Progress: {msg:<{table_width-15}} │"
                stdscr.addstr(y_pos, 0, progress_msg_line[:table_width], curses.color_pair(5))
                y_pos += 1
            
            # Add warnings for this execution
            exec_warnings = [w for w in warnings if w.get('unit_id') == exec_detail.unit_id]
            for warning in exec_warnings[:1]:  # Show only first warning per execution
                if y_pos < height - 4:
                    warning_msg = warning.get('message', 'Unknown warning')[:table_width-20]
                    warning_line = f"│     │ ⚠️  {warning_msg:<{table_width-15}} │"
                    stdscr.addstr(y_pos, 0, warning_line[:table_width], curses.color_pair(2))
                    y_pos += 1
            
            # Separator line between executions  
            if i < len(executions[:max_rows]) - 1 and y_pos < height - 4:
                stdscr.addstr(y_pos, 0, "├" + "─" * (table_width - 2) + "┤")
                y_pos += 1
        
        # Table bottom
        if y_pos < height - 3:
            stdscr.addstr(y_pos, 0, "└" + "─" * (table_width - 2) + "┘")
            y_pos += 1
        
        return y_pos + 1
    
    def draw_errors(self, stdscr, errors: List[ErrorSummary], y_pos: int) -> int:
        """Draw enhanced error section."""
        height, width = stdscr.getmaxyx()
        
        if errors and y_pos < height - 2:
            error_title = "❌ Recent Errors"
            stdscr.addstr(y_pos, 2, error_title, curses.A_BOLD | curses.color_pair(1))
            y_pos += 1
            
            for i, error in enumerate(errors[:3]):  # Show only first 3 errors
                if y_pos >= height - 1:
                    break
                    
                error_prefix = f"  🔸 Unit {error.unit_id[:8]}:"
                error_msg = error.error_message
                max_msg_len = width - len(error_prefix) - 5
                
                if len(error_msg) > max_msg_len:
                    error_msg = error_msg[:max_msg_len-3] + "..."
                
                full_error = f"{error_prefix} {error_msg}"
                stdscr.addstr(y_pos, 0, full_error[:width-1], curses.color_pair(1))
                y_pos += 1
        
        return y_pos
    
    def draw_footer(self, stdscr, y_pos: int):
        """Draw enhanced footer with controls."""
        height, width = stdscr.getmaxyx()
        
        if y_pos < height - 1:
            # Controls and timestamp
            controls = "⌨️  Controls: [Q]uit | [C]ancel Work"
            timestamp = datetime.now().strftime("🕐 %H:%M:%S")
            refresh_info = "🔄 Auto-refresh: 2s"
            
            # Center the footer info
            footer_content = f"{controls} | {refresh_info} | {timestamp}"
            
            if len(footer_content) <= width - 1:
                stdscr.addstr(height - 1, (width - len(footer_content)) // 2, footer_content, curses.color_pair(5))
            else:
                # Truncate if too long
                stdscr.addstr(height - 1, 0, footer_content[:width-1], curses.color_pair(5))
    
    def run(self, stdscr):
        """Main monitoring loop."""
        # Setup curses
        curses.curs_set(0)  # Hide cursor
        stdscr.nodelay(1)   # Non-blocking input
        stdscr.timeout(100) # 100ms timeout for getch()
        
        # Setup colors
        curses.start_color()
        curses.init_pair(1, curses.COLOR_RED, curses.COLOR_BLACK)      # Errors
        curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)   # Warnings
        curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)    # Success/Progress
        curses.init_pair(4, curses.COLOR_CYAN, curses.COLOR_BLACK)     # Titles
        curses.init_pair(5, curses.COLOR_BLUE, curses.COLOR_BLACK)     # Info text
        
        last_refresh = 0
        
        while self.running:
            current_time = time.time()
            
            # Check for user input
            key = stdscr.getch()
            if key == ord('q') or key == 27:  # 'q' or ESC
                break
            elif key == ord('c') or key == ord('C'):  # 'c' to cancel
                # Show enhanced confirmation dialog
                stdscr.clear()
                stdscr.addstr(0, 0, "┌" + "─" * 60 + "┐")
                stdscr.addstr(1, 0, "│" + " " * 60 + "│")
                stdscr.addstr(2, 0, f"│  ⚠️  Cancel work '{self.work_id[:20]}{'...' if len(self.work_id) > 20 else ''}'?" + " " * (35 - len(self.work_id[:20])) + "│")
                stdscr.addstr(3, 0, "│" + " " * 60 + "│")
                stdscr.addstr(4, 0, "│  Press [Y] to confirm cancellation or any other key    │")
                stdscr.addstr(5, 0, "│  to continue monitoring...                             │")
                stdscr.addstr(6, 0, "│" + " " * 60 + "│")
                stdscr.addstr(7, 0, "└" + "─" * 60 + "┘")
                stdscr.refresh()
                
                confirm_key = stdscr.getch()
                while confirm_key == -1:  # Wait for input
                    time.sleep(0.1)
                    confirm_key = stdscr.getch()
                
                if confirm_key == ord('y') or confirm_key == ord('Y'):
                    try:
                        if self.work_service.cancel(self.work_id):
                            stdscr.clear()
                            stdscr.addstr(0, 0, "✅ Work cancelled successfully!", curses.color_pair(3))
                            stdscr.addstr(1, 0, "Monitor will exit in 2 seconds...")
                            stdscr.refresh()
                            time.sleep(2)
                            break
                        else:
                            stdscr.clear()
                            stdscr.addstr(0, 0, "❌ Failed to cancel work", curses.color_pair(1))
                            stdscr.addstr(1, 0, "Press any key to continue...")
                            stdscr.refresh()
                            stdscr.getch()
                    except Exception as e:
                        stdscr.clear()
                        stdscr.addstr(0, 0, f"❌ Error cancelling work: {e}", curses.color_pair(1))
                        stdscr.addstr(1, 0, "Press any key to continue...")
                        stdscr.refresh()
                        stdscr.getch()
            
            # Refresh data periodically
            if current_time - last_refresh >= self.refresh_interval:
                try:
                    # Clear screen
                    stdscr.clear()
                    
                    # Get fresh data
                    progress = self.queries.get_work_progress_summary(self.work_id)
                    if not progress:
                        stdscr.addstr(0, 0, f"Work '{self.work_id}' not found!", curses.color_pair(1))
                        stdscr.addstr(1, 0, "Press 'q' to quit")
                        stdscr.refresh()
                        time.sleep(1)
                        continue
                    
                    # Check work status - exit if no longer running
                    work_status = self.queries.get_work_status(self.work_id)
                    if work_status and work_status not in ['running', 'queued']:
                        # Mapear valores para tela final
                        finished = progress.global_execution["Finished"]
                        total = progress.global_execution["Total"]
                        overall_pct = progress.global_progress
                        stdscr.clear()
                        stdscr.addstr(0, 0, "┌" + "─" * 70 + "┐")
                        stdscr.addstr(1, 0, "│" + " " * 70 + "│")
                        stdscr.addstr(2, 0, f"│  🏁 Work completed with status: {work_status:<20}         │", curses.color_pair(3))
                        stdscr.addstr(3, 0, "│" + " " * 70 + "│")
                        stdscr.addstr(4, 0, f"│  📊 Final progress: {finished:,}/{total:,} executions completed" + " " * (70 - 45 - len(f"{finished:,}/{total:,}")) + "│")
                        stdscr.addstr(5, 0, f"│  🎯 Overall success rate: {overall_pct:.1%}" + " " * (70 - 30 - len(f"{overall_pct:.1%}")) + "│")
                        stdscr.addstr(6, 0, "│" + " " * 70 + "│")
                        stdscr.addstr(7, 0, "│  Monitor will exit in 3 seconds... (press any key to exit now)  │")
                        stdscr.addstr(8, 0, "│" + " " * 70 + "│")
                        stdscr.addstr(9, 0, "└" + "─" * 70 + "┘")
                        stdscr.refresh()
                        
                        # Wait 3 seconds or until key press
                        for _ in range(30):  # 30 * 0.1s = 3s
                            if stdscr.getch() != -1:  # Key pressed
                                break
                            time.sleep(0.1)
                        break
                    
                    executions = self.queries.get_running_executions_detail(self.work_id, limit=10)
                    errors = self.queries.get_error_summary(self.work_id, limit=5)
                    warnings = self.queries.get_execution_warnings(self.work_id, limit=10)
                    
                    # Draw interface
                    y_pos = 0
                    y_pos = self.draw_header(stdscr, progress, y_pos)
                    y_pos = self.draw_execution_table(stdscr, executions, warnings, y_pos)
                    y_pos = self.draw_errors(stdscr, errors, y_pos)
                    self.draw_footer(stdscr, y_pos)
                    
                    # Refresh screen
                    stdscr.refresh()
                    last_refresh = current_time
                    
                except Exception as e:
                    stdscr.clear()
                    stdscr.addstr(0, 0, f"Error: {str(e)}", curses.color_pair(1))
                    stdscr.addstr(1, 0, "Press 'q' to quit")
                    stdscr.refresh()
                    time.sleep(1)
            
            # Small sleep to prevent high CPU usage
            time.sleep(0.1)
    
    def start(self):
        """Start the monitor."""
        # Setup signal handler
        signal.signal(signal.SIGINT, self.signal_handler)
        
        try:
            # Verify work exists
            if not self.queries.work_exists(self.work_id):
                print(f"Error: Work '{self.work_id}' not found!")
                return 1
            
            # Start curses interface
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
    """CLI entry point - for backward compatibility."""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description="CSPBench Progress Monitor")
    parser.add_argument("work_id", help="Work ID to monitor")
    parser.add_argument("--db", default="data/work_manager.db", help="Database path (deprecated, auto-detected)")
    
    args = parser.parse_args()
    
    print("Note: Auto-detecting state database from work_service")
    
    try:
        # Start monitor
        monitor = ProgressMonitor(args.work_id)
        return monitor.start()
    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
