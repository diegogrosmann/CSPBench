"""Advanced terminal display with live hierarchical monitoring."""

import sys
import os
from typing import Optional, Any, Dict, List
from datetime import datetime
import threading
import time
import logging


class LiveTerminalDisplay:
    """
    Advanced terminal display with live hierarchical monitoring.

    Shows execution progress in a structured hierarchical format with
    real-time updates and progress bars.
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize live terminal display.

        Args:
            verbose: Whether to show detailed progress
        """
        self._verbose = verbose
        self._logger = logging.getLogger(__name__)
        # Compact mode: show only progress lines (no headers/footers)
        # Default now is non-compact (hierarchical display). Set CSPB_COMPACT_DISPLAY=1 to enable compact.
        self._compact = os.getenv("CSPB_COMPACT_DISPLAY", "0") not in {"0", "false", "False"}

        # Terminal capabilities
        self._terminal_width = self._get_terminal_width()

        # State tracking
        self._batch_info: Dict[str, Any] = {}
        self._executions: Dict[str, Any] = {}
        self._current_execution: Optional[str] = None
        self._current_dataset: Optional[str] = None
        self._current_config: Optional[str] = None
        self._algorithms: Dict[str, Any] = {}
        self._start_time: Optional[datetime] = None
        self._last_update: datetime = datetime.now()

        # Tracks which algorithm line is currently open (being updated in-place)
        self._open_line_key: Optional[str] = None

        # Display control
        self._display_lock = threading.Lock()
        self._header_printed = False
        # Track what was last displayed to show only changes
        self._last_display_state: Dict[str, Any] = {}

    def _get_terminal_width(self) -> int:
        """Get terminal width or default to 80."""
        try:
            return os.get_terminal_size().columns
        except:
            return 80

    def handle_event(self, event: Any) -> None:
        """
        Handle events from the progress broker (compatibility method).

        Args:
            event: Progress event to process
        """
        self.on_event(event)

    def on_event(self, event) -> None:
        """
        Handle events from the progress broker.

        Args:
            event: Progress event to process
        """
        with self._display_lock:
            try:
                event_type = type(event).__name__

                # Unified display support
                if event_type == "DisplayEvent":
                    self._handle_display_event(event)
                    return

                if event_type == "TaskStartedEvent":
                    self._handle_task_started(event)
                elif event_type == "ExecutionStartedEvent":
                    self._handle_execution_started(event)
                elif event_type == "ExecutionProgressEvent":
                    self._handle_execution_progress(event)
                elif event_type == "ExecutionFinishedEvent":
                    self._handle_execution_finished(event)
                elif event_type == "AlgorithmProgressEvent":
                    self._handle_algorithm_progress(event)
                elif event_type == "TaskFinishedEvent":
                    self._handle_task_finished(event)
                elif event_type == "ErrorEvent":
                    self._handle_error(event)

            except Exception as e:
                self._logger.error(f"Error handling event {event_type}: {e}")


    def _handle_display_event(self, event) -> None:
        """Render unified DisplayEvent; for optimization and analysis, mirror hierarchical display."""
        phase = getattr(event, "phase", None)
        phase_name = str(getattr(phase, "value", phase)).lower()

        # Tratar optimization e analysis (sensibilidade) como hierarquia
        if phase_name not in {"optimization", "analysis"}:
            if self._compact:
                return
            print(f"📝 {getattr(event, 'message', '')}")
            return

        # Extract optimization context
        dataset = getattr(event, "dataset_id", "-")
        algo = getattr(event, "algorithm_name", "-")
        task_id = getattr(event, "task_id", "-")
        trial = int(getattr(event, "trial_no", 0) or 0)
        rep = int(getattr(event, "rep_idx", 0) or 0)
        progress = float(getattr(event, "progress", 0.0) or 0.0)
        msg = getattr(event, "message", "") or ""
        payload = getattr(event, "payload", {}) or {}

        total_trials = int(payload.get("total_trials") or 0)
        repetitions = int(payload.get("repetitions") or 0)
        # Optional hierarchical counters provided by executor
        task_index = int(payload.get("task_index", 1))
        total_tasks = int(payload.get("total_tasks", 1))
        task_name = payload.get("task_name") or task_id or "Task"
        ds_index = int(payload.get("dataset_index", 1))
        total_ds = int(payload.get("total_datasets", 1))
        cfg_index = int(payload.get("config_index", 1))
        total_cfgs = int(payload.get("total_configs", 1))

        # Skip cosmetic-noise events
        if msg == "algorithm skipped":
            return

        # Initialize header when first optimization event arrives
        if not self._header_printed and not self._compact:
            self._print_header()

        # Maintain a structure mirroring processing for optimization
        exec_name = f"Optimization: {task_name}"
        if exec_name not in self._executions:
            self._executions[exec_name] = {
                "name": exec_name,
                "status": "running",
                "datasets": {},
                "metadata": {"type": "optimization"},
                "current_dataset": None,
                "current_config": None,
                "index": task_index or 1,
                "total_executions": total_tasks or 1,
            }

        self._current_execution = exec_name

        # Dataset level
        execution = self._executions[exec_name]
        if dataset not in execution["datasets"]:
            execution["datasets"][dataset] = {
                "name": dataset,
                "configs": {},
                "current_config": None,
                # Use indices from payload when available for consistent counters
                "index": ds_index,
                "total_datasets": total_ds,
            }
        d = execution["datasets"][dataset]
        # Keep indices updated if later events provide better totals
        d["index"] = ds_index or d.get("index", 1)
        d["total_datasets"] = total_ds or d.get("total_datasets", 1)

        # Treat algorithm as a config entry under a pseudo-config "Default"
        cfg_name = "Default Configuration"
        if cfg_name not in d["configs"]:
            d["configs"][cfg_name] = {
                "name": cfg_name,
                "algorithms": {},
                "progress": 0,
                "index": cfg_index,
                "total_configs": total_cfgs,
            }
        d["current_config"] = cfg_name
        self._current_dataset = dataset
        self._current_config = cfg_name

        cfg = d["configs"][cfg_name]

        # Initialize algorithm tracking
        # Ignore algorithms with zero trials to avoid 0/0 bars
        if total_trials <= 0 and msg in {"trial start", "trial rep", "trial end"}:
            return

        if algo not in cfg["algorithms"]:
            cfg["algorithms"][algo] = {
                "name": algo,
                "progress": 0,
                "current_run": 0,
                "total_runs": total_trials,
                "status": "running",
                "message": "",
            }

        algo_info = cfg["algorithms"][algo]

        # Update current_run based on trial index
        if msg == "trial end" and trial > 0:
            # trial indices start at 1
            algo_info["current_run"] = max(algo_info.get("current_run", 0), trial)

        # Determine percent by trials completed
        # Update total_runs if the incoming total_trials is higher/better info
        if total_trials and total_trials > int(algo_info.get("total_runs", 0)):
            algo_info["total_runs"] = int(total_trials)
        total_runs = max(1, int(algo_info.get("total_runs") or total_trials or 1))
        current_runs = int(algo_info.get("current_run", 0))
        pct = (current_runs / total_runs) * 100.0

        # Completion status
        if current_runs >= total_runs:
            algo_info["status"] = "completed"
            pct = 100.0

        # Store message
        algo_info["message"] = msg

        # Compact vs non-compact rendering
        if self._compact:
            # Single line similar to processing output when finishing a trial rep
            bar = self._create_progress_bar(pct)
            print(
                f"      {'✅' if pct >= 100 else '🔄'} {algo} ({current_runs}/{total_runs}): [{bar}] {pct:5.1f}%"
            )
        else:
            # Reuse hierarchical incremental printer
            self._refresh_display()

    def _handle_task_started(self, event) -> None:
        """Handle TaskStartedEvent."""
        self._start_time = datetime.now()
        self._batch_info = {
            "name": event.task_name,
            "type": event.task_type.value,
            "metadata": event.metadata or {},
        }

        # Print initial header unless compact mode
        if not self._compact:
            self._print_header()

    def _handle_execution_started(self, event) -> None:
        """Handle ExecutionStartedEvent."""
        self._current_execution = event.execution_name
        metadata = getattr(event, "metadata", {})

        exec_key = self._current_execution or "Default Execution"
        if exec_key not in self._executions:
            self._executions[exec_key] = {
                "name": event.execution_name,
                "status": "running",
                "datasets": {},
                "metadata": metadata,
                "current_dataset": None,
                "current_config": None,
                # Extract counter information from metadata
                "index": metadata.get("index", 1),
                "total_executions": metadata.get("total_executions", 1),
            }

        if not self._compact:
            self._refresh_display()

    def _handle_execution_progress(self, event) -> None:
        """Handle ExecutionProgressEvent."""
        if not self._current_execution:
            return

        context = event.context or {}

        # Extract hierarchical information from context
        dataset_name = context.get("dataset_name", event.item_name)
        config_name = context.get("algorithm_config_name", "Default")

        # Extract counter information from context
        dataset_index = context.get("dataset_index", 1)
        total_datasets = context.get("total_datasets", 1)
        config_index = context.get("algorithm_config_index", 1)
        total_configs = context.get("total_algorithm_configs", 1)

        # Update current tracking
        self._current_dataset = dataset_name
        self._current_config = config_name

        # Update execution structure
        exec_key = self._current_execution or "Default Execution"
        execution = self._executions.get(exec_key)
        if execution is None:
            self._executions[exec_key] = {
                "name": exec_key,
                "status": "running",
                "datasets": {},
                "metadata": {},
                "current_dataset": None,
                "current_config": None,
                "index": 1,
                "total_executions": 1,
            }
            execution = self._executions[exec_key]

        if dataset_name not in execution["datasets"]:
            execution["datasets"][dataset_name] = {
                "name": dataset_name,
                "configs": {},
                "current_config": None,
                "index": dataset_index,
                "total_datasets": total_datasets,
            }

        dataset = execution["datasets"][dataset_name]
        if config_name not in dataset["configs"]:
            dataset["configs"][config_name] = {
                "name": config_name,
                "algorithms": {},
                "progress": 0,
                "index": config_index,
                "total_configs": total_configs,
            }

        # Update progress
        dataset["configs"][config_name]["progress"] = event.progress_percent
        dataset["current_config"] = config_name
        execution["current_dataset"] = dataset_name

        if self._compact:
            # Compact mode: print a single concise line per update
            try:
                print(
                    f"[execution] ds={dataset_name} cfg={config_name} p={event.progress_percent:.0f}%"
                )
                self._last_update = datetime.now()
            except Exception:
                pass
        else:
            self._refresh_display()

    def _handle_algorithm_progress(self, event) -> None:
        """Handle AlgorithmProgressEvent."""
        # Para o display funcionar, precisamos inferir dataset e config do contexto
        # Por enquanto, vamos usar valores padrão até termos essa informação

        if not self._current_execution:
            # Se não há execução atual, usa a primeira disponível ou cria uma padrão
            if self._executions:
                self._current_execution = list(self._executions.keys())[0]
            else:
                self._current_execution = "Default Execution"
                self._executions["Default Execution"] = {
                    "name": "Default Execution",
                    "status": "running",
                    "datasets": {},
                    "metadata": {},
                    "current_dataset": None,
                    "current_config": None,
                }

        # Usar dataset padrão se não estiver definido
        dataset_name = self._current_dataset or "default_dataset"
        config_name = self._current_config or "default_config"

        exec_key = self._current_execution or "Default Execution"
        execution = self._executions[exec_key]

        # Ensure dataset exists
        if dataset_name not in execution["datasets"]:
            execution["datasets"][dataset_name] = {"name": dataset_name, "configs": {}}

        dataset = execution["datasets"][dataset_name]

        # Ensure config exists
        if config_name not in dataset["configs"]:
            dataset["configs"][config_name] = {"name": config_name, "algorithms": {}}

        config = dataset["configs"][config_name]

        # Extract run information from message if available
        current_run = 1
        total_runs = 1

        # Try to parse run info from message like "Completed run 3/5"
        if event.message:
            import re

            run_match = re.search(r"run (\d+)/(\d+)", event.message)
            if run_match:
                current_run = int(run_match.group(1))
                total_runs = int(run_match.group(2))

        # Also check context for run information
        context = getattr(event, "context", {}) or {}
        current_run = context.get(
            "current_repetition", context.get("current_run", current_run)
        )
        total_runs = context.get(
            "total_repetitions", context.get("total_runs", total_runs)
        )

        # Update or create algorithm entry (aggregated)
        if event.algorithm_name not in config["algorithms"]:
            config["algorithms"][event.algorithm_name] = {
                "name": event.algorithm_name,
                "progress": 0,
                "current_run": 0,
                "total_runs": total_runs,
                "status": "running",
                "message": event.message,
            }

        # Update algorithm progress - only track the highest completed run
        algo_info = config["algorithms"][event.algorithm_name]

        # If this run is completed (100%), update current_run
        if event.progress_percent >= 100 and current_run > algo_info["current_run"]:
            algo_info["current_run"] = current_run

        # Update total_runs if we have better info
        if total_runs > algo_info["total_runs"]:
            algo_info["total_runs"] = total_runs

        # Update status based on completion
        if algo_info["current_run"] >= algo_info["total_runs"]:
            algo_info["status"] = "completed"
        else:
            algo_info["status"] = "running"

        # Update message
        algo_info["message"] = event.message

        if self._compact:
            # Compact mode: show a minimal per-algorithm progress line
            try:
                print(
                    f"[algorithm] alg={event.algorithm_name} runs={algo_info['current_run']}/{algo_info['total_runs']} status={algo_info['status']}"
                )
                self._last_update = datetime.now()
            except Exception:
                pass
        else:
            self._refresh_display()

    def _handle_execution_finished(self, event) -> None:
        """Handle ExecutionFinishedEvent."""
        exec_key = self._current_execution or "Default Execution"
        if exec_key in self._executions:
            self._executions[exec_key]["status"] = (
                "completed" if event.success else "failed"
            )

        if not self._compact:
            self._refresh_display()

    def _handle_task_finished(self, event) -> None:
        """Handle TaskFinishedEvent."""
        # Print final summary
        if not self._compact:
            self._print_final_summary(event)

    def _handle_error(self, event) -> None:
        """Handle ErrorEvent."""
        print(f"\\n❌ Error ({event.error_type}): {event.error_message}")

    def _print_header(self) -> None:
        """Print the monitoring header."""
        if self._header_printed:
            return

        self._header_printed = True

        print("=" * 60)
        print("🚀 CSPBench - Execution Monitoring")
        print("=" * 60)

        batch_name = self._batch_info.get("name", "Unknown Batch")
        start_time = (
            self._start_time.strftime("%Y-%m-%d %H:%M:%S")
            if self._start_time
            else "Unknown"
        )

        print(f"📋 Batch: {batch_name}")
        print(f"⏰ Started: {start_time}")
        print("=" * 60)
        print()

    def _refresh_display(self) -> None:
        """Refresh the live display showing only changes."""
        if not self._header_printed:
            self._print_status_header()
            self._header_printed = True
        # Show only incremental changes
        self._print_execution_changes()

    def _print_status_header(self) -> None:
        """Print the status header once."""
        print("\n" + "=" * 60)
        print("📋 LIVE EXECUTION STATUS")
        print("=" * 60)

    def _print_execution_changes(self) -> None:
        """Print only the changes since last update."""
        for exec_name, execution in self._executions.items():
            if exec_name not in self._last_display_state:
                # New execution - include counter information
                exec_status = "🔄" if execution.get("status") == "running" else "✅"
                exec_index = execution.get("index", 1)
                total_executions = execution.get("total_executions", 1)
                print(
                    f"\n{exec_status} Task: {exec_name} ({exec_index}/{total_executions})"
                )
                self._last_display_state[exec_name] = {"datasets": {}}

            exec_state = self._last_display_state[exec_name]

            if execution.get("datasets"):
                for dataset_name, dataset in execution["datasets"].items():
                    if dataset_name not in exec_state["datasets"]:
                        # Ensure any open algorithm line is closed before a new header
                        if self._open_line_key is not None:
                            print()
                            self._open_line_key = None
                        # New dataset - include counter information
                        dataset_index = dataset.get("index", 1)
                        total_datasets = dataset.get("total_datasets", 1)
                        print(
                            f"  📊 Dataset: {dataset_name} ({dataset_index}/{total_datasets})"
                        )
                        exec_state["datasets"][dataset_name] = {"configs": {}}

                    dataset_state = exec_state["datasets"][dataset_name]

                    if dataset.get("configs"):
                        for config_name, config in dataset["configs"].items():
                            if config_name not in dataset_state["configs"]:
                                # Ensure any open algorithm line is closed before a new header
                                if self._open_line_key is not None:
                                    print()
                                    self._open_line_key = None
                                # New configuration - include counter information
                                config_index = config.get("index", 1)
                                total_configs = config.get("total_configs", 1)
                                print(
                                    f"    ⚙️  Configuration: {config_name} ({config_index}/{total_configs})"
                                )
                                dataset_state["configs"][config_name] = {
                                    "algorithms": {}
                                }

                            config_state = dataset_state["configs"][config_name]

                            if config.get("algorithms"):
                                for algo_name, algo in config["algorithms"].items():
                                    # Calculate aggregated progress safely
                                    total_runs = max(0, int(algo.get("total_runs", 0)))
                                    current_runs = max(0, int(algo.get("current_run", 0)))
                                    completed_runs = (
                                        total_runs if total_runs > 0 and current_runs >= total_runs else current_runs
                                    )
                                    progress_percent = (
                                        (completed_runs / total_runs) * 100 if total_runs > 0 else 0
                                    )

                                    # Use a stable key per algorithm (avoid duplication when total_runs changes)
                                    algo_key = algo_name
                                    full_key = f"{exec_name}|{dataset_name}|{config_name}|{algo_name}"

                                    # Last printed progress for this algorithm
                                    last_progress = config_state["algorithms"].get(algo_key, -1)

                                    # If first time seeing the algorithm, print initial 0% line (without newline)
                                    if algo_key not in config_state["algorithms"]:
                                        # Ensure any previously open line is finalized
                                        if self._open_line_key is not None and self._open_line_key != full_key:
                                            print()
                                            self._open_line_key = None
                                        init_bar = self._create_progress_bar(0)
                                        print(
                                            f"      ⏳ {algo_name} (0/{total_runs}): [{init_bar}]   0.0%",
                                            end="",
                                            flush=True,
                                        )
                                        self._open_line_key = full_key
                                        config_state["algorithms"][algo_key] = 0

                                    # Update in-place while running; print newline on completion
                                    if 0 <= last_progress < 100 and progress_percent > last_progress and progress_percent < 100:
                                        status_icon = self._get_algorithm_status_icon_by_progress(progress_percent)
                                        progress_bar = self._create_progress_bar(progress_percent)
                                        run_info = f"({completed_runs}/{total_runs})"
                                        line_text = f"      {status_icon} {algo_name} {run_info}: [{progress_bar}] {progress_percent:5.1f}%"
                                        if self._open_line_key == full_key:
                                            self._overwrite_line(line_text)
                                        else:
                                            print("\n" + line_text, end="", flush=True)
                                            self._open_line_key = full_key
                                        config_state["algorithms"][algo_key] = progress_percent

                                    if progress_percent >= 100 and last_progress < 100:
                                        status_icon = "✅"
                                        progress_bar = self._create_progress_bar(100)
                                        run_info = f"({total_runs}/{total_runs})"
                                        line_text = f"      {status_icon} {algo_name} {run_info}: [{progress_bar}] 100.0%"
                                        if self._open_line_key == full_key:
                                            self._overwrite_line(line_text)
                                            print()
                                        else:
                                            print("\n" + line_text)
                                        self._open_line_key = None
                                        config_state["algorithms"][algo_key] = 100

    def _get_algorithm_status_icon(self, status: str) -> str:
        """Get status icon for algorithm."""
        icons = {"running": "🔄", "completed": "✅", "failed": "❌", "pending": "⏳"}
        return icons.get(status, "❓")

    def _get_algorithm_status_icon_by_progress(self, progress_percent: float) -> str:
        """Get status icon based on progress percentage."""
        if progress_percent >= 100:
            return "✅"
        elif progress_percent > 0:
            return "🔄"
        else:
            return "⏳"

    def _create_progress_bar(self, percentage: float, width: int = 20) -> str:
        """Create a progress bar string."""
        filled = int(percentage / 100 * width)
        bar = "█" * filled + "░" * (width - filled)
        return bar

    def _overwrite_line(self, text: str) -> None:
        """Overwrite the current terminal line, clearing residual characters."""
        try:
            sys.stdout.write("\r\033[K" + text)
            sys.stdout.flush()
        except Exception:
            # Fallback if ANSI not supported
            print("\r" + text, end="", flush=True)

    def _print_final_summary(self, event) -> None:
        """Print final execution summary."""
        print("✅ Task completed successfully!")

        if self._start_time:
            duration = datetime.now() - self._start_time
            total_seconds = int(duration.total_seconds())
            minutes = total_seconds // 60
            seconds = total_seconds % 60
            print(f"⏰ Total time: {minutes}m {seconds}s")

        print("=" * 60)

        # Print results summary if available
        results = event.results or {}
        if "successful" in results and "failed" in results:
            print("✅ Batch completed:")
            print(f"   Successful: {results['successful']}")
            print(f"   Failed: {results['failed']}")
        elif "processed" in results:
            print(f"✅ Batch completed: {results['processed']} items processed")
        else:
            print("✅ Batch completed successfully!")
