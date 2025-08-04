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
        
        # Terminal capabilities
        self._terminal_width = self._get_terminal_width()
        
        # State tracking
        self._batch_info = {}
        self._executions = {}
        self._current_execution = None
        self._current_dataset = None
        self._current_config = None
        self._algorithms = {}
        self._start_time = None
        self._last_update = datetime.now()
        
        # Display control
        self._display_lock = threading.Lock()
        self._header_printed = False
        self._last_display_state = {}  # Track what was last displayed to show only changes
        
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
    
    def _handle_task_started(self, event) -> None:
        """Handle TaskStartedEvent."""
        self._start_time = datetime.now()
        self._batch_info = {
            'name': event.task_name,
            'type': event.task_type.value,
            'metadata': event.metadata or {}
        }
        
        # Print initial header
        self._print_header()
    
    def _handle_execution_started(self, event) -> None:
        """Handle ExecutionStartedEvent."""
        self._current_execution = event.execution_name
        metadata = getattr(event, 'metadata', {})
        
        if self._current_execution not in self._executions:
            self._executions[self._current_execution] = {
                'name': event.execution_name,
                'status': 'running',
                'datasets': {},
                'metadata': metadata,
                'current_dataset': None,
                'current_config': None
            }
        
        self._refresh_display()
    
    def _handle_execution_progress(self, event) -> None:
        """Handle ExecutionProgressEvent."""
        if not self._current_execution:
            return
            
        context = event.context or {}
        
        # Extract hierarchical information from context
        dataset_name = context.get('dataset_name', event.item_name)
        config_name = context.get('algorithm_config_name', 'Default')
        
        # Update current tracking
        self._current_dataset = dataset_name
        self._current_config = config_name
        
        # Update execution structure
        execution = self._executions[self._current_execution]
        if dataset_name not in execution['datasets']:
            execution['datasets'][dataset_name] = {
                'name': dataset_name,
                'configs': {},
                'current_config': None
            }
        
        dataset = execution['datasets'][dataset_name]
        if config_name not in dataset['configs']:
            dataset['configs'][config_name] = {
                'name': config_name,
                'algorithms': {},
                'progress': 0
            }
        
        # Update progress
        dataset['configs'][config_name]['progress'] = event.progress_percent
        dataset['current_config'] = config_name
        execution['current_dataset'] = dataset_name
        
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
                self._executions[self._current_execution] = {
                    'name': self._current_execution,
                    'status': 'running',
                    'datasets': {},
                    'metadata': {},
                    'current_dataset': None,
                    'current_config': None
                }
        
        # Usar dataset padrão se não estiver definido
        dataset_name = self._current_dataset or "default_dataset"
        config_name = self._current_config or "default_config"
        
        execution = self._executions[self._current_execution]
        
        # Ensure dataset exists
        if dataset_name not in execution['datasets']:
            execution['datasets'][dataset_name] = {
                'name': dataset_name,
                'configs': {}
            }
        
        dataset = execution['datasets'][dataset_name]
        
        # Ensure config exists
        if config_name not in dataset['configs']:
            dataset['configs'][config_name] = {
                'name': config_name,
                'algorithms': {}
            }
        
        config = dataset['configs'][config_name]
        
        # Extract run information from message if available
        current_run = 1
        total_runs = 1
        
        # Try to parse run info from message like "Completed run 3/5" 
        if event.message:
            import re
            run_match = re.search(r'run (\d+)/(\d+)', event.message)
            if run_match:
                current_run = int(run_match.group(1))
                total_runs = int(run_match.group(2))
        
        # Also check context for run information
        context = getattr(event, 'context', {}) or {}
        current_run = context.get('current_repetition', context.get('current_run', current_run))
        total_runs = context.get('total_repetitions', context.get('total_runs', total_runs))
        
        config['algorithms'][event.algorithm_name] = {
            'name': event.algorithm_name,
            'progress': event.progress_percent,
            'current_run': current_run,
            'total_runs': total_runs,
            'status': 'running' if event.progress_percent < 100 else 'completed',
            'message': event.message
        }
        
        self._refresh_display()
    
    def _handle_execution_finished(self, event) -> None:
        """Handle ExecutionFinishedEvent."""
        if self._current_execution in self._executions:
            self._executions[self._current_execution]['status'] = 'completed' if event.success else 'failed'
        
        self._refresh_display()
    
    def _handle_task_finished(self, event) -> None:
        """Handle TaskFinishedEvent."""
        # Print final summary
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
        
        batch_name = self._batch_info.get('name', 'Unknown Batch')
        start_time = self._start_time.strftime("%Y-%m-%d %H:%M:%S") if self._start_time else "Unknown"
        
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
                # New execution
                exec_status = "🔄" if execution.get('status') == 'running' else "✅"
                print(f"\n{exec_status} Execution: {exec_name}")
                self._last_display_state[exec_name] = {'datasets': {}}
            
            exec_state = self._last_display_state[exec_name]
            
            if execution.get('datasets'):
                for dataset_name, dataset in execution['datasets'].items():
                    if dataset_name not in exec_state['datasets']:
                        # New dataset
                        print(f"  📊 Dataset: {dataset_name}")
                        exec_state['datasets'][dataset_name] = {'configs': {}}
                    
                    dataset_state = exec_state['datasets'][dataset_name]
                    
                    if dataset.get('configs'):
                        for config_name, config in dataset['configs'].items():
                            if config_name not in dataset_state['configs']:
                                # New configuration
                                print(f"    ⚙️  Configuration: {config_name}")
                                dataset_state['configs'][config_name] = {'algorithms': {}}
                            
                            config_state = dataset_state['configs'][config_name]
                            
                            if config.get('algorithms'):
                                for algo_name, algo in config['algorithms'].items():
                                    algo_key = f"{algo_name}_{algo['current_run']}_{algo['total_runs']}"
                                    
                                    if algo_key not in config_state['algorithms'] or config_state['algorithms'][algo_key] != algo['progress']:
                                        # New algorithm or progress change
                                        status_icon = self._get_algorithm_status_icon(algo['status'])
                                        progress_bar = self._create_progress_bar(algo['progress'])
                                        run_info = f"(Run {algo['current_run']}/{algo['total_runs']})"
                                        
                                        print(f"      {status_icon} {algo_name} {run_info}: [{progress_bar}] {algo['progress']:5.1f}%")
                                        
                                        if algo.get('message'):
                                            print(f"        💬 {algo['message']}")
                                        
                                        config_state['algorithms'][algo_key] = algo['progress']
    
    def _get_algorithm_status_icon(self, status: str) -> str:
        """Get status icon for algorithm."""
        icons = {
            'running': '🔄',
            'completed': '✅',
            'failed': '❌',
            'pending': '⏳'
        }
        return icons.get(status, '❓')
    
    def _create_progress_bar(self, progress: float, width: int = 20) -> str:
        """Create a progress bar string."""
        filled = int(progress / 100 * width)
        bar = "█" * filled + "░" * (width - filled)
        return bar
    
    def _print_final_summary(self, event) -> None:
        """Print final execution summary."""
        print("\\n✅ Task completed successfully!")
        
        if self._start_time:
            duration = datetime.now() - self._start_time
            total_seconds = int(duration.total_seconds())
            minutes = total_seconds // 60
            seconds = total_seconds % 60
            print(f"⏰ Total time: {minutes}m {seconds}s")
        
        print("=" * 60)
        
        # Print results summary if available
        results = event.results or {}
        if 'successful' in results and 'failed' in results:
            print("✅ Batch completed:")
            print(f"   Successful: {results['successful']}")
            print(f"   Failed: {results['failed']}")
        elif 'processed' in results:
            print(f"✅ Batch completed: {results['processed']} items processed")
        else:
            print("✅ Batch completed successfully!")
