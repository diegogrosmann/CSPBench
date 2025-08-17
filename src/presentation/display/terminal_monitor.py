"""Monitor de linha de comando moderno para CSPBench.

Este monitor oferece uma experiência de acompanhamento otimizada para execuções
no terminal, com exibição hierárquica e atualização em tempo real.
"""

import logging
import os
import sys
import threading
from collections import defaultdict
from datetime import datetime, timedelta
from typing import Any, Dict, Optional, Set

# Constantes para formatação
class Colors:
    """Códigos ANSI para cores no terminal."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    
    # Cores básicas
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    
    # Cores brilhantes
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'

class Symbols:
    """Símbolos Unicode para diferentes estados."""
    RUNNING = "🔄"
    COMPLETED = "✅"
    FAILED = "❌"
    PENDING = "⏳"
    WARNING = "⚠️"
    INFO = "ℹ️"
    START = "🚀"
    FINISH = "🏁"
    PROGRESS = "📊"
    TIME = "⏰"
    DATASET = "🗂️"
    ALGORITHM = "⚙️"
    REPETITION = "🔁"


class TerminalMonitor:
    """
    Monitor moderno para linha de comando.
    
    Oferece uma experiência otimizada de acompanhamento com:
    - Exibição hierárquica clara
    - Atualização em tempo real
    - Suporte a cores
    - Informações de progresso detalhadas
    """
    
    def __init__(self, verbose: bool = True, compact: bool = False, colors: bool = None):
        """
        Inicializa o monitor de terminal.
        
        Args:
            verbose: Se deve mostrar informações detalhadas
            compact: Se deve usar modo compacto (menos informações)
            colors: Se deve usar cores (auto-detecta se None)
        """
        self.verbose = verbose
        self.compact = compact
        self.colors_enabled = colors if colors is not None else self._detect_color_support()
        
        # Estado interno
        self._logger = logging.getLogger(__name__)
        self._lock = threading.RLock()
        self._start_time: Optional[datetime] = None
        
        # Controle de exibição
        self._terminal_width = self._get_terminal_width()
        self._last_line_length = 0
        self._header_shown = False
        
        # Estado da execução
        self._task_info: Dict[str, Any] = {}
        self._executions: Dict[str, Dict[str, Any]] = {}
        self._current_execution: Optional[str] = None
        self._active_items: Set[str] = set()
        
        # Cache para otimização
        self._cached_display: Dict[str, str] = {}
        self._last_update = datetime.now()
    
    def _detect_color_support(self) -> bool:
        """Detecta se o terminal suporta cores."""
        try:
            return (
                hasattr(sys.stdout, 'isatty') and sys.stdout.isatty() and
                os.getenv('TERM', '').lower() not in ('dumb', '') and
                os.getenv('NO_COLOR') is None
            )
        except Exception:
            return False
    
    def _get_terminal_width(self) -> int:
        """Obtém a largura do terminal."""
        try:
            return os.get_terminal_size().columns
        except Exception:
            return 80
    
    def _colorize(self, text: str, color: str) -> str:
        """Aplica cor ao texto se cores estão habilitadas."""
        if not self.colors_enabled:
            return text
        return f"{color}{text}{Colors.RESET}"
    
    def _clear_line(self) -> None:
        """Limpa a linha atual do terminal."""
        if self._last_line_length > 0:
            sys.stdout.write('\r' + ' ' * self._last_line_length + '\r')
            self._last_line_length = 0
    
    def _print_line(self, text: str, end: str = '\n') -> None:
        """Imprime uma linha no terminal."""
        with self._lock:
            self._clear_line()
            print(text, end=end, flush=True)
            if end == '':
                self._last_line_length = len(text)
    
    def _update_line(self, text: str) -> None:
        """Atualiza a linha atual (para indicadores de progresso)."""
        with self._lock:
            # Throttle updates para evitar flooding
            now = datetime.now()
            if (now - self._last_update).total_seconds() < 0.1:
                return
            
            self._clear_line()
            sys.stdout.write(text)
            sys.stdout.flush()
            self._last_line_length = len(text)
            self._last_update = now
    
    def _create_progress_bar(self, percentage: float, width: int = 30) -> str:
        """Cria uma barra de progresso visual."""
        filled = int(width * percentage / 100)
        empty = width - filled
        
        if self.colors_enabled:
            bar = (
                self._colorize("█" * filled, Colors.BRIGHT_GREEN) +
                self._colorize("░" * empty, Colors.DIM)
            )
        else:
            bar = "█" * filled + "░" * empty
        
        return f"[{bar}]"
    
    def _format_duration(self, duration: timedelta) -> str:
        """Formata duração de forma legível."""
        total_seconds = int(duration.total_seconds())
        hours = total_seconds // 3600
        minutes = (total_seconds % 3600) // 60
        seconds = total_seconds % 60
        
        if hours > 0:
            return f"{hours}h{minutes:02d}m{seconds:02d}s"
        elif minutes > 0:
            return f"{minutes}m{seconds:02d}s"
        else:
            return f"{seconds}s"
    
    def _print_header(self) -> None:
        """Imprime o cabeçalho do monitor."""
        if self._header_shown:
            return
        
        header_text = "CSPBench - Monitor de Execução"
        separator = "=" * len(header_text)
        
        self._print_line(self._colorize(separator, Colors.BRIGHT_BLUE))
        self._print_line(self._colorize(f"{Symbols.START} {header_text}", Colors.BRIGHT_BLUE))
        self._print_line(self._colorize(separator, Colors.BRIGHT_BLUE))
        
        if self._start_time:
            start_time = self._start_time.strftime("%Y-%m-%d %H:%M:%S")
            self._print_line(f"{Symbols.TIME} Iniciado em: {start_time}")
        
        self._print_line("")
        self._header_shown = True
    
    def handle_event(self, event: Any) -> None:
        """
        Manipula eventos de progresso.
        
        Args:
            event: Evento de progresso recebido
        """
        with self._lock:
            try:
                event_type = type(event).__name__
                
                # Mapeia eventos para handlers
                handler_map = {
                    'TaskStartedEvent': self._handle_task_started,
                    'TaskFinishedEvent': self._handle_task_finished,
                    'ExecutionStartedEvent': self._handle_execution_started,
                    'ExecutionProgressEvent': self._handle_execution_progress,
                    'ExecutionFinishedEvent': self._handle_execution_finished,
                    'ExperimentStartedEvent': self._handle_execution_started,  # Alias
                    'ExperimentProgressEvent': self._handle_execution_progress,  # Alias
                    'ExperimentFinishedEvent': self._handle_execution_finished,  # Alias
                    'AlgorithmProgressEvent': self._handle_algorithm_progress,
                    'ErrorEvent': self._handle_error,
                    'WarningEvent': self._handle_warning,
                    'DisplayEvent': self._handle_display_event,
                }
                
                handler = handler_map.get(event_type, self._handle_generic)
                handler(event)
                
            except Exception as e:
                self._logger.error(f"Erro ao processar evento {type(event).__name__}: {e}")
    
    def _handle_task_started(self, event) -> None:
        """Manipula início de tarefa."""
        self._start_time = datetime.now()
        
        self._task_info = {
            'name': event.task_name,
            'type': getattr(event.task_type, 'value', str(event.task_type)),
            'metadata': getattr(event, 'metadata', {}),
            'status': 'running'
        }
        
        if not self.compact:
            self._print_header()
            task_type = self._task_info['type'].title()
            self._print_line(
                f"{Symbols.START} Iniciando {self._colorize(task_type, Colors.BRIGHT_YELLOW)}: "
                f"{self._colorize(self._task_info['name'], Colors.WHITE)}"
            )
            
            if self.verbose and self._task_info['metadata']:
                for key, value in self._task_info['metadata'].items():
                    self._print_line(f"  {key}: {value}")
            
            self._print_line("")
    
    def _handle_task_finished(self, event) -> None:
        """Manipula conclusão de tarefa."""
        success = getattr(event, 'success', True)
        
        if success:
            symbol = Symbols.COMPLETED
            color = Colors.BRIGHT_GREEN
            status = "CONCLUÍDA"
        else:
            symbol = Symbols.FAILED
            color = Colors.BRIGHT_RED
            status = "FALHOU"
        
        self._task_info['status'] = 'completed' if success else 'failed'
        
        # Limpa linha de progresso se existir
        self._clear_line()
        
        duration = ""
        if self._start_time:
            elapsed = datetime.now() - self._start_time
            duration = f" em {self._format_duration(elapsed)}"
        
        self._print_line(
            f"\n{symbol} Tarefa {self._colorize(status, color)}: "
            f"{self._colorize(self._task_info['name'], Colors.WHITE)}{duration}"
        )
        
        # Mostra erro se houver
        error_msg = getattr(event, 'error_message', None)
        if error_msg:
            self._print_line(f"  Erro: {self._colorize(error_msg, Colors.RED)}")
        
        # Mostra resultados se disponíveis
        results = getattr(event, 'results', {})
        if results:
            self._print_line("\n" + self._colorize("Resultados:", Colors.BRIGHT_CYAN))
            for key, value in results.items():
                self._print_line(f"  {key}: {value}")
        
        self._print_line("")
    
    def _handle_execution_started(self, event) -> None:
        """Manipula início de execução."""
        execution_name = getattr(event, 'execution_name', 
                                getattr(event, 'experiment_name', 'Execução'))
        
        self._current_execution = execution_name
        
        metadata = getattr(event, 'metadata', {})
        total_items = getattr(event, 'total_items', 0)
        
        self._executions[execution_name] = {
            'name': execution_name,
            'status': 'running',
            'start_time': datetime.now(),
            'total_items': total_items,
            'completed_items': 0,
            'failed_items': 0,
            'current_item': '',
            'progress': 0.0,
            'metadata': metadata,
            'datasets': defaultdict(dict),
            'algorithms': defaultdict(dict)
        }
        
        if not self.compact:
            index = metadata.get('index', 1)
            total = metadata.get('total_executions', 1)
            counter = f"({index}/{total})" if total > 1 else ""
            
            self._print_line(
                f"{Symbols.PROGRESS} Execução: {self._colorize(execution_name, Colors.BRIGHT_CYAN)} {counter}"
            )
            
            if total_items > 0:
                self._print_line(f"  Total de itens: {total_items}")
    
    def _handle_execution_progress(self, event) -> None:
        """Manipula progresso de execução."""
        if not self._current_execution:
            return
        
        execution = self._executions.get(self._current_execution, {})
        if not execution:
            return
        
        # Atualiza informações de progresso
        execution['progress'] = getattr(event, 'progress_percent', 0.0)
        execution['current_item'] = getattr(event, 'item_name', '')
        
        context = getattr(event, 'context', {})
        dataset_name = context.get('dataset_name', execution['current_item'])
        algorithm_name = context.get('algorithm_name', '')
        config_name = context.get('algorithm_config_name', '')
        
        # Atualiza estruturas hierárquicas
        if dataset_name:
            execution['datasets'][dataset_name].update({
                'status': 'running',
                'progress': execution['progress']
            })
        
        if algorithm_name:
            execution['algorithms'][algorithm_name].update({
                'status': 'running',
                'progress': execution['progress'],
                'dataset': dataset_name,
                'config': config_name
            })
        
        # Exibe progresso
        if self.compact:
            # Modo compacto: linha única
            progress_bar = self._create_progress_bar(execution['progress'], 20)
            status = (
                f"{progress_bar} {execution['progress']:5.1f}% "
                f"- {dataset_name[:20]}{'...' if len(dataset_name) > 20 else ''}"
            )
            if algorithm_name:
                status += f" | {algorithm_name[:15]}{'...' if len(algorithm_name) > 15 else ''}"
            
            self._update_line(status)
        else:
            # Modo detalhado: exibição hierárquica
            self._display_hierarchical_progress(execution)
    
    def _handle_execution_finished(self, event) -> None:
        """Manipula conclusão de execução."""
        execution_name = getattr(event, 'execution_name',
                                getattr(event, 'experiment_name', self._current_execution))
        
        if not execution_name or execution_name not in self._executions:
            return
        
        execution = self._executions[execution_name]
        success = getattr(event, 'success', True)
        
        execution['status'] = 'completed' if success else 'failed'
        execution['end_time'] = datetime.now()
        
        # Limpa linha de progresso
        self._clear_line()
        
        # Calcula duração
        duration = ""
        if 'start_time' in execution:
            elapsed = execution['end_time'] - execution['start_time']
            duration = f" em {self._format_duration(elapsed)}"
        
        # Exibe resultado
        if success:
            symbol = Symbols.COMPLETED
            color = Colors.BRIGHT_GREEN
            status_text = "CONCLUÍDA"
        else:
            symbol = Symbols.FAILED
            color = Colors.BRIGHT_RED
            status_text = "FALHOU"
        
        self._print_line(
            f"\n{symbol} Execução {self._colorize(status_text, color)}: "
            f"{self._colorize(execution_name, Colors.WHITE)}{duration}"
        )
        
        # Mostra estatísticas finais
        if self.verbose and not self.compact:
            self._print_execution_summary(execution)
    
    def _handle_algorithm_progress(self, event) -> None:
        """Manipula progresso de algoritmo."""
        algorithm_name = getattr(event, 'algorithm_name', 'Unknown')
        progress = getattr(event, 'progress_percent', 0.0)
        message = getattr(event, 'message', '')
        
        if not self._current_execution:
            return
        
        execution = self._executions.get(self._current_execution, {})
        if not execution:
            return
        
        # Atualiza estado do algoritmo
        execution['algorithms'][algorithm_name].update({
            'progress': progress,
            'message': message,
            'status': 'completed' if progress >= 100 else 'running'
        })
        
        if self.verbose and not self.compact:
            # Exibe progresso detalhado do algoritmo
            progress_bar = self._create_progress_bar(progress, 20)
            status = f"    {Symbols.ALGORITHM} {algorithm_name}: {progress_bar} {progress:5.1f}%"
            if message:
                status += f" - {message[:30]}{'...' if len(message) > 30 else ''}"
            
            self._print_line(status)
    
    def _handle_display_event(self, event) -> None:
        """Manipula eventos de display unificados."""
        phase = getattr(event, 'phase', None)
        if phase:
            phase_name = str(getattr(phase, 'value', phase)).lower()
        else:
            phase_name = 'unknown'
        
        dataset = getattr(event, 'dataset_id', '')
        algorithm = getattr(event, 'algorithm_name', '')
        message = getattr(event, 'message', '')
        payload = getattr(event, 'payload', {})
        
        # Para eventos de otimização e análise
        if phase_name in ['optimization', 'analysis']:
            trial = int(getattr(event, 'trial_no', 0) or 0)
            total_trials = int(payload.get('total_trials', 0) or 0)
            
            if total_trials > 0 and message == 'trial end':
                progress = (trial / total_trials) * 100.0
                self._update_algorithm_progress(algorithm, progress, f"Trial {trial}/{total_trials}")
        
        # Para eventos de processamento
        elif phase_name == 'processing':
            rep = int(getattr(event, 'rep_idx', 0) or 0)
            repetitions = int(payload.get('repetitions', 0) or 0)
            
            if repetitions > 0 and message == 'repetition':
                progress = (rep / repetitions) * 100.0
                self._update_algorithm_progress(algorithm, progress, f"Rep {rep}/{repetitions}")
    
    def _update_algorithm_progress(self, algorithm: str, progress: float, message: str) -> None:
        """Atualiza progresso de um algoritmo específico."""
        if not self._current_execution:
            return
        
        execution = self._executions.get(self._current_execution, {})
        if not execution:
            return
        
        execution['algorithms'][algorithm].update({
            'progress': progress,
            'message': message,
            'status': 'completed' if progress >= 100 else 'running'
        })
        
        if self.compact:
            # Atualização compacta
            progress_bar = self._create_progress_bar(progress, 15)
            status = f"{algorithm[:20]}: {progress_bar} {progress:5.1f}%"
            self._update_line(status)
        elif self.verbose:
            # Atualização detalhada
            symbol = Symbols.COMPLETED if progress >= 100 else Symbols.RUNNING
            progress_bar = self._create_progress_bar(progress, 20)
            status = f"    {symbol} {algorithm}: {progress_bar} {progress:5.1f}% - {message}"
            self._print_line(status)
    
    def _handle_error(self, event) -> None:
        """Manipula eventos de erro."""
        error_type = getattr(event, 'error_type', 'Error')
        error_message = getattr(event, 'error_message', 'Unknown error')
        
        self._clear_line()
        self._print_line(
            f"{Symbols.FAILED} {self._colorize(error_type, Colors.BRIGHT_RED)}: "
            f"{self._colorize(error_message, Colors.RED)}"
        )
    
    def _handle_warning(self, event) -> None:
        """Manipula eventos de aviso."""
        warning_message = getattr(event, 'message', 'Unknown warning')
        
        if self.verbose:
            self._print_line(
                f"{Symbols.WARNING} {self._colorize('Aviso', Colors.BRIGHT_YELLOW)}: "
                f"{warning_message}"
            )
    
    def _handle_generic(self, event) -> None:
        """Manipula eventos genéricos."""
        if self.verbose:
            event_type = type(event).__name__
            message = getattr(event, 'message', 'No message')
            self._print_line(
                f"{Symbols.INFO} {self._colorize(event_type, Colors.CYAN)}: {message}"
            )
    
    def _display_hierarchical_progress(self, execution: Dict[str, Any]) -> None:
        """Exibe progresso em formato hierárquico."""
        # Para evitar spam, só atualiza a cada segundo
        now = datetime.now()
        if (now - self._last_update).total_seconds() < 1.0:
            return
        
        # Limpa e reconstrói exibição
        self._clear_line()
        
        progress_bar = self._create_progress_bar(execution['progress'], 30)
        main_status = (
            f"  {Symbols.PROGRESS} Progresso geral: {progress_bar} "
            f"{execution['progress']:5.1f}%"
        )
        
        if execution['current_item']:
            main_status += f" - {execution['current_item']}"
        
        self._print_line(main_status)
        
        # Mostra datasets ativos
        active_datasets = {k: v for k, v in execution['datasets'].items() 
                          if v.get('status') == 'running'}
        
        if active_datasets and len(active_datasets) <= 3:  # Evita spam
            for dataset_name, dataset_info in active_datasets.items():
                dataset_progress = dataset_info.get('progress', 0.0)
                dataset_bar = self._create_progress_bar(dataset_progress, 20)
                self._print_line(
                    f"    {Symbols.DATASET} {dataset_name}: {dataset_bar} "
                    f"{dataset_progress:5.1f}%"
                )
        
        self._last_update = now
    
    def _print_execution_summary(self, execution: Dict[str, Any]) -> None:
        """Imprime resumo da execução."""
        total_datasets = len(execution['datasets'])
        total_algorithms = len(execution['algorithms'])
        
        completed_algorithms = sum(
            1 for alg in execution['algorithms'].values()
            if alg.get('status') == 'completed'
        )
        
        self._print_line(f"  Datasets processados: {total_datasets}")
        self._print_line(f"  Algoritmos executados: {completed_algorithms}/{total_algorithms}")
        
        if execution.get('total_items', 0) > 0:
            self._print_line(f"  Itens concluídos: {execution['completed_items']}/{execution['total_items']}")
    
    def finish(self) -> None:
        """Finaliza o monitor e limpa a tela."""
        with self._lock:
            self._clear_line()
            
            if self._start_time and not self.compact:
                total_duration = datetime.now() - self._start_time
                self._print_line(
                    f"\n{Symbols.FINISH} Monitor finalizado após "
                    f"{self._format_duration(total_duration)}"
                )

    # ========================================================================
    # INTERFACE METHODS FOR MONITOR PROTOCOL COMPATIBILITY
    # ========================================================================
    
    def pipeline_started(self, name: str, total_tasks: int) -> None:
        """Implementa Monitor.pipeline_started para compatibilidade."""
        self._start_time = datetime.now()
        self._print_header()
        self._print_line(
            f"{Symbols.START} Pipeline iniciado: {self._colorize(name, Colors.BOLD)}"
        )
        self._print_line(f"  Total de tarefas: {total_tasks}")
        self._print_line("")
    
    def pipeline_finished(self, success: bool, summary: Optional[Dict[str, Any]] = None) -> None:
        """Implementa Monitor.pipeline_finished para compatibilidade."""
        if self._start_time:
            duration = datetime.now() - self._start_time
            duration_str = self._format_duration(duration)
        else:
            duration_str = "desconhecido"
            
        if success:
            self._print_line(
                f"\n{Symbols.COMPLETED} Pipeline concluído com sucesso em {duration_str}"
            )
        else:
            self._print_line(
                f"\n{Symbols.FAILED} Pipeline falhou após {duration_str}"
            )
        
        if summary and not self.compact:
            self._print_line("  Resumo:")
            for key, value in summary.items():
                self._print_line(f"    {key}: {value}")
    
    def task_started(self, task_id: str, meta: Dict[str, Any]) -> None:
        """Implementa Monitor.task_started para compatibilidade."""
        task_name = meta.get('name', task_id)
        task_type = meta.get('type', 'task')
        self._print_line(
            f"{Symbols.RUNNING} Iniciando {task_type}: {self._colorize(task_name, Colors.CYAN)}"
        )
        if meta.get('description') and self.verbose:
            self._print_line(f"  {meta['description']}")
    
    def task_finished(self, task_id: str, summary: Dict[str, Any]) -> None:
        """Implementa Monitor.task_finished para compatibilidade."""
        status = summary.get('status', 'completed')
        if status == 'ok' or status == 'completed':
            symbol = Symbols.COMPLETED
            color = Colors.GREEN
        else:
            symbol = Symbols.FAILED  
            color = Colors.RED
            
        self._print_line(
            f"{symbol} Tarefa finalizada: {self._colorize(task_id, color)}"
        )
        
        if self.verbose and summary:
            for key, value in summary.items():
                if key != 'status':
                    self._print_line(f"  {key}: {value}")
    
    def unit_started(self, unit_id: str, meta: Dict[str, Any]) -> None:
        """Implementa Monitor.unit_started para compatibilidade."""
        if self.verbose:
            unit_type = meta.get('type', 'unit')
            self._print_line(f"  {Symbols.RUNNING} {unit_type}: {unit_id}")
    
    def unit_progress(self, unit_id: str, progress: float, message: Optional[str] = None) -> None:
        """Implementa Monitor.unit_progress para compatibilidade."""
        if not self.compact:
            progress_bar = self._create_progress_bar(progress, 20)
            line = f"    {progress_bar} {progress:5.1f}%"
            if message:
                line += f" - {message}"
            self._update_line(line)
    
    def unit_finished(self, unit_id: str, result: Dict[str, Any]) -> None:
        """Implementa Monitor.unit_finished para compatibilidade."""
        if self.verbose:
            status = result.get('status', 'completed')
            symbol = Symbols.COMPLETED if status == 'ok' else Symbols.FAILED
            self._print_line(f"  {symbol} Concluído: {unit_id}")
    
    def log(self, level: str, message: str, ctx: Optional[Dict[str, Any]] = None) -> None:
        """Implementa Monitor.log para compatibilidade."""
        # Mapear níveis para símbolos e cores
        level_map = {
            'debug': (Symbols.INFO, Colors.DIM),
            'info': (Symbols.INFO, Colors.CYAN),
            'warning': (Symbols.WARNING, Colors.YELLOW),
            'error': (Symbols.FAILED, Colors.RED),
            'critical': (Symbols.FAILED, Colors.BRIGHT_RED)
        }
        
        symbol, color = level_map.get(level.lower(), (Symbols.INFO, Colors.WHITE))
        
        # Mostrar logs de acordo com verbosidade
        if level.lower() in ['error', 'critical'] or self.verbose:
            formatted_message = self._colorize(message, color)
            self._print_line(f"{symbol} {formatted_message}")
            
            if ctx and self.verbose:
                for key, value in ctx.items():
                    self._print_line(f"    {key}: {value}")
    
    def error(self, unit_id: Optional[str], exc: Exception, ctx: Optional[Dict[str, Any]] = None) -> None:
        """Implementa Monitor.error para compatibilidade."""
        unit_info = f" (unit: {unit_id})" if unit_id else ""
        self._print_line(
            f"{Symbols.FAILED} {self._colorize('ERRO', Colors.BRIGHT_RED)}{unit_info}: {exc}"
        )
        
        if ctx and self.verbose:
            self._print_line("  Contexto:")
            for key, value in ctx.items():
                self._print_line(f"    {key}: {value}")
