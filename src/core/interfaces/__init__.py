"""
Interfaces padronizadas para CSP-BLFGA.

Este módulo define os protocolos e interfaces principais para garantir
compatibilidade e padronização entre componentes.
"""

from .algorithm import IAlgorithm, Result
from .console import IConsole, TaskSlot
from .console_impl import CursesConsoleAdapter, SimpleConsole, create_console
from .console_listener import ConsoleExecutionListener, create_console_listener

# Novas interfaces de execução
from .execution_console import IConsole as IExecutionConsole
from .execution_controller import ExecutionController, create_execution_controller
from .executor import IExecutor, TaskHandle, TaskStatus
from .factory import ExecutorFactory, create_auto_executor, create_executor
from .impl.execution_console_factory import create_console as create_execution_console
from .impl.execution_curses_console import CursesConsole as ExecutionCursesConsole
from .impl.execution_simple_console import SimpleConsole as ExecutionSimpleConsole
from .scheduler_logger import (
    SchedulerLogger,
    get_scheduler_logger,
    log_task_end,
    log_task_start,
)
from .task_result import TaskResult

__all__ = [
    "IAlgorithm",
    "Result",
    "IConsole",
    "TaskSlot",
    "SimpleConsole",
    "CursesConsoleAdapter",
    "create_console",
    "ConsoleExecutionListener",
    "create_console_listener",
    "ExecutionController",
    "create_execution_controller",
    "IExecutor",
    "TaskHandle",
    "TaskStatus",
    "create_executor",
    "create_auto_executor",
    "ExecutorFactory",
    "TaskResult",
    "SchedulerLogger",
    "get_scheduler_logger",
    "log_task_start",
    "log_task_end",
    "IExecutionConsole",
    "create_execution_console",
    "ExecutionSimpleConsole",
    "ExecutionCursesConsole",
]
