"""
Sistema de escalonamento de execução para CSP-BLFGA.

Este módulo implementa o ExecutionScheduler que gerencia a fila de tarefas
e controla a execução baseada em recursos disponíveis.
"""

from .executor import SchedulerExecutor, create_scheduler_executor
from .resource_monitor import (
    ProcessWatcher,
    ResourceChecker,
    ResourceMonitor,
    ResourceSnapshot,
    TaskEvent,
    TaskEventData,
)
from .scheduler import ExecutionScheduler, ScheduledTask, TaskState

__all__ = [
    "ExecutionScheduler",
    "TaskState",
    "ScheduledTask",
    "SchedulerExecutor",
    "create_scheduler_executor",
    "ResourceMonitor",
    "ResourceChecker",
    "ProcessWatcher",
    "ResourceSnapshot",
    "TaskEvent",
    "TaskEventData",
]
