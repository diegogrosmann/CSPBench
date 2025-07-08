"""
Scheduler avançado com fila FIFO absoluta e controle de recursos.

Este módulo implementa um scheduler que:
- Mantém uma fila FIFO absoluta para todas as tarefas
- Controla automaticamente o número de workers baseado nos recursos
- Verifica recursos antes de iniciar tarefas
- Aplica delay entre inícios de tarefas
- Garante ordem de execução rigorosa
- Monitora processos filhos com timeout
- Usa TaskResult padronizado e logging estruturado
"""

import logging
import multiprocessing
import os
import queue
import threading
import time
import traceback
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Set
from uuid import uuid4

from ..interfaces.scheduler_logger import get_scheduler_logger

# Importar TaskResult e logging estruturado
from ..interfaces.task_result import TaskResult
from .resource_monitor import (
    ProcessWatcher,
    ResourceChecker,
    ResourceMonitor,
    TaskEvent,
    TaskEventData,
)

logger = logging.getLogger(__name__)
scheduler_logger = get_scheduler_logger()

logger = logging.getLogger(__name__)


class TaskState(Enum):
    """Estados possíveis de uma tarefa."""

    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class ScheduledTask:
    """Informações sobre uma tarefa no scheduler."""

    task_id: str
    func: Callable
    args: tuple
    kwargs: dict
    future: Future
    state: TaskState = TaskState.QUEUED
    submitted_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    finished_at: Optional[float] = None
    worker_id: Optional[str] = None
    error: Optional[Exception] = None
    process: Optional[multiprocessing.Process] = None
    algorithm_name: Optional[str] = None
    timeout_limit: float = 300.0
    metadata: Dict[str, Any] = field(default_factory=dict)


class ExecutionScheduler:
    """
    Scheduler avançado com fila FIFO absoluta e controle de recursos.

    Características:
    - Fila FIFO rigorosamente respeitada
    - Controle automático de workers baseado em recursos
    - Verificação de recursos antes de iniciar tarefas
    - Delay configurável entre inícios
    - Mantém pelo menos 1 tarefa ativa se houver tarefas pendentes
    """

    def __init__(
        self,
        start_delay: float = 2.0,
        cpu_threshold: float = 85.0,
        mem_threshold: float = 300.0,
        min_active_tasks: int = 1,
        timeout: float = 300.0,
    ):
        """
        Inicializa o scheduler.

        Args:
            start_delay: Delay em segundos entre inícios de tarefas
            cpu_threshold: Limite de CPU em %
            mem_threshold: Limite de memória em MB
            min_active_tasks: Número mínimo de tarefas ativas
            timeout: Timeout padrão por tarefa em segundos
        """
        self.start_delay = start_delay
        self.min_active_tasks = min_active_tasks
        self.timeout = timeout

        # Usar novo sistema de monitoramento
        self.resource_checker = ResourceChecker(cpu_threshold, mem_threshold)
        self.process_watcher = ProcessWatcher()
        self.resource_monitor = ResourceMonitor()

        # Fila FIFO absoluta
        self.task_queue = queue.Queue()
        self.scheduled_tasks: Dict[str, ScheduledTask] = {}
        self.active_tasks: Set[str] = set()
        self.completed_tasks: Set[str] = set()

        # Controle de workers
        self.max_workers = os.cpu_count() or 4  # Fallback para 4
        self.executor = ThreadPoolExecutor(max_workers=self.max_workers)

        # Controle de estado
        self.running = False
        self.scheduler_thread: Optional[threading.Thread] = None
        self.event_handler_thread: Optional[threading.Thread] = None
        self.lock = threading.Lock()

        # Métricas
        self.total_submitted = 0
        self.total_completed = 0
        self.last_start_time = 0.0

        logger.info("Scheduler iniciado com %s workers máximos", self.max_workers)

    def submit(
        self,
        func: Callable,
        *args,
        algorithm_name: Optional[str] = None,
        timeout: Optional[float] = None,
        metadata: Optional[Dict[str, Any]] = None,
        **kwargs,
    ) -> str:
        """
        Submete uma tarefa para execução.

        Args:
            func: Função a ser executada
            *args: Argumentos posicionais
            algorithm_name: Nome do algoritmo (para logging)
            timeout: Timeout específico para esta tarefa
            metadata: Metadados adicionais
            **kwargs: Argumentos nomeados

        Returns:
            str: ID da tarefa
        """
        task_id = str(uuid4())
        future = Future()

        task = ScheduledTask(
            task_id=task_id,
            func=func,
            args=args,
            kwargs=kwargs,
            future=future,
            algorithm_name=algorithm_name or getattr(func, "__name__", "unknown"),
            timeout_limit=timeout or self.timeout,
            metadata=metadata or {},
        )

        with self.lock:
            self.scheduled_tasks[task_id] = task
            self.task_queue.put(task_id)
            self.total_submitted += 1

        logger.debug("Tarefa %s adicionada à fila", task_id)

        # Iniciar scheduler se não estiver rodando
        if not self.running:
            self.start()

        return task_id

    def start(self):
        """Inicia o scheduler e todos os componentes de monitoramento."""
        if self.running:
            return

        self.running = True

        # Iniciar componentes de monitoramento
        self.resource_monitor.start()
        self.process_watcher.start()

        # Iniciar threads do scheduler
        self.scheduler_thread = threading.Thread(
            target=self._scheduler_loop, daemon=True
        )
        self.scheduler_thread.start()

        self.event_handler_thread = threading.Thread(
            target=self._event_handler_loop, daemon=True
        )
        self.event_handler_thread.start()

        logger.info("ExecutionScheduler iniciado")

    def stop(self):
        """Para o scheduler e aguarda conclusão das tarefas ativas."""
        self.running = False

        if self.scheduler_thread:
            self.scheduler_thread.join()

        if self.event_handler_thread:
            self.event_handler_thread.join()

        # Parar componentes de monitoramento
        self.process_watcher.stop()
        self.resource_monitor.stop()

        self.executor.shutdown(wait=True)
        logger.info("ExecutionScheduler parado")

    def _scheduler_loop(self):
        """Loop principal do scheduler."""
        while self.running:
            try:
                # Verificar se pode iniciar nova tarefa
                should_start = self._should_start_new_task()

                if should_start and not self.task_queue.empty():
                    self._start_next_task()

                # Limpar tarefas completadas
                self._cleanup_completed_tasks()

                # Aguardar antes da próxima verificação
                time.sleep(0.1)

            except Exception as e:
                logger.error("Erro no loop do scheduler: %s", e)
                time.sleep(1.0)

    def _should_start_new_task(self) -> bool:
        """
        Determina se deve iniciar uma nova tarefa.

        Returns:
            bool: True se deve iniciar nova tarefa
        """
        with self.lock:
            active_count = len(self.active_tasks)

            # Sempre manter pelo menos min_active_tasks se houver tarefas
            if active_count < self.min_active_tasks and not self.task_queue.empty():
                return True

            # Verificar se já atingiu o máximo de workers
            if active_count >= self.max_workers:
                return False

            # Verificar delay entre inícios
            current_time = time.time()
            if current_time - self.last_start_time < self.start_delay:
                return False

            # Verificar recursos do sistema
            ok_to_start, reason = self.resource_checker.ok_to_start()
            if not ok_to_start:
                logger.debug("Recursos não permitem nova tarefa: %s", reason)
                return False

            return True

    def _start_next_task(self):
        """Inicia a próxima tarefa da fila."""
        try:
            task_id = self.task_queue.get_nowait()
            task = self.scheduled_tasks[task_id]

            with self.lock:
                self.active_tasks.add(task_id)
                task.state = TaskState.RUNNING
                task.started_at = time.time()
                self.last_start_time = task.started_at

            logger.info("Iniciando tarefa %s", task_id)

            # Executar tarefa
            future = self.executor.submit(self._execute_task, task)
            task.future = future

        except queue.Empty:
            pass
        except Exception as e:
            logger.error("Erro ao iniciar tarefa: %s", e)

    def _execute_task(self, task: ScheduledTask) -> TaskResult:
        """
        Executa uma tarefa específica com logging estruturado.

        Args:
            task: Tarefa a ser executada

        Returns:
            TaskResult: Resultado padronizado da tarefa
        """
        start_time = time.time()

        # Log início da tarefa
        scheduler_logger.log_task_start(
            task.task_id, task.algorithm_name or "unknown", task.metadata
        )

        try:
            # Executar a função
            result = task.func(*task.args, **task.kwargs)
            execution_time = time.time() - start_time

            # Processar resultado baseado no tipo retornado
            if isinstance(result, dict) and "center" in result and "distance" in result:
                # Resultado de algoritmo CSP
                task_result = TaskResult.success_result(
                    distance=float(result["distance"]),
                    center=str(result["center"]),
                    time=execution_time,
                    metadata={**task.metadata, "original_result": result},
                )
            elif isinstance(result, tuple) and len(result) >= 2:
                # Resultado em formato (center, distance, metadata)
                center, distance = result[0], result[1]
                metadata = result[2] if len(result) > 2 else {}
                task_result = TaskResult.success_result(
                    distance=float(distance),
                    center=str(center),
                    time=execution_time,
                    metadata={**task.metadata, **metadata},
                )
            else:
                # Resultado genérico
                task_result = TaskResult.success_result(
                    distance=0.0,
                    center=str(result),
                    time=execution_time,
                    metadata={
                        **task.metadata,
                        "result_type": type(result).__name__,
                        "raw_result": str(result),
                    },
                )

            # Marcar como concluída
            with self.lock:
                task.state = TaskState.COMPLETED
                task.finished_at = time.time()
                self.active_tasks.discard(task.task_id)
                self.completed_tasks.add(task.task_id)
                self.total_completed += 1

            # Log fim da tarefa
            scheduler_logger.log_task_end(
                task.task_id,
                task.algorithm_name or "unknown",
                task_result,
                {"worker_id": task.worker_id},
            )

            logger.info(
                "Tarefa %s concluída: %s", task.task_id, task_result.get_summary()
            )
            return task_result

        except Exception as e:
            execution_time = time.time() - start_time
            error_message = str(e)
            error_traceback = traceback.format_exc()

            # Verificar se é timeout
            if (
                "timeout" in error_message.lower()
                or execution_time >= task.timeout_limit
            ):
                task_result = TaskResult.timeout_result(
                    time=execution_time,
                    timeout_limit=task.timeout_limit,
                    metadata={**task.metadata, "error_type": type(e).__name__},
                )
            else:
                task_result = TaskResult.failure_result(
                    error=error_message,
                    time=execution_time,
                    traceback=error_traceback,
                    metadata={**task.metadata, "error_type": type(e).__name__},
                )

            logger.error(
                "Erro na tarefa %s: %s", task.task_id, task_result.get_summary()
            )

            with self.lock:
                task.state = TaskState.FAILED
                task.error = e
                task.finished_at = time.time()
                self.active_tasks.discard(task.task_id)
                self.completed_tasks.add(task.task_id)

            # Log fim da tarefa com erro
            scheduler_logger.log_task_end(
                task.task_id,
                task.algorithm_name or "unknown",
                task_result,
                {"worker_id": task.worker_id},
            )

            return task_result

    def _cleanup_completed_tasks(self):
        """Remove tarefas antigas da memória."""
        with self.lock:
            # Manter apenas as últimas 100 tarefas completadas
            if len(self.completed_tasks) > 100:
                # Remover as 50 mais antigas
                tasks_to_remove = []
                for task_id in list(self.completed_tasks):
                    if len(tasks_to_remove) >= 50:
                        break
                    tasks_to_remove.append(task_id)

                for task_id in tasks_to_remove:
                    self.completed_tasks.remove(task_id)
                    self.scheduled_tasks.pop(task_id, None)

    def get_task(self, task_id: str) -> Optional[ScheduledTask]:
        """
        Retorna informações sobre uma tarefa.

        Args:
            task_id: ID da tarefa

        Returns:
            Optional[ScheduledTask]: Informações da tarefa ou None se não encontrada
        """
        return self.scheduled_tasks.get(task_id)

    def get_task_result(self, task_id: str, timeout: Optional[float] = None) -> Any:
        """
        Retorna o resultado de uma tarefa.

        Args:
            task_id: ID da tarefa
            timeout: Timeout em segundos

        Returns:
            Any: Resultado da tarefa

        Raises:
            KeyError: Se a tarefa não existe
            TimeoutError: Se o timeout for atingido
        """
        task = self.scheduled_tasks.get(task_id)
        if not task:
            raise KeyError(f"Tarefa {task_id} não encontrada")

        return task.future.result(timeout=timeout)

    def cancel_task(self, task_id: str) -> bool:
        """
        Cancela uma tarefa.

        Args:
            task_id: ID da tarefa

        Returns:
            bool: True se a tarefa foi cancelada
        """
        task = self.scheduled_tasks.get(task_id)
        if not task:
            return False

        # Tentar cancelar o Future
        cancelled = task.future.cancel()

        if cancelled:
            with self.lock:
                task.state = TaskState.CANCELLED
                self.active_tasks.discard(task_id)
                self.completed_tasks.add(task_id)

            logger.info("Tarefa %s cancelada", task_id)

        return cancelled

    def get_status(self) -> Dict[str, Any]:
        """
        Retorna o status atual do scheduler.

        Returns:
            Dict[str, Any]: Status do scheduler
        """
        with self.lock:
            # Obter snapshot de recursos
            resource_snapshot = self.resource_monitor.snapshot()
            ok_to_start, reason = self.resource_checker.ok_to_start()

            return {
                "running": self.running,
                "total_submitted": self.total_submitted,
                "total_completed": self.total_completed,
                "queue_size": self.task_queue.qsize(),
                "active_tasks": len(self.active_tasks),
                "max_workers": self.max_workers,
                "start_delay": self.start_delay,
                "can_start_task": ok_to_start,
                "resource_status": reason,
                "resources": {
                    "cpu_percent": resource_snapshot.cpu_percent,
                    "load_avg": resource_snapshot.load_avg,
                    "mem_available_mb": resource_snapshot.mem_available_mb,
                    "timestamp": resource_snapshot.timestamp,
                },
                "process_watcher": {
                    "active_processes": self.process_watcher.get_active_count(),
                    "monitored_tasks": self.process_watcher.get_active_tasks(),
                },
            }

    def get_active_tasks(self) -> List[str]:
        """
        Retorna lista de IDs das tarefas ativas.

        Returns:
            List[str]: IDs das tarefas ativas
        """
        with self.lock:
            return list(self.active_tasks)

    def get_queued_tasks(self) -> List[str]:
        """
        Retorna lista de IDs das tarefas na fila.

        Returns:
            List[str]: IDs das tarefas na fila
        """
        with self.lock:
            queued = []
            for task_id, task in self.scheduled_tasks.items():
                if task.state == TaskState.QUEUED:
                    queued.append(task_id)
            return queued

    def wait_for_completion(self, timeout: Optional[float] = None):
        """
        Aguarda até que todas as tarefas sejam concluídas.

        Args:
            timeout: Timeout em segundos

        Raises:
            TimeoutError: Se o timeout for atingido
        """
        start_time = time.time()

        while True:
            with self.lock:
                if self.task_queue.empty() and len(self.active_tasks) == 0:
                    break

            if timeout and (time.time() - start_time) > timeout:
                raise TimeoutError("Timeout aguardando conclusão das tarefas")

            time.sleep(0.1)

    def _event_handler_loop(self):
        """Loop para processar eventos do ProcessWatcher."""
        while self.running:
            try:
                # Obter próximo evento (timeout de 1 segundo)
                event_data = self.process_watcher.get_event(timeout=1.0)

                if event_data is None:
                    continue

                # Processar evento
                self._handle_task_event(event_data)

            except Exception as e:
                logger.error("Erro no event handler loop: %s", e)
                time.sleep(1.0)

    def _handle_task_event(self, event_data: TaskEventData):
        """
        Processa evento de tarefa.

        Args:
            event_data: Dados do evento
        """
        task_id = event_data.task_id
        task = self.scheduled_tasks.get(task_id)

        if not task:
            logger.warning("Evento para tarefa desconhecida: %s", task_id)
            return

        with self.lock:
            if event_data.event == TaskEvent.TASK_FINISHED:
                # Tarefa completou com sucesso
                task.state = TaskState.COMPLETED
                task.finished_at = event_data.timestamp
                self.active_tasks.discard(task_id)
                self.completed_tasks.add(task_id)
                self.total_completed += 1

                logger.info("Tarefa %s completada via processo", task_id)

            elif event_data.event == TaskEvent.TASK_TIMEOUT:
                # Tarefa excedeu timeout
                task.state = TaskState.FAILED
                task.error = TimeoutError(f"Tarefa excedeu timeout de {self.timeout}s")
                task.finished_at = event_data.timestamp
                self.active_tasks.discard(task_id)
                self.completed_tasks.add(task_id)

                # Marcar Future como erro
                if not task.future.done():
                    task.future.set_exception(task.error)

                logger.warning("Tarefa %s excedeu timeout", task_id)

            elif event_data.event == TaskEvent.TASK_ERROR:
                # Tarefa falhou
                task.state = TaskState.FAILED
                task.error = event_data.error or Exception(
                    f"Processo terminou com código {event_data.exit_code}"
                )
                task.finished_at = event_data.timestamp
                self.active_tasks.discard(task_id)
                self.completed_tasks.add(task_id)

                # Marcar Future como erro
                if not task.future.done():
                    task.future.set_exception(task.error)

                logger.error("Tarefa %s falhou: %s", task_id, task.error)

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.stop()
