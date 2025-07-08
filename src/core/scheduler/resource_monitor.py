"""
Monitoramento avançado de recursos do sistema.

Este módulo implementa um ResourceMonitor singleton leve que monitora
recursos do sistema e um watcher de processos para controlar execuções.
"""

import logging
import multiprocessing
import os
import queue
import signal
import threading
import time
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import psutil

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ResourceSnapshot:
    """Snapshot imutável dos recursos do sistema."""

    cpu_percent: float
    load_avg: float
    mem_available_mb: float
    timestamp: float

    def is_healthy(
        self, cpu_threshold: float = 85.0, mem_threshold: float = 300.0
    ) -> Tuple[bool, str]:
        """
        Verifica se os recursos estão saudáveis.

        Args:
            cpu_threshold: Limite de CPU em %
            mem_threshold: Limite de memória em MB

        Returns:
            Tuple[bool, str]: (True/False, razão)
        """
        if self.cpu_percent >= cpu_threshold:
            return False, f"CPU alto: {self.cpu_percent:.1f}% >= {cpu_threshold}%"

        if self.mem_available_mb <= mem_threshold:
            return (
                False,
                f"Memória baixa: {self.mem_available_mb:.1f}MB <= {mem_threshold}MB",
            )

        return True, "Recursos OK"


class TaskEvent(Enum):
    """Eventos de tarefas."""

    TASK_FINISHED = "task_finished"
    TASK_TIMEOUT = "task_timeout"
    TASK_ERROR = "task_error"


@dataclass
class TaskEventData:
    """Dados de evento de tarefa."""

    event: TaskEvent
    task_id: str
    process_id: int
    exit_code: Optional[int] = None
    error: Optional[Exception] = None
    traceback: Optional[str] = None
    timestamp: float = None

    def __post_init__(self):
        if self.timestamp is None:
            self.timestamp = time.time()


class ResourceMonitor:
    """
    Monitor singleton de recursos do sistema.

    Coleta métricas de CPU, memória e load average a cada MONITOR_INTERVAL.
    """

    _instance: Optional["ResourceMonitor"] = None
    _lock = threading.Lock()

    MONITOR_INTERVAL = 1.0  # 1 segundo

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if hasattr(self, "_initialized"):
            return

        self._initialized = True
        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._current_snapshot: Optional[ResourceSnapshot] = None
        self._snapshot_lock = threading.Lock()

        # Verificar se psutil está disponível
        try:
            import psutil

            self._psutil_available = True
        except ImportError:
            self._psutil_available = False
            logger.warning("psutil não disponível, monitoramento limitado")

    def start(self):
        """Inicia o monitoramento de recursos."""
        if self._running:
            return

        self._running = True
        self._thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self._thread.start()

        logger.info("ResourceMonitor iniciado")

    def stop(self):
        """Para o monitoramento de recursos."""
        if not self._running:
            return

        self._running = False
        if self._thread:
            self._thread.join()

        logger.info("ResourceMonitor parado")

    def snapshot(self) -> ResourceSnapshot:
        """
        Retorna snapshot atual dos recursos.

        Returns:
            ResourceSnapshot: Estado atual dos recursos
        """
        with self._snapshot_lock:
            if self._current_snapshot is None:
                # Coletar métricas imediatamente se não há snapshot
                self._current_snapshot = self._collect_metrics()
            return self._current_snapshot

    def _monitor_loop(self):
        """Loop principal de monitoramento."""
        while self._running:
            try:
                snapshot = self._collect_metrics()

                with self._snapshot_lock:
                    self._current_snapshot = snapshot

                time.sleep(self.MONITOR_INTERVAL)

            except Exception as e:
                logger.error("Erro no monitoramento: %s", e)
                time.sleep(self.MONITOR_INTERVAL)

    def _collect_metrics(self) -> ResourceSnapshot:
        """Coleta métricas atuais do sistema."""
        timestamp = time.time()

        # CPU percentage
        cpu_percent = 0.0
        if self._psutil_available:
            try:
                cpu_percent = psutil.cpu_percent(interval=None)
            except Exception:
                # Fallback para load average
                load_avg = os.getloadavg()[0]
                cpu_count = os.cpu_count()
                cpu_percent = min(100.0, (load_avg / cpu_count) * 100)
        else:
            # Usar load average como proxy
            load_avg = os.getloadavg()[0]
            cpu_count = os.cpu_count()
            cpu_percent = min(100.0, (load_avg / cpu_count) * 100)

        # Load average
        load_avg = os.getloadavg()[0]

        # Memória disponível
        mem_available_mb = 1024.0  # Default fallback
        if self._psutil_available:
            try:
                memory = psutil.virtual_memory()
                mem_available_mb = memory.available / (1024 * 1024)
            except Exception:
                pass

        return ResourceSnapshot(
            cpu_percent=cpu_percent,
            load_avg=load_avg,
            mem_available_mb=mem_available_mb,
            timestamp=timestamp,
        )


class ResourceChecker:
    """Verifica se é seguro iniciar novas tarefas."""

    def __init__(self, cpu_threshold: float = 85.0, mem_threshold: float = 300.0):
        self.cpu_threshold = cpu_threshold
        self.mem_threshold = mem_threshold
        self.monitor = ResourceMonitor()

        # Garantir que o monitor está rodando
        self.monitor.start()

    def ok_to_start(self) -> Tuple[bool, str]:
        """
        Verifica se é seguro iniciar uma nova tarefa.

        Returns:
            Tuple[bool, str]: (pode_iniciar, razão)
        """
        snapshot = self.monitor.snapshot()
        return snapshot.is_healthy(self.cpu_threshold, self.mem_threshold)

    def can_start_task(self) -> bool:
        """
        Método de compatibilidade.

        Returns:
            bool: True se pode iniciar nova tarefa
        """
        ok, _ = self.ok_to_start()
        return ok


@dataclass
class WatchedProcess:
    """Processo sendo monitorado."""

    task_id: str
    process: multiprocessing.Process
    started_at: float
    timeout: float
    error_queue: Optional[multiprocessing.Queue] = None


class ProcessWatcher:
    """
    Watcher de processos ativos.

    Monitora processos filhos e gera eventos quando terminam, falham ou excedem timeout.
    """

    def __init__(self, check_interval: float = 1.0):
        self.check_interval = check_interval
        self.watched_processes: Dict[str, WatchedProcess] = {}
        self.event_queue = queue.Queue()

        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._lock = threading.Lock()

        # Context para multiprocessing
        self._mp_context = multiprocessing.get_context("spawn")

    def start(self):
        """Inicia o watcher de processos."""
        if self._running:
            return

        self._running = True
        self._thread = threading.Thread(target=self._watch_loop, daemon=True)
        self._thread.start()

        logger.info("ProcessWatcher iniciado")

    def stop(self):
        """Para o watcher de processos."""
        if not self._running:
            return

        self._running = False
        if self._thread:
            self._thread.join()

        # Terminar processos ainda ativos
        with self._lock:
            for watched in self.watched_processes.values():
                if watched.process.is_alive():
                    watched.process.terminate()
                    watched.process.join(timeout=5)
                    if watched.process.is_alive():
                        watched.process.kill()

        logger.info("ProcessWatcher parado")

    def watch_process(
        self,
        task_id: str,
        process: multiprocessing.Process,
        timeout: float,
        error_queue: Optional[multiprocessing.Queue] = None,
    ):
        """
        Adiciona processo para monitoramento.

        Args:
            task_id: ID da tarefa
            process: Processo a ser monitorado
            timeout: Timeout em segundos
            error_queue: Queue para capturar erros
        """
        with self._lock:
            self.watched_processes[task_id] = WatchedProcess(
                task_id=task_id,
                process=process,
                started_at=time.time(),
                timeout=timeout,
                error_queue=error_queue,
            )

        logger.debug("Processo %s adicionado ao watcher", task_id)

    def unwatch_process(self, task_id: str):
        """Remove processo do monitoramento."""
        with self._lock:
            self.watched_processes.pop(task_id, None)

        logger.debug("Processo %s removido do watcher", task_id)

    def get_event(self, timeout: Optional[float] = None) -> Optional[TaskEventData]:
        """
        Obtém próximo evento da fila.

        Args:
            timeout: Timeout em segundos

        Returns:
            TaskEventData ou None se timeout
        """
        try:
            return self.event_queue.get(timeout=timeout)
        except queue.Empty:
            return None

    def _watch_loop(self):
        """Loop principal de monitoramento."""
        while self._running:
            try:
                current_time = time.time()
                completed_tasks = []

                with self._lock:
                    for task_id, watched in self.watched_processes.items():
                        # Verificar se processo terminou
                        if not watched.process.is_alive():
                            exit_code = watched.process.exitcode
                            error = None
                            traceback_str = None

                            # Tentar capturar erro da queue
                            if watched.error_queue:
                                try:
                                    error_data = watched.error_queue.get_nowait()
                                    if isinstance(error_data, dict):
                                        error = error_data.get("error")
                                        traceback_str = error_data.get("traceback")
                                except queue.Empty:
                                    pass

                            # Determinar tipo de evento
                            if exit_code == 0:
                                event = TaskEvent.TASK_FINISHED
                            else:
                                event = TaskEvent.TASK_ERROR

                            event_data = TaskEventData(
                                event=event,
                                task_id=task_id,
                                process_id=watched.process.pid,
                                exit_code=exit_code,
                                error=error,
                                traceback=traceback_str,
                            )

                            self.event_queue.put(event_data)
                            completed_tasks.append(task_id)

                            logger.debug(
                                f"Processo {task_id} terminou com código {exit_code}"
                            )

                        # Verificar timeout
                        elif current_time - watched.started_at > watched.timeout:
                            logger.warning(
                                f"Processo {task_id} excedeu timeout de {watched.timeout}s"
                            )

                            # Terminar processo
                            try:
                                watched.process.terminate()
                                watched.process.join(timeout=5)
                                if watched.process.is_alive():
                                    watched.process.kill()
                            except Exception as e:
                                logger.error(
                                    f"Erro ao terminar processo {task_id}: {e}"
                                )

                            # Gerar evento de timeout
                            event_data = TaskEventData(
                                event=TaskEvent.TASK_TIMEOUT,
                                task_id=task_id,
                                process_id=watched.process.pid,
                                exit_code=None,
                            )

                            self.event_queue.put(event_data)
                            completed_tasks.append(task_id)

                # Remover processos completados
                with self._lock:
                    for task_id in completed_tasks:
                        self.watched_processes.pop(task_id, None)

                time.sleep(self.check_interval)

            except Exception as e:
                logger.error("Erro no watcher loop: %s", e)
                time.sleep(self.check_interval)

    def get_active_count(self) -> int:
        """Retorna número de processos ativos."""
        with self._lock:
            return len(self.watched_processes)

    def get_active_tasks(self) -> List[str]:
        """Retorna lista de IDs de tarefas ativas."""
        with self._lock:
            return list(self.watched_processes.keys())


# Função para executar tarefa em processo separado
def _execute_task_in_process(func, args, kwargs, error_queue):
    """
    Executa tarefa em processo separado com captura de erros.

    Args:
        func: Função a ser executada
        args: Argumentos posicionais
        kwargs: Argumentos nomeados
        error_queue: Queue para reportar erros
    """
    try:
        return func(*args, **kwargs)
    except Exception as e:
        import traceback

        error_data = {"error": e, "traceback": traceback.format_exc()}
        error_queue.put(error_data)
        raise
