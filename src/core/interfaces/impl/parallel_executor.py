"""
Executor paralelo que implementa IExecutor.

Wrapper que adapta o SchedulerExecutor para compatibilidade com código legado.
"""

import logging
import threading
import time
from typing import Dict, Union

from ..algorithm import IAlgorithm, Result
from ..executor import IExecutor, TaskHandle, TaskStatus, create_task_handle

logger = logging.getLogger(__name__)


class ParallelExecutor(IExecutor):
    """
    Executor paralelo que implementa IExecutor.

    Wrapper que adapta o SchedulerExecutor para compatibilidade com código legado.
    """

    def __init__(self, max_workers: int | None = None, timeout: int = 300):
        """
        Inicializa o executor paralelo.

        Args:
            max_workers: Número máximo de workers (ignorado - usa SchedulerExecutor)
            timeout: Timeout por tarefa em segundos
        """
        self.max_workers = max_workers or 4
        self.timeout = timeout

        # Importação tardia para evitar importação circular
        from ...scheduler.executor import SchedulerExecutor

        self._scheduler_executor = SchedulerExecutor(start_delay=0.1, timeout=timeout)
        self._started = False
        self._shutdown = False
        self._lock = threading.Lock()
        self._tasks = {}

        logger.info(
            f"ParallelExecutor inicializado com timeout={timeout}s (usando SchedulerExecutor)"
        )

    def __enter__(self):
        """Context manager entry."""
        self._scheduler_executor.__enter__()
        self._started = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self._scheduler_executor.__exit__(exc_type, exc_val, exc_tb)
        self._started = False

    def submit(
        self, algorithm_instance: IAlgorithm, *, priority: int = 0
    ) -> TaskHandle:
        """
        Submete um algoritmo para execução paralela.

        Args:
            algorithm_instance: Instância do algoritmo
            priority: Prioridade da tarefa (menor valor = maior prioridade)

        Returns:
            TaskHandle: Handle da tarefa
        """
        if self._shutdown:
            raise RuntimeError("Executor foi encerrado")

        # Usar o SchedulerExecutor diretamente
        handle = self._scheduler_executor.submit(algorithm_instance, priority=priority)

        logger.info(
            f"Tarefa {handle.task_id} submetida para execução paralela (prioridade: {priority})"
        )
        return handle

    def poll(self, handle: TaskHandle) -> TaskStatus:
        """
        Verifica o status de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            TaskStatus: Status atual da tarefa
        """
        return self._scheduler_executor.poll(handle)

    def result(self, handle: TaskHandle) -> Union[Result, Exception]:
        """
        Obtém o resultado de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            Result ou Exception: Resultado da execução

        Raises:
            RuntimeError: Se a tarefa ainda não foi finalizada
        """
        return self._scheduler_executor.result(handle)

    def cancel(self, handle: TaskHandle) -> bool:
        """
        Cancela uma tarefa em execução.

        Args:
            handle: Handle da tarefa

        Returns:
            bool: True se a tarefa foi cancelada com sucesso
        """
        return self._scheduler_executor.cancel(handle)

    def shutdown(self, wait: bool = True) -> None:
        """
        Encerra o executor e limpa recursos.

        Args:
            wait: Se True, aguarda todas as tarefas terminarem
        """
        self._shutdown = True
        self._scheduler_executor.shutdown(wait=wait)
        logger.info("ParallelExecutor encerrado")

    def _monitor_task(self, handle: TaskHandle) -> None:
        """
        Monitora uma tarefa em thread separada.

        Args:
            handle: Handle da tarefa
        """

        def monitor():
            while not self._shutdown:
                status = self.poll(handle)
                if status != TaskStatus.RUNNING:
                    break
                time.sleep(0.1)  # Verificar a cada 100ms

        thread = threading.Thread(target=monitor, daemon=True)
        thread.start()

    def get_active_tasks(self) -> int:
        """
        Retorna o número de tarefas ativas.

        Returns:
            int: Número de tarefas em execução
        """
        with self._lock:
            return sum(
                1
                for task in self._tasks.values()
                if task["handle"]._status == TaskStatus.RUNNING
            )

    def get_task_info(self, handle: TaskHandle) -> dict:
        """
        Obtém informações detalhadas sobre uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            dict: Informações da tarefa
        """
        with self._lock:
            if handle.task_id not in self._tasks:
                return {}

            task = self._tasks[handle.task_id]
            return {
                "task_id": handle.task_id,
                "status": handle._status.value,
                "priority": task["priority"],
                "submitted_at": task["submitted_at"],
                "started_at": task["started_at"],
                "finished_at": task["finished_at"],
                "duration": (
                    task["finished_at"] - task["started_at"]
                    if task["started_at"] and task["finished_at"]
                    else None
                ),
            }

    def wait_for_completion(
        self, handles: list[TaskHandle], timeout: int | None = None
    ) -> list[TaskHandle]:
        """
        Aguarda múltiplas tarefas terminarem.

        Args:
            handles: Lista de handles para aguardar
            timeout: Timeout total em segundos

        Returns:
            list[TaskHandle]: Lista de handles que terminaram
        """
        start_time = time.time()
        completed = []

        while handles and not self._shutdown:
            for handle in handles[:]:  # Cópia para modificar durante iteração
                status = self.poll(handle)
                if status != TaskStatus.RUNNING:
                    completed.append(handle)
                    handles.remove(handle)

            # Verificar timeout
            if timeout and (time.time() - start_time) > timeout:
                break

            time.sleep(0.1)  # Verificar a cada 100ms

        return completed
