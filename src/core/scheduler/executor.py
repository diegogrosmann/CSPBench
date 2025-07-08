"""
Wrapper que implementa IExecutor usando ExecutionScheduler.

Este módulo adapta o ExecutionScheduler para a interface IExecutor,
permitindo uso transparente do novo sistema de escalonamento com TaskResult.
"""

import logging
import time
from typing import Any, Union

try:
    from ..interfaces import (
        IAlgorithm,
        IExecutor,
        Result,
        TaskHandle,
        TaskResult,
        TaskStatus,
    )
    from .scheduler import ExecutionScheduler, TaskState
except ImportError:
    # Para execução direta (testes)
    import os
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    from src.core.interfaces import (
        IAlgorithm,
        IExecutor,
        Result,
        TaskHandle,
        TaskResult,
        TaskStatus,
    )
    from src.core.scheduler.scheduler import ExecutionScheduler, TaskState

logger = logging.getLogger(__name__)


class SchedulerExecutor(IExecutor):
    """
    Wrapper que implementa IExecutor usando ExecutionScheduler.

    Permite usar o novo sistema de escalonamento com a interface
    padrão IExecutor, mantendo compatibilidade com código existente.
    """

    def __init__(self, start_delay: float = 2.0, timeout: int = 300):
        """
        Inicializa o executor baseado em scheduler.

        Args:
            start_delay: Delay entre início de tarefas em segundos
            timeout: Timeout por tarefa em segundos
        """
        self.scheduler = ExecutionScheduler(start_delay=start_delay)
        self.timeout = timeout
        self._started = False

        logger.info(
            f"SchedulerExecutor inicializado com delay={start_delay}s, timeout={timeout}s"
        )

    def __enter__(self):
        """Context manager entry."""
        self.scheduler.start()
        self._started = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.scheduler.stop()
        self._started = False

    def submit(
        self, algorithm_instance: IAlgorithm, *, priority: int = 0
    ) -> TaskHandle:
        """
        Submete um algoritmo para execução.

        Args:
            algorithm_instance: Instância do algoritmo
            priority: Prioridade da tarefa (ignorada - FIFO absoluta)

        Returns:
            TaskHandle: Handle da tarefa
        """
        if not self._started:
            self.scheduler.start()
            self._started = True

        # Extrair nome do algoritmo para logging
        algorithm_name = getattr(
            algorithm_instance,
            "name",
            getattr(algorithm_instance.__class__, "__name__", "unknown"),
        )

        # Criar wrapper para executar o algoritmo
        def execute_algorithm():
            try:
                return algorithm_instance.run()
            except Exception as e:
                logger.error("Erro na execução do algoritmo %s: %s", algorithm_name, e)
                raise

        # Submeter com metadados
        task_id = self.scheduler.submit(
            execute_algorithm,
            algorithm_name=algorithm_name,
            timeout=self.timeout,
            metadata={
                "algorithm_type": type(algorithm_instance).__name__,
                "is_deterministic": getattr(
                    algorithm_instance, "is_deterministic", False
                ),
            },
        )

        return TaskHandle(task_id)

    def poll(self, handle: TaskHandle) -> TaskStatus:
        """
        Verifica o status de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            TaskStatus: Status atual da tarefa
        """
        task = self.scheduler.get_task(handle.task_id)

        if not task:
            return TaskStatus.ERROR

        if task.state == TaskState.QUEUED:
            return TaskStatus.QUEUED  # Mudança: usar QUEUED para tarefas em fila
        elif task.state == TaskState.RUNNING:
            return TaskStatus.RUNNING  # Só para tarefas realmente executando
        elif task.state == TaskState.COMPLETED:
            return TaskStatus.DONE
        elif task.state == TaskState.FAILED:
            return TaskStatus.ERROR
        elif task.state == TaskState.CANCELLED:
            return TaskStatus.ERROR
        else:
            return TaskStatus.QUEUED  # Mudança: usar QUEUED como padrão

    def result(self, handle: TaskHandle) -> Union[Result, Exception]:
        """
        Obtém o resultado de uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            Result ou Exception: Resultado da execução
        """
        try:
            task_result = self.scheduler.get_task_result(
                handle.task_id, timeout=self.timeout
            )

            # Converter TaskResult para Result (interface legada)
            if isinstance(task_result, TaskResult):
                if task_result.success:
                    # Incluir tempo no metadata
                    result_metadata = task_result.metadata.copy()
                    result_metadata["execution_time"] = task_result.time
                    result_metadata["tempo"] = task_result.time  # Compatibilidade
                    result_metadata["time"] = task_result.time  # Compatibilidade

                    return Result(
                        center=task_result.center or "",
                        distance=task_result.distance or 0.0,
                        metadata=result_metadata,
                    )
                else:
                    # Retornar exceção para falhas
                    return Exception(task_result.error or "Falha na execução")
            else:
                # Resultado legado, retornar como está
                return task_result

        except Exception as e:
            return e

    def cancel(self, handle: TaskHandle) -> bool:
        """
        Cancela uma tarefa.

        Args:
            handle: Handle da tarefa

        Returns:
            bool: True se a tarefa foi cancelada
        """
        return self.scheduler.cancel_task(handle.task_id)

    def shutdown(self, wait: bool = True) -> None:
        """
        Encerra o executor.

        Args:
            wait: Se True, aguarda tarefas terminarem
        """
        if self._started:
            if wait:
                # Aguardar conclusão de todas as tarefas
                try:
                    self.scheduler.wait_for_completion(timeout=self.timeout)
                except Exception as e:
                    logger.warning("Timeout aguardando conclusão das tarefas: %s", e)

            self.scheduler.stop()
            self._started = False

    def get_stats(self) -> dict:
        """
        Retorna estatísticas do scheduler.

        Returns:
            dict: Estatísticas atuais
        """
        return self.scheduler.get_status()

    def get_active_tasks_count(self) -> int:
        """
        Retorna o número de tarefas ativas.

        Returns:
            int: Número de tarefas em execução
        """
        stats = self.scheduler.get_status()
        return stats.get("active_tasks", 0)

    def get_queue_size(self) -> int:
        """
        Retorna o tamanho da fila de tarefas.

        Returns:
            int: Número de tarefas na fila
        """
        stats = self.scheduler.get_status()
        return stats.get("queue_size", 0)


def create_scheduler_executor(
    start_delay: float = 2.0, timeout: int = 300
) -> SchedulerExecutor:
    """
    Cria um executor baseado no novo sistema de escalonamento.

    Args:
        start_delay: Delay entre início de tarefas em segundos
        timeout: Timeout por tarefa em segundos

    Returns:
        SchedulerExecutor: Executor com escalonamento
    """
    return SchedulerExecutor(start_delay=start_delay, timeout=timeout)
