"""
Factory para criação de executores baseado em configuração.

Facilita a criação de executores sequenciais ou paralelos baseado
em parâmetros ou configuração do sistema.
"""

import logging
import multiprocessing
from typing import Optional

from .executor import IExecutor
from .impl import ParallelExecutor

logger = logging.getLogger(__name__)


class ExecutorFactory:
    """
    Factory para criação de executores de algoritmos.

    Permite criar executores sequenciais ou paralelos baseado em
    configuração ou detecção automática de recursos.
    """

    @staticmethod
    def create_executor(
        executor_type: str = "auto",
        max_workers: Optional[int] = None,
        timeout: int = 300,
        **kwargs,
    ) -> IExecutor:
        """
        Cria um executor baseado nos parâmetros fornecidos.

        Args:
            executor_type: Tipo de executor ("sequential", "parallel", "scheduler", "auto")
            max_workers: Número máximo de workers (para executor paralelo)
            timeout: Timeout por tarefa em segundos
            **kwargs: Argumentos adicionais para o executor

        Returns:
            IExecutor: Instância do executor criado
        """
        if executor_type == "sequential":
            logger.info(
                "Criando executor sequencial (usando ParallelExecutor com 1 worker)"
            )
            return ParallelExecutor(max_workers=1, timeout=timeout, **kwargs)

        elif executor_type == "parallel":
            logger.info("Criando executor paralelo com %s workers", max_workers)
            return ParallelExecutor(max_workers=max_workers, timeout=timeout, **kwargs)

        elif executor_type == "scheduler":
            logger.info("Criando executor com scheduler avançado")
            from ..scheduler.executor import SchedulerExecutor

            start_delay = kwargs.pop(
                "start_delay", 0.1
            )  # Reduzir delay para paralelismo
            return SchedulerExecutor(start_delay=start_delay, timeout=timeout)

        elif executor_type == "auto":
            # Decidir automaticamente baseado no sistema
            cpu_count = multiprocessing.cpu_count()

            if cpu_count >= 4:
                # Para sistemas com 4+ CPUs, usar scheduler avançado
                logger.info(
                    f"Detectado {cpu_count} CPUs, usando executor com scheduler"
                )
                from ..scheduler.executor import SchedulerExecutor

                start_delay = kwargs.pop(
                    "start_delay", 0.1
                )  # Reduzir delay para paralelismo
                return SchedulerExecutor(start_delay=start_delay, timeout=timeout)
            elif cpu_count >= 2 and max_workers != 1:
                logger.info("Detectado %s CPUs, usando executor paralelo", cpu_count)
                return ParallelExecutor(
                    max_workers=max_workers, timeout=timeout, **kwargs
                )
            else:
                logger.info("Usando executor com 1 worker (comportamento sequencial)")
                return ParallelExecutor(max_workers=1, timeout=timeout, **kwargs)

        else:
            raise ValueError(f"Tipo de executor inválido: {executor_type}")

    @staticmethod
    def create_sequential_executor(timeout: int = 300, **kwargs) -> ParallelExecutor:
        """
        Cria um executor sequencial (ParallelExecutor com 1 worker).

        Args:
            timeout: Timeout por tarefa em segundos
            **kwargs: Argumentos adicionais

        Returns:
            ParallelExecutor: Executor com 1 worker (comportamento sequencial)
        """
        return ParallelExecutor(max_workers=1, timeout=timeout, **kwargs)

    @staticmethod
    def create_parallel_executor(
        max_workers: Optional[int] = None, timeout: int = 300, **kwargs
    ) -> ParallelExecutor:
        """
        Cria um executor paralelo.

        Args:
            max_workers: Número máximo de workers
            timeout: Timeout por tarefa em segundos
            **kwargs: Argumentos adicionais

        Returns:
            ParallelExecutor: Executor paralelo
        """
        return ParallelExecutor(max_workers=max_workers, timeout=timeout, **kwargs)

    @staticmethod
    def recommend_executor_type(
        num_algorithms: int,
        system_load: float = 0.0,
        available_memory_mb: float = 1024.0,
    ) -> str:
        """
        Recomenda o tipo de executor baseado nas condições do sistema.

        Args:
            num_algorithms: Número de algoritmos a executar
            system_load: Carga atual do sistema (0.0 a 1.0)
            available_memory_mb: Memória disponível em MB

        Returns:
            str: Tipo de executor recomendado ("sequential" ou "parallel")
        """
        cpu_count = multiprocessing.cpu_count()

        # Verificar se vale a pena usar paralelismo
        if num_algorithms == 1:
            return "sequential"

        # Verificar recursos do sistema
        if system_load > 0.8:
            logger.info("Alta carga do sistema, recomendando executor com 1 worker")
            return "sequential"

        if available_memory_mb < 512:
            logger.info("Pouca memória disponível, recomendando executor com 1 worker")
            return "sequential"

        # Verificar se há múltiplos CPUs
        if cpu_count >= 2 and num_algorithms > 1:
            logger.info(
                f"Condições favoráveis para paralelismo: {cpu_count} CPUs, {num_algorithms} algoritmos"
            )
            return "parallel"

        return "sequential"

    @staticmethod
    def calculate_optimal_workers(
        num_algorithms: int,
        available_memory_mb: float = 1024.0,
        max_workers: Optional[int] = None,
    ) -> int:
        """
        Calcula o número ótimo de workers baseado nos recursos.

        Args:
            num_algorithms: Número de algoritmos a executar
            available_memory_mb: Memória disponível em MB
            max_workers: Limite máximo de workers

        Returns:
            int: Número ótimo de workers
        """
        cpu_count = multiprocessing.cpu_count()

        # Limite baseado no número de algoritmos
        workers_by_tasks = min(num_algorithms, cpu_count)

        # Limite baseado na memória (assumindo ~256MB por worker)
        workers_by_memory = max(1, int(available_memory_mb / 256))

        # Usar o menor dos limites
        optimal_workers = min(workers_by_tasks, workers_by_memory)

        # Aplicar limite máximo se fornecido
        if max_workers is not None:
            optimal_workers = min(optimal_workers, max_workers)

        # Garantir pelo menos 1 worker
        optimal_workers = max(1, optimal_workers)

        logger.info(
            f"Número ótimo de workers: {optimal_workers} "
            f"(CPUs: {cpu_count}, Algoritmos: {num_algorithms}, "
            f"Memória: {available_memory_mb}MB)"
        )

        return optimal_workers


# Funções de conveniência
def create_executor(
    executor_type: str = "auto",
    max_workers: Optional[int] = None,
    timeout: int = 300,
    **kwargs,
) -> IExecutor:
    """
    Função de conveniência para criar executores.

    Args:
        executor_type: Tipo de executor ("sequential", "parallel", "auto")
        max_workers: Número máximo de workers
        timeout: Timeout por tarefa em segundos
        **kwargs: Argumentos adicionais

    Returns:
        IExecutor: Instância do executor
    """
    return ExecutorFactory.create_executor(
        executor_type=executor_type, max_workers=max_workers, timeout=timeout, **kwargs
    )


def create_auto_executor(
    num_algorithms: int, timeout: int = 300, **kwargs
) -> IExecutor:
    """
    Cria um executor automaticamente baseado no número de algoritmos.

    Args:
        num_algorithms: Número de algoritmos a executar
        timeout: Timeout por tarefa em segundos
        **kwargs: Argumentos adicionais

    Returns:
        IExecutor: Executor apropriado
    """
    executor_type = ExecutorFactory.recommend_executor_type(num_algorithms)

    if executor_type == "parallel":
        max_workers = ExecutorFactory.calculate_optimal_workers(num_algorithms)
        return ExecutorFactory.create_parallel_executor(
            max_workers=max_workers, timeout=timeout, **kwargs
        )
    else:
        return ExecutorFactory.create_sequential_executor(timeout=timeout, **kwargs)
