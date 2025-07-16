"""Interfaces e tipos para sistema de monitoramento."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from threading import RLock
from typing import Any, Dict, List, Optional


class TaskType(Enum):
    """Tipos de tarefa suportadas."""

    EXECUTION = "execution"
    OPTIMIZATION = "optimization"
    SENSITIVITY = "sensitivity"


class TaskStatus(Enum):
    """Status de uma tarefa."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ItemStatus(Enum):
    """Status de um item individual (repetição, trial, amostra)."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ExecutionLevel(Enum):
    """Níveis hierárquicos de execução."""

    EXECUTION = "execution"
    DATASET = "dataset"
    ALGORITHM = "algorithm"
    REPETITION = "repetition"
    TRIAL = "trial"
    SAMPLE = "sample"


@dataclass
class HierarchicalContext:
    """Contexto hierárquico para localizar item na execução."""

    execution_id: Optional[str] = None
    dataset_id: Optional[str] = None
    algorithm_id: Optional[str] = None
    repetition_id: Optional[str] = None
    trial_id: Optional[str] = None
    sample_id: Optional[str] = None

    def get_path(self) -> str:
        """Retorna caminho hierárquico como string."""
        parts = []
        if self.execution_id:
            parts.append(f"exec:{self.execution_id}")
        if self.dataset_id:
            parts.append(f"dataset:{self.dataset_id}")
        if self.algorithm_id:
            parts.append(f"algo:{self.algorithm_id}")
        if self.repetition_id:
            parts.append(f"rep:{self.repetition_id}")
        if self.trial_id:
            parts.append(f"trial:{self.trial_id}")
        if self.sample_id:
            parts.append(f"sample:{self.sample_id}")
        return "/".join(parts)

    def get_level_id(self, level: ExecutionLevel) -> Optional[str]:
        """Retorna ID do nível específico."""
        level_map = {
            ExecutionLevel.EXECUTION: self.execution_id,
            ExecutionLevel.DATASET: self.dataset_id,
            ExecutionLevel.ALGORITHM: self.algorithm_id,
            ExecutionLevel.REPETITION: self.repetition_id,
            ExecutionLevel.TRIAL: self.trial_id,
            ExecutionLevel.SAMPLE: self.sample_id,
        }
        return level_map.get(level)


@dataclass
class TaskSpecificData:
    """Dados específicos por tipo de tarefa - thread-safe e simplificado."""

    # Progresso hierárquico simplificado
    total_executions: int = 0
    current_execution_index: int = 0
    current_execution_name: str = ""

    total_datasets: int = 0
    current_dataset_index: int = 0
    current_dataset_name: str = ""

    total_algorithms: int = 0
    current_algorithm_index: int = 0
    current_algorithm_name: str = ""

    # Estado consolidado por algoritmo
    algorithm_data: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # Resultados consolidados
    best_distance: Optional[float] = None
    best_algorithm: str = ""
    current_best_result: Optional[Dict[str, Any]] = None

    # Dados específicos por tipo de tarefa
    task_specific: Dict[str, Any] = field(default_factory=dict)

    # Lock para thread-safety
    _lock: RLock = field(default_factory=RLock, init=False, repr=False)

    def get_execution_progress(self) -> float:
        """Calcula progresso geral da execução - thread-safe."""
        with self._lock:
            if self.total_executions == 0:
                return 0.0
            return self.current_execution_index / self.total_executions

    def get_dataset_progress(self) -> float:
        """Calcula progresso dos datasets - thread-safe."""
        with self._lock:
            if self.total_datasets == 0:
                return 0.0
            return self.current_dataset_index / self.total_datasets

    def get_algorithm_progress(self) -> float:
        """Calcula progresso dos algoritmos - thread-safe."""
        with self._lock:
            if self.total_algorithms == 0:
                return 0.0
            return self.current_algorithm_index / self.total_algorithms

    def update_algorithm_data(self, algorithm_name: str, **data) -> None:
        """Atualiza dados de algoritmo de forma thread-safe."""
        with self._lock:
            if algorithm_name not in self.algorithm_data:
                self.algorithm_data[algorithm_name] = {}
            self.algorithm_data[algorithm_name].update(data)

    def get_algorithm_data(self, algorithm_name: str) -> Dict[str, Any]:
        """Obtém dados de algoritmo de forma thread-safe."""
        with self._lock:
            return self.algorithm_data.get(algorithm_name, {}).copy()

    def update_best_result(
        self, algorithm: str, distance: float, result: Dict[str, Any]
    ) -> None:
        """Atualiza melhor resultado se necessário - thread-safe."""
        with self._lock:
            if self.best_distance is None or distance > self.best_distance:
                self.best_distance = distance
                self.best_algorithm = algorithm
                self.current_best_result = result.copy()

    def set_task_specific(self, key: str, value: Any) -> None:
        """Define valor específico da tarefa - thread-safe."""
        with self._lock:
            self.task_specific[key] = value

    def get_task_specific(self, key: str, default: Any = None) -> Any:
        """Obtém valor específico da tarefa - thread-safe."""
        with self._lock:
            return self.task_specific.get(key, default)


@dataclass
class TaskItem:
    """Representa um item individual de trabalho (repetição, trial, amostra)."""

    item_id: str
    item_type: str  # "repetition", "trial", "sample", etc.
    status: ItemStatus
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    progress: float = 0.0
    message: str = ""
    result: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    worker_id: Optional[str] = None
    context: Optional[HierarchicalContext] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def elapsed_time(self) -> Optional[timedelta]:
        """Retorna tempo decorrido da execução."""
        if self.start_time is None:
            return None
        end = self.end_time or datetime.now()
        return end - self.start_time

    @property
    def is_active(self) -> bool:
        """Verifica se item está ativo (rodando)."""
        return self.status == ItemStatus.RUNNING

    @property
    def duration(self) -> Optional[timedelta]:
        """Retorna duração total da execução (para itens finalizados)."""
        if self.start_time is None or self.end_time is None:
            return None
        return self.end_time - self.start_time

    @property
    def hierarchical_path(self) -> str:
        """Retorna caminho hierárquico completo."""
        if self.context:
            return self.context.get_path()
        return self.item_id

    @property
    def is_finished(self) -> bool:
        """Verifica se o item terminou (sucesso ou falha)."""
        return self.status in (
            ItemStatus.COMPLETED,
            ItemStatus.FAILED,
            ItemStatus.CANCELLED,
        )


@dataclass
class TaskProgress:
    """Estrutura thread-safe para progresso de qualquer tipo de tarefa."""

    # Identificação da tarefa
    task_type: TaskType
    task_name: str
    task_id: str = ""

    # Progresso geral
    total_items: int = 0
    completed_items: int = 0
    failed_items: int = 0
    progress_percentage: float = 0.0

    # Item atual em processamento
    current_item: str = ""
    current_item_progress: float = 0.0
    current_item_message: str = ""

    # Timing
    start_time: Optional[datetime] = None
    current_time: datetime = field(default_factory=datetime.now)
    estimated_completion: Optional[datetime] = None

    # Status e controle
    status: TaskStatus = TaskStatus.PENDING
    is_active: bool = True

    # Dados específicos por tipo de tarefa (uso delegado)
    task_data: TaskSpecificData = field(default_factory=TaskSpecificData)

    # Itens individuais sendo processados - thread-safe
    items: Dict[str, TaskItem] = field(default_factory=dict)
    active_items: List[str] = field(default_factory=list)

    # Mensagens e feedback
    messages: List[str] = field(default_factory=list)
    callback_info: str = ""

    # Lock para thread-safety
    _lock: RLock = field(default_factory=RLock, init=False, repr=False)

    @property
    def elapsed_time(self) -> Optional[timedelta]:
        """Retorna tempo decorrido desde o início - thread-safe."""
        with self._lock:
            if self.start_time is None:
                return None
            return datetime.now() - self.start_time

    @property
    def running_items_count(self) -> int:
        """Retorna número de itens em execução - thread-safe."""
        with self._lock:
            return len(self.active_items)

    def update_progress(self) -> None:
        """Atualiza progresso baseado nos itens - thread-safe."""
        with self._lock:
            if self.total_items > 0:
                self.progress_percentage = (
                    self.completed_items / self.total_items * 100.0
                )

            # Atualiza contadores
            self.completed_items = sum(
                1 for item in self.items.values() if item.status == ItemStatus.COMPLETED
            )
            self.failed_items = sum(
                1 for item in self.items.values() if item.status == ItemStatus.FAILED
            )

            # Atualiza lista de itens ativos
            self.active_items = [
                item_id
                for item_id, item in self.items.items()
                if item.status == ItemStatus.RUNNING
            ]

    def add_message(self, message: str) -> None:
        """Adiciona mensagem com timestamp - thread-safe."""
        with self._lock:
            timestamp = datetime.now().strftime("%H:%M:%S")
            self.messages.append(f"[{timestamp}] {message}")

    def add_item(self, item: TaskItem) -> None:
        """Adiciona item de forma thread-safe."""
        with self._lock:
            self.items[item.item_id] = item
            if item.status == ItemStatus.RUNNING:
                self.active_items.append(item.item_id)

    def update_item(self, item_id: str, **updates) -> None:
        """Atualiza item de forma thread-safe."""
        with self._lock:
            if item_id in self.items:
                item = self.items[item_id]
                for key, value in updates.items():
                    if hasattr(item, key):
                        setattr(item, key, value)

                # Atualizar lista de itens ativos
                if (
                    item.status == ItemStatus.RUNNING
                    and item_id not in self.active_items
                ):
                    self.active_items.append(item_id)
                elif item.status != ItemStatus.RUNNING and item_id in self.active_items:
                    self.active_items.remove(item_id)

    def get_summary(self) -> Dict[str, Any]:
        """Retorna resumo da execução - thread-safe."""
        with self._lock:
            return {
                "task_type": self.task_type.value,
                "task_name": self.task_name,
                "progress": self.progress_percentage,
                "completed": self.completed_items,
                "failed": self.failed_items,
                "total": self.total_items,
                "running": self.running_items_count,
                "status": self.status.value,
                "elapsed_time": str(self.elapsed_time) if self.elapsed_time else None,
                "current_item": self.current_item,
                "task_specific": {
                    "best_distance": self.task_data.best_distance,
                    "best_algorithm": self.task_data.best_algorithm,
                    "current_dataset": self.task_data.current_dataset_name,
                    "current_algorithm": self.task_data.current_algorithm_name,
                },
            }

    def get_algorithm_summary(self, algorithm_name: str) -> Dict[str, Any]:
        """Retorna resumo de um algoritmo específico."""
        algorithm_items = [
            item
            for item in self.items.values()
            if item.context and item.context.algorithm_id == algorithm_name
        ]

        completed = sum(
            1 for item in algorithm_items if item.status == ItemStatus.COMPLETED
        )
        failed = sum(1 for item in algorithm_items if item.status == ItemStatus.FAILED)
        running = sum(
            1 for item in algorithm_items if item.status == ItemStatus.RUNNING
        )

        best_result = None
        best_distance = None

        for item in algorithm_items:
            if item.result and "max_distance" in item.result:
                distance = item.result["max_distance"]
                if best_distance is None or distance > best_distance:
                    best_distance = distance
                    best_result = item.result

        # Calcular progresso baseado nos itens do algoritmo
        total_items = len(algorithm_items)
        progress = (completed / total_items * 100.0) if total_items > 0 else 0.0

        # Determinar status baseado no estado dos itens
        if running > 0:
            status = ItemStatus.RUNNING
        elif failed > 0 and completed == 0:
            status = ItemStatus.FAILED
        elif completed == total_items and total_items > 0:
            status = ItemStatus.COMPLETED
        else:
            status = ItemStatus.PENDING

        # Obter dados adicionais do algoritmo
        algorithm_data = self.task_data.get_algorithm_data(algorithm_name)

        return {
            "algorithm": algorithm_name,
            "total_items": total_items,
            "completed": completed,
            "failed": failed,
            "running": running,
            "progress": progress,
            "status": status.value,
            "best_distance": best_distance,
            "best_result": best_result,
            "algorithm_data": algorithm_data,
        }


class MonitoringInterface(ABC):
    """Interface simplificada e thread-safe para sistemas de monitoramento."""

    @abstractmethod
    def start_task(
        self, task_type: TaskType, task_name: str, config: Dict[str, Any]
    ) -> None:
        """
        Inicia monitoramento de uma tarefa.

        Args:
            task_type: Tipo da tarefa (execution, optimization, sensitivity)
            task_name: Nome da tarefa
            config: Configuração da tarefa
        """
        ...

    @abstractmethod
    def start_item(
        self,
        item_id: str,
        item_type: str = "repetition",
        context: Optional[HierarchicalContext] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Inicia um item individual (repetição, trial, amostra).

        Args:
            item_id: ID único do item
            item_type: Tipo do item (repetition, trial, sample, etc.)
            context: Contexto hierárquico
            metadata: Metadados opcionais
        """
        ...

    @abstractmethod
    def update_item(
        self,
        item_id: str,
        progress: float,
        message: str = "",
        context: Optional[HierarchicalContext] = None,
    ) -> None:
        """
        Atualiza um item individual (repetição, trial, amostra).

        Args:
            item_id: ID único do item
            progress: Progresso (0.0 a 1.0)
            message: Mensagem de status
            context: Contexto hierárquico
        """
        ...

    @abstractmethod
    def finish_item(
        self,
        item_id: str,
        success: bool = True,
        result: Optional[Dict[str, Any]] = None,
        error: Optional[str] = None,
    ) -> None:
        """
        Finaliza um item individual.

        Args:
            item_id: ID único do item
            success: Se foi executado com sucesso
            result: Resultado da execução
            error: Mensagem de erro se falhou
        """
        ...

    @abstractmethod
    def update_hierarchy(
        self,
        level: ExecutionLevel,
        level_id: str,
        progress: float,
        message: str = "",
        data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Atualiza progresso hierárquico (execução, dataset, algoritmo).

        Args:
            level: Nível hierárquico
            level_id: ID do nível
            progress: Progresso (0.0 a 1.0)
            message: Mensagem de status
            data: Dados adicionais específicos do nível
        """
        ...

    @abstractmethod
    def algorithm_callback(
        self,
        algorithm_name: str,
        progress: float,
        message: str,
        item_id: Optional[str] = None,
    ) -> None:
        """
        Callback direto do algoritmo durante execução.

        Args:
            algorithm_name: Nome do algoritmo
            progress: Progresso do algoritmo (0.0 a 1.0)
            message: Mensagem do algoritmo
            item_id: ID do item sendo executado (se aplicável)
        """
        ...

    @abstractmethod
    def finish_task(
        self,
        success: bool = True,
        final_results: Optional[Dict[str, Any]] = None,
        error_message: str = "",
    ) -> None:
        """
        Finaliza monitoramento da tarefa.

        Args:
            success: Se a tarefa foi concluída com sucesso
            final_results: Resultados finais consolidados
            error_message: Mensagem de erro se falhou
        """
        ...

    @abstractmethod
    def get_summary(self) -> Dict[str, Any]:
        """Retorna resumo consolidado da execução atual."""
        ...

    @abstractmethod
    def show_error(self, error: str) -> None:
        """Exibe erro ao usuário."""
        ...

    @abstractmethod
    def stop(self) -> None:
        """Para monitoramento."""
        ...

    @abstractmethod
    def close(self) -> None:
        """Fecha sistema de monitoramento."""
        ...


@dataclass
class TaskConfiguration:
    """Configuração simplificada para inicialização de tarefas."""

    task_type: TaskType
    task_name: str
    total_items: int
    task_specific_config: Dict[str, Any] = field(default_factory=dict)
    timeout: Optional[int] = None
    parallel_execution: bool = False
    max_workers: Optional[int] = None


def create_task_progress(
    task_type: TaskType, task_name: str, total_items: int = 0, **task_data
) -> TaskProgress:
    """
    Factory function para criar TaskProgress com validação.

    Args:
        task_type: Tipo da tarefa
        task_name: Nome da tarefa
        total_items: Total de itens a processar
        **task_data: Dados específicos da tarefa

    Returns:
        TaskProgress: Objeto de progresso inicializado

    Raises:
        ValueError: Se parâmetros obrigatórios estão ausentes
    """
    if not task_name:
        raise ValueError("task_name é obrigatório")

    if total_items < 0:
        raise ValueError("total_items deve ser >= 0")

    # Criar TaskSpecificData com dados fornecidos
    specific_data = TaskSpecificData(**task_data)

    return TaskProgress(
        task_type=task_type,
        task_name=task_name,
        total_items=total_items,
        task_data=specific_data,
        start_time=datetime.now(),
        status=TaskStatus.PENDING,
    )


def create_task_item(
    item_id: str,
    item_type: str,
    context: Optional[HierarchicalContext] = None,
    metadata: Optional[Dict[str, Any]] = None,
) -> TaskItem:
    """
    Factory function para criar TaskItem com validação.

    Args:
        item_id: ID único do item
        item_type: Tipo do item (repetition, trial, sample, etc.)
        context: Contexto hierárquico do item
        metadata: Metadados opcionais

    Returns:
        TaskItem: Objeto de item inicializado

    Raises:
        ValueError: Se parâmetros obrigatórios estão ausentes
    """
    if not item_id:
        raise ValueError("item_id é obrigatório")

    if not item_type:
        raise ValueError("item_type é obrigatório")

    return TaskItem(
        item_id=item_id,
        item_type=item_type,
        status=ItemStatus.PENDING,
        context=context,
        metadata=metadata or {},
    )


def create_hierarchical_context(
    execution_id: Optional[str] = None,
    dataset_id: Optional[str] = None,
    algorithm_id: Optional[str] = None,
    repetition_id: Optional[str] = None,
    trial_id: Optional[str] = None,
    sample_id: Optional[str] = None,
) -> HierarchicalContext:
    """
    Factory function para criar contexto hierárquico.

    Args:
        execution_id: ID da execução
        dataset_id: ID do dataset
        algorithm_id: ID do algoritmo
        repetition_id: ID da repetição
        trial_id: ID do trial (otimização)
        sample_id: ID da amostra (sensibilidade)

    Returns:
        HierarchicalContext: Contexto hierárquico inicializado
    """
    return HierarchicalContext(
        execution_id=execution_id,
        dataset_id=dataset_id,
        algorithm_id=algorithm_id,
        repetition_id=repetition_id,
        trial_id=trial_id,
        sample_id=sample_id,
    )


def create_task_configuration(
    task_type: TaskType, task_name: str, total_items: int, **config
) -> TaskConfiguration:
    """
    Factory function para criar configuração de tarefa.

    Args:
        task_type: Tipo da tarefa
        task_name: Nome da tarefa
        total_items: Total de itens a processar
        **config: Configurações adicionais

    Returns:
        TaskConfiguration: Configuração inicializada

    Raises:
        ValueError: Se parâmetros obrigatórios estão ausentes
    """
    if not task_name:
        raise ValueError("task_name é obrigatório")

    if total_items < 0:
        raise ValueError("total_items deve ser >= 0")

    return TaskConfiguration(
        task_type=task_type,
        task_name=task_name,
        total_items=total_items,
        task_specific_config=config,
    )
