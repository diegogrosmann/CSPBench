"""
Resultado padronizado de execução de tarefas.

Define o TaskResult que encapsula o resultado de uma tarefa executada,
incluindo informações de sucesso, erro, tempo de execução e dados do algoritmo.
"""

import time as time_module
from dataclasses import dataclass
from typing import Any, Dict, Optional


@dataclass(frozen=True)
class TaskResult:
    """
    Resultado padronizado de execução de tarefa.

    Encapsula todas as informações relevantes sobre a execução de uma tarefa,
    incluindo sucesso/falha, resultados do algoritmo, tempo e erros.

    Attributes:
        success: Se a tarefa foi executada com sucesso
        distance: Distância encontrada pelo algoritmo (None se falhou)
        center: String centro encontrada (None se falhou)
        time: Tempo de execução em segundos
        error: Mensagem de erro (None se sucesso)
        traceback: Stack trace completo (None se sucesso)
        metadata: Metadados adicionais da execução
    """

    success: bool
    distance: Optional[float]
    center: Optional[str]
    time: float
    error: Optional[str]
    traceback: Optional[str]
    metadata: Dict[str, Any]

    @classmethod
    def success_result(
        cls,
        distance: float,
        center: str,
        time: float,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> "TaskResult":
        """
        Cria um TaskResult para execução bem-sucedida.

        Args:
            distance: Distância encontrada
            center: String centro encontrada
            time: Tempo de execução
            metadata: Metadados adicionais

        Returns:
            TaskResult: Resultado de sucesso
        """
        return cls(
            success=True,
            distance=distance,
            center=center,
            time=time,
            error=None,
            traceback=None,
            metadata=metadata or {},
        )

    @classmethod
    def failure_result(
        cls,
        error: str,
        time: float,
        traceback: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> "TaskResult":
        """
        Cria um TaskResult para execução falhada.

        Args:
            error: Mensagem de erro
            time: Tempo de execução até falha
            traceback: Stack trace completo
            metadata: Metadados adicionais

        Returns:
            TaskResult: Resultado de falha
        """
        return cls(
            success=False,
            distance=None,
            center=None,
            time=time,
            error=error,
            traceback=traceback,
            metadata=metadata or {},
        )

    @classmethod
    def timeout_result(
        cls,
        time: float,
        timeout_limit: float,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> "TaskResult":
        """
        Cria um TaskResult para timeout.

        Args:
            time: Tempo de execução até timeout
            timeout_limit: Limite de timeout configurado
            metadata: Metadados adicionais

        Returns:
            TaskResult: Resultado de timeout
        """
        timeout_metadata = metadata or {}
        timeout_metadata.update(
            {"timeout_limit": timeout_limit, "timeout_occurred": True}
        )

        return cls(
            success=False,
            distance=None,
            center=None,
            time=time,
            error=f"Timeout após {timeout_limit}s",
            traceback=None,
            metadata=timeout_metadata,
        )

    def to_dict(self) -> Dict[str, Any]:
        """
        Converte o TaskResult para dicionário.

        Returns:
            Dict: Representação em dicionário
        """
        return {
            "success": self.success,
            "distance": self.distance,
            "center": self.center,
            "time": self.time,
            "error": self.error,
            "traceback": self.traceback,
            "metadata": self.metadata,
        }

    def to_json_log(self, task_id: str, algorithm_name: str) -> Dict[str, Any]:
        """
        Converte para formato de log JSON estruturado.

        Args:
            task_id: ID da tarefa
            algorithm_name: Nome do algoritmo

        Returns:
            Dict: Dados para log JSON
        """
        return {
            "timestamp": time_module.time(),
            "task_id": task_id,
            "algorithm": algorithm_name,
            "success": self.success,
            "distance": self.distance,
            "center": self.center,
            "execution_time": self.time,
            "error": self.error,
            "has_traceback": self.traceback is not None,
            "metadata": self.metadata,
        }

    def get_status_symbol(self) -> str:
        """
        Retorna símbolo visual do status da tarefa.

        Returns:
            str: Símbolo (✓ para sucesso, ❌ para erro, ⏰ para timeout)
        """
        if self.success:
            return "✓"
        elif self.error and "timeout" in self.error.lower():
            return "⏰"
        else:
            return "❌"

    def get_summary(self) -> str:
        """
        Retorna resumo textual do resultado.

        Returns:
            str: Resumo para exibição
        """
        symbol = self.get_status_symbol()

        if self.success:
            return f"{symbol} Sucesso: distância={self.distance:.2f}, tempo={self.time:.2f}s"
        else:
            return f"{symbol} Falha: {self.error}, tempo={self.time:.2f}s"

    def __str__(self) -> str:
        """Representação string do TaskResult."""
        return self.get_summary()
