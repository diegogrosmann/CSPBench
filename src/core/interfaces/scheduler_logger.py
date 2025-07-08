"""
Sistema de log estruturado para tarefas do ExecutionScheduler.

Implementa logging JSON estruturado para todas as execuções de tarefas,
permitindo análise detalhada e auditoria das execuções.
"""

import json
import logging
import os
from datetime import datetime
from pathlib import Path
from threading import Lock
from typing import Any, Dict, Optional

from .task_result import TaskResult


class SchedulerLogger:
    """
    Logger estruturado para tarefas do ExecutionScheduler.

    Grava logs em formato JSON-line para facilitar análise e parsing.
    Cada linha representa uma execução de tarefa completa.
    """

    def __init__(self, log_dir: str = "outputs/logs"):
        """
        Inicializa o logger.

        Args:
            log_dir: Diretório para arquivos de log
        """
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.log_file = self.log_dir / "scheduler.log"
        self._lock = Lock()

        # Configurar logger para arquivo apenas (sem console)
        self.logger = logging.getLogger("scheduler")
        self.logger.setLevel(logging.INFO)

        # Evitar duplicação de handlers
        if not self.logger.handlers:
            # Adicionar apenas handler de arquivo, não console
            file_handler = logging.FileHandler(self.log_file)
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
            # Evitar propagação para loggers pais (que podem ter handlers de console)
            self.logger.propagate = False

    def log_task_start(
        self,
        task_id: str,
        algorithm_name: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Loga início de execução de tarefa.

        Args:
            task_id: ID da tarefa
            algorithm_name: Nome do algoritmo
            metadata: Metadados adicionais
        """
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "event": "task_start",
            "task_id": task_id,
            "algorithm": algorithm_name,
            "metadata": metadata or {},
        }

        self._write_log_entry(log_entry)

    def log_task_end(
        self,
        task_id: str,
        algorithm_name: str,
        result: TaskResult,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Loga fim de execução de tarefa.

        Args:
            task_id: ID da tarefa
            algorithm_name: Nome do algoritmo
            result: Resultado da execução
            metadata: Metadados adicionais
        """
        # Criar entrada de log estruturada
        combined_metadata = result.metadata.copy()
        combined_metadata.update(metadata or {})

        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "event": "task_end",
            "task_id": task_id,
            "algorithm": algorithm_name,
            "success": result.success,
            "distance": result.distance,
            "center": result.center,
            "execution_time": result.time,
            "error": result.error,
            "has_traceback": result.traceback is not None,
            "metadata": combined_metadata,
        }

        self._write_log_entry(log_entry)

    def log_scheduler_event(self, event_type: str, details: Dict[str, Any]) -> None:
        """
        Loga evento do scheduler.

        Args:
            event_type: Tipo do evento
            details: Detalhes do evento
        """
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "event": f"scheduler_{event_type}",
            **details,
        }

        self._write_log_entry(log_entry)

    def log_resource_snapshot(
        self,
        cpu_percent: float,
        memory_mb: float,
        load_avg: float,
        active_tasks: int,
        queue_size: int,
    ) -> None:
        """
        Loga snapshot de recursos.

        Args:
            cpu_percent: Uso de CPU em %
            memory_mb: Memória disponível em MB
            load_avg: Load average
            active_tasks: Tarefas ativas
            queue_size: Tamanho da fila
        """
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "event": "resource_snapshot",
            "cpu_percent": cpu_percent,
            "memory_mb": memory_mb,
            "load_avg": load_avg,
            "active_tasks": active_tasks,
            "queue_size": queue_size,
        }

        self._write_log_entry(log_entry)

    def _write_log_entry(self, entry: Dict[str, Any]) -> None:
        """
        Escreve entrada de log no arquivo.

        Args:
            entry: Entrada de log
        """
        with self._lock:
            try:
                with open(self.log_file, "a", encoding="utf-8") as f:
                    json.dump(entry, f, ensure_ascii=False, default=str)
                    f.write("\n")
            except Exception as e:
                # Erro ao escrever log - apenas continua silenciosamente
                pass

    def get_log_stats(self) -> Dict[str, Any]:
        """
        Retorna estatísticas dos logs.

        Returns:
            Dict: Estatísticas dos logs
        """
        if not self.log_file.exists():
            return {
                "total_entries": 0,
                "task_starts": 0,
                "task_ends": 0,
                "successes": 0,
                "failures": 0,
                "file_size": 0,
            }

        stats = {
            "total_entries": 0,
            "task_starts": 0,
            "task_ends": 0,
            "successes": 0,
            "failures": 0,
            "file_size": self.log_file.stat().st_size,
        }

        try:
            with open(self.log_file, "r", encoding="utf-8") as f:
                for line in f:
                    try:
                        entry = json.loads(line.strip())
                        stats["total_entries"] += 1

                        if entry.get("event") == "task_start":
                            stats["task_starts"] += 1
                        elif entry.get("event") == "task_end":
                            stats["task_ends"] += 1
                            if entry.get("success"):
                                stats["successes"] += 1
                            else:
                                stats["failures"] += 1

                    except json.JSONDecodeError:
                        continue

        except Exception as e:
            # Erro ao ler estatísticas - apenas continua silenciosamente
            pass

        return stats

    def clear_logs(self) -> None:
        """Remove todos os logs."""
        try:
            if self.log_file.exists():
                self.log_file.unlink()
        except Exception as e:
            # Erro ao limpar logs - apenas continua silenciosamente
            pass


# Singleton global do logger
_scheduler_logger: Optional[SchedulerLogger] = None


def get_scheduler_logger() -> SchedulerLogger:
    """
    Retorna instância singleton do SchedulerLogger.

    Returns:
        SchedulerLogger: Instância do logger
    """
    global _scheduler_logger
    if _scheduler_logger is None:
        _scheduler_logger = SchedulerLogger()
    return _scheduler_logger


def log_task_start(
    task_id: str, algorithm_name: str, metadata: Optional[Dict[str, Any]] = None
) -> None:
    """Função conveniente para log de início de tarefa."""
    get_scheduler_logger().log_task_start(task_id, algorithm_name, metadata)


def log_task_end(
    task_id: str,
    algorithm_name: str,
    result: TaskResult,
    metadata: Optional[Dict[str, Any]] = None,
) -> None:
    """Função conveniente para log de fim de tarefa."""
    get_scheduler_logger().log_task_end(task_id, algorithm_name, result, metadata)
