"""
Implementação SimpleConsole para a nova interface IConsole.

Console baseado em prints padrão que implementa a interface IConsole
para execução de tarefas com eventos padronizados.
"""

import logging
from typing import Any, Dict

from ..execution_console import IConsole
from ..task_result import TaskResult

logger = logging.getLogger(__name__)


class SimpleConsole:
    """
    Console simples baseado em prints padrão.

    Implementa IConsole usando prints padrão do Python,
    adequado para execução em terminais simples ou logs.
    """

    def __init__(self, verbose: bool = True):
        """
        Inicializa o console simples.

        Args:
            verbose: Se True, mostra mensagens detalhadas
        """
        self.verbose = verbose
        self.tasks: Dict[str, Dict[str, Any]] = {}
        logger.info("SimpleConsole inicializado")

    def on_task_start(self, task_id: str, meta: Dict[str, Any]) -> None:
        """
        Mostra o início de uma tarefa.

        Args:
            task_id: ID único da tarefa
            meta: Metadados da tarefa
        """
        self.tasks[task_id] = {"meta": meta, "start_time": None, "progress": 0.0}

        algorithm_name = meta.get("algorithm_name", "Unknown")
        if self.verbose:
            print(f"🚀 Iniciando: {algorithm_name} (ID: {task_id})")
        else:
            print(f"🚀 {algorithm_name}")

    def on_task_progress(self, task_id: str, pct: float, msg: str = "") -> None:
        """
        Atualiza o progresso de uma tarefa.

        Args:
            task_id: ID único da tarefa
            pct: Porcentagem de progresso
            msg: Mensagem opcional
        """
        if task_id not in self.tasks:
            return

        self.tasks[task_id]["progress"] = pct

        if self.verbose:
            algorithm_name = self.tasks[task_id]["meta"].get(
                "algorithm_name", "Unknown"
            )
            progress_msg = f"📊 {algorithm_name}: {pct:.1f}%"
            if msg:
                progress_msg += f" - {msg}"
            print(progress_msg)

    def on_task_finish(self, task_id: str, result: TaskResult) -> None:
        """
        Mostra o fim de uma tarefa.

        Args:
            task_id: ID único da tarefa
            result: Resultado da execução
        """
        if task_id not in self.tasks:
            return

        task_info = self.tasks[task_id]
        algorithm_name = task_info["meta"].get("algorithm_name", "Unknown")

        if result.success:
            print(
                f"✅ Concluído: {algorithm_name} - "
                f"Distância: {result.distance:.2f}, "
                f"Tempo: {result.time:.2f}s"
            )
        else:
            error_msg = result.error or "Erro desconhecido"
            print(f"❌ Falhou: {algorithm_name} - {error_msg}")

        # Remover tarefa da lista
        del self.tasks[task_id]

    def cleanup(self) -> None:
        """
        Limpa recursos do console.
        """
        if self.verbose:
            print("🧹 Limpando SimpleConsole...")
        self.tasks.clear()
        logger.info("SimpleConsole finalizado")
