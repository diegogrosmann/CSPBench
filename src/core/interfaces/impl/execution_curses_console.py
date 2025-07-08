"""
Implementação CursesConsole para a nova interface IConsole.

Console baseado em ncurses que implementa a interface IConsole
para execução de tarefas com interface visual avançada.
"""

import curses
import logging
from typing import Any, Dict, Optional

from ..execution_console import IConsole
from ..task_result import TaskResult

logger = logging.getLogger(__name__)


class CursesConsole:
    """
    Console baseado em ncurses para interface visual.

    Implementa IConsole usando ncurses para criar um painel
    visual com progresso em tempo real das tarefas.
    """

    def __init__(self, stdscr: Optional[Any] = None):
        """
        Inicializa o console curses.

        Args:
            stdscr: Tela curses (opcional, será criada se não fornecida)
        """
        self.stdscr = stdscr
        self.tasks: Dict[str, Dict[str, Any]] = {}
        self.task_order = []  # Para manter ordem de exibição
        self.y_offset = 2  # Offset para início das tarefas
        self.initialized = False

        if self.stdscr:
            self._init_colors()
            self.initialized = True

        logger.info("CursesConsole inicializado")

    def _init_colors(self) -> None:
        """Inicializa as cores do curses."""
        try:
            curses.start_color()
            curses.use_default_colors()

            # Definir pares de cores
            curses.init_pair(1, curses.COLOR_GREEN, -1)  # Verde
            curses.init_pair(2, curses.COLOR_YELLOW, -1)  # Amarelo
            curses.init_pair(3, curses.COLOR_RED, -1)  # Vermelho
            curses.init_pair(4, curses.COLOR_BLUE, -1)  # Azul
            curses.init_pair(5, curses.COLOR_CYAN, -1)  # Ciano

        except curses.error:
            logger.warning("Não foi possível inicializar cores curses")

    def on_task_start(self, task_id: str, meta: Dict[str, Any]) -> None:
        """
        Mostra o início de uma tarefa.

        Args:
            task_id: ID único da tarefa
            meta: Metadados da tarefa
        """
        self.tasks[task_id] = {
            "meta": meta,
            "progress": 0.0,
            "status": "running",
            "message": "Iniciando...",
        }

        if task_id not in self.task_order:
            self.task_order.append(task_id)

        self._update_display()

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
        if msg:
            self.tasks[task_id]["message"] = msg

        self._update_display()

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
        task_info["status"] = "completed" if result.success else "failed"
        task_info["progress"] = 100.0

        if result.success:
            task_info["message"] = (
                f"Concluído - Dist: {result.distance:.2f}, Tempo: {result.time:.2f}s"
            )
        else:
            task_info["message"] = f"Falhou - {result.error or 'Erro desconhecido'}"

        self._update_display()

    def _update_display(self) -> None:
        """Atualiza a tela curses."""
        if not self.initialized or not self.stdscr:
            return

        try:
            self.stdscr.clear()

            # Título
            self.stdscr.addstr(
                0,
                0,
                "🚀 Execução de Algoritmos CSP",
                curses.color_pair(4) | curses.A_BOLD,
            )
            self.stdscr.addstr(1, 0, "=" * 60)

            # Tarefas
            for i, task_id in enumerate(self.task_order):
                if task_id not in self.tasks:
                    continue

                task_info = self.tasks[task_id]
                y_pos = self.y_offset + i

                # Nome do algoritmo
                algorithm_name = task_info["meta"].get("algorithm_name", "Unknown")
                status = task_info["status"]
                progress = task_info["progress"]
                message = task_info["message"]

                # Cor baseada no status
                if status == "running":
                    color = curses.color_pair(2)  # Amarelo
                    status_icon = "🔄"
                elif status == "completed":
                    color = curses.color_pair(1)  # Verde
                    status_icon = "✅"
                else:  # failed
                    color = curses.color_pair(3)  # Vermelho
                    status_icon = "❌"

                # Linha da tarefa
                task_line = (
                    f"{status_icon} {algorithm_name:<15} [{progress:6.1f}%] {message}"
                )
                self.stdscr.addstr(y_pos, 0, task_line[:70], color)

                # Barra de progresso
                if status == "running":
                    progress_bar = self._create_progress_bar(progress)
                    self.stdscr.addstr(y_pos, 75, progress_bar, color)

            # Instruções
            instructions_y = self.y_offset + len(self.task_order) + 2
            self.stdscr.addstr(
                instructions_y,
                0,
                "Pressione Ctrl+C para cancelar",
                curses.color_pair(5),
            )

            self.stdscr.refresh()

        except curses.error as e:
            logger.error("Erro ao atualizar display curses: %s", e)

    def _create_progress_bar(self, progress: float, width: int = 20) -> str:
        """
        Cria uma barra de progresso visual.

        Args:
            progress: Progresso (0.0 a 100.0)
            width: Largura da barra

        Returns:
            str: Barra de progresso
        """
        filled = int((progress / 100.0) * width)
        empty = width - filled
        return f"[{'█' * filled}{'░' * empty}]"

    def cleanup(self) -> None:
        """
        Limpa recursos do console.
        """
        if self.initialized and self.stdscr:
            try:
                self.stdscr.clear()
                self.stdscr.refresh()
            except curses.error:
                pass

        self.tasks.clear()
        self.task_order.clear()
        logger.info("CursesConsole finalizado")
