"""Monitor simples para terminal."""

from datetime import datetime
from typing import Any, Dict, Optional

from .interfaces import (
    ExecutionLevel,
    HierarchicalContext,
    MonitoringInterface,
    TaskType,
)


class SimpleMonitor(MonitoringInterface):
    """Monitor simples que exibe progresso no terminal."""

    def __init__(self):
        self.task_type: Optional[TaskType] = None
        self.task_name: str = ""
        self.start_time: Optional[datetime] = None
        self.header_printed: bool = False

        # Controle hierárquico
        self.current_config_id: Optional[str] = None
        self.current_config_index: int = 0
        self.total_configs: int = 0

        self.current_dataset_id: Optional[str] = None
        self.current_dataset_index: int = 0
        self.total_datasets: int = 0

        # Estado dos algoritmos do dataset atual
        self.current_algorithms: Dict[str, Dict[str, Any]] = {}
        self.algorithm_order: list[str] = []

    def start_task(
        self, task_type: TaskType, task_name: str, config: Dict[str, Any]
    ) -> None:
        """Inicia monitoramento de uma tarefa."""
        self.task_type = task_type
        self.task_name = task_name
        self.start_time = datetime.now()
        self.header_printed = False
        self._print_header()

    def _print_header(self) -> None:
        """Imprime cabeçalho do monitor."""
        if self.task_type == TaskType.EXECUTION:
            print("CSPBench - Monitoramento de Execução")
        elif self.task_type == TaskType.OPTIMIZATION:
            print("CSPBench - Monitoramento de Otimização")
        elif self.task_type == TaskType.SENSITIVITY:
            print("CSPBench - Análise de Sensibilidade")
        else:
            print("CSPBench - Monitoramento")

        print("=" * 50)
        print(f"📋\tBatch: {self.task_name}")
        if self.start_time:
            print(f"⏰\tIniciado: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print()
        self.header_printed = True

    def update_hierarchy(
        self,
        level: ExecutionLevel,
        level_id: str,
        progress: float,
        message: str = "",
        data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Atualiza progresso hierárquico."""
        if level == ExecutionLevel.EXECUTION:
            # Nova configuração
            if level_id != self.current_config_id:
                self.current_config_id = level_id
                if data:
                    self.current_config_index = data.get("config_index", 1)
                    self.total_configs = data.get("total_configs", 1)

                # Reset dataset tracking quando nova configuração
                self.current_dataset_id = None

        elif level == ExecutionLevel.DATASET:
            # Novo dataset
            if level_id != self.current_dataset_id:
                self.current_dataset_id = level_id
                if data:
                    self.current_dataset_index = data.get("dataset_index", 1)
                    self.total_datasets = data.get("total_datasets", 1)
                    total_algorithms = data.get("total_algorithms", 0)

                # Reset algoritmos quando novo dataset
                self.current_algorithms = {}
                self.algorithm_order = []

                # Exibir informações do dataset
                print(
                    f"\tConfiguração: ({self.current_config_index}/{self.total_configs})"
                )
                print(
                    f"🗂️\tDataset: {level_id} ({self.current_dataset_index}/{self.total_datasets})"
                )
                print(f"🧠\tAlgoritmos: {total_algorithms} total")
                print()
                print("Progresso por algoritmo:")

    def update_item(
        self,
        item_id: str,
        progress: float,
        message: str = "",
        context: Optional[HierarchicalContext] = None,
    ) -> None:
        """Atualiza um item individual (algoritmo)."""
        algorithm_name = item_id  # item_id é o nome do algoritmo

        if context and context.algorithm_id:
            algorithm_name = context.algorithm_id

        # Primeiro algoritmo - adiciona à lista de ordem
        if algorithm_name not in self.current_algorithms:
            self.algorithm_order.append(algorithm_name)
            self.current_algorithms[algorithm_name] = {
                "progress": 0.0,
                "completed": False,
                "current_rep": 0,
                "total_rep": 1,
                "name_printed": False,
            }

        # Atualizar status
        algo_status = self.current_algorithms[algorithm_name]
        old_progress = algo_status["progress"]
        algo_status["progress"] = progress
        algo_status["completed"] = progress >= 100.0

        # Extrair informações de repetições se disponível no context
        if context and context.repetition_id:
            # Tentar extrair números de repetição do ID se for formatado como "1/10"
            rep_parts = context.repetition_id.split("/")
            if len(rep_parts) == 2:
                try:
                    algo_status["current_rep"] = int(rep_parts[0])
                    algo_status["total_rep"] = int(rep_parts[1])
                except ValueError:
                    pass

        # Decidir se deve imprimir/atualizar
        significant_change = (
            not algo_status["name_printed"]
            or abs(progress - old_progress) >= 10.0  # Mudança de 10% ou mais
            or algo_status["completed"]
        )

        if significant_change:
            if not algo_status["name_printed"]:
                print(f"• {algorithm_name}:")
                algo_status["name_printed"] = True

            # Criar linha de progresso
            progress_bar = self._create_progress_bar(progress)

            if algo_status["completed"]:
                if algo_status["total_rep"] == 1:
                    status_line = f"  {progress_bar} ✅ 100%"
                else:
                    status_line = f"  {progress_bar} ✅ 100%"
            else:
                if algo_status["total_rep"] > 1:
                    rep_info = (
                        f"({algo_status['current_rep']}/{algo_status['total_rep']}) "
                    )
                else:
                    rep_info = ""
                status_line = f"  {progress_bar} {rep_info}{progress:.0f}%"

            print(status_line)

    def algorithm_callback(
        self,
        algorithm_name: str,
        progress: float,
        message: str,
        item_id: Optional[str] = None,
    ) -> None:
        """Callback direto do algoritmo durante execução."""
        # Atualizar progresso interno através do algoritmo callback
        context = HierarchicalContext(algorithm_id=algorithm_name)
        self.update_item(item_id or algorithm_name, progress, message, context)

    def finish_item(
        self,
        item_id: str,
        success: bool = True,
        result: Optional[Dict[str, Any]] = None,
        error: Optional[str] = None,
    ) -> None:
        """Finaliza um item individual."""
        # Marcar como 100% se bem sucedido
        if success:
            context = HierarchicalContext(algorithm_id=item_id)
            self.update_item(item_id, 100.0, "Concluído", context)

    def finish_task(
        self,
        success: bool = True,
        final_results: Optional[Dict[str, Any]] = None,
        error_message: str = "",
    ) -> None:
        """Finaliza monitoramento da tarefa."""
        print()
        if success:
            print("✅ Monitoramento concluído!")
        else:
            print(f"❌ Erro: {error_message}")

        if self.start_time:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            print(f"⏰ Tempo total: {elapsed/60:.0f}m {elapsed%60:.0f}s")
        print("=" * 50)

    def _create_progress_bar(self, progress: float) -> str:
        """Cria uma barra de progresso ASCII."""
        return "━" * 32  # Barra completa como no exemplo

    def get_summary(self) -> Dict[str, Any]:
        """Retorna resumo consolidado da execução atual."""
        return {
            "task_type": self.task_type.value if self.task_type else "unknown",
            "task_name": self.task_name,
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "current_config": self.current_config_id,
            "current_dataset": self.current_dataset_id,
            "algorithms": dict(self.current_algorithms),
        }

    def show_error(self, error: str) -> None:
        """Exibe erro ao usuário."""
        print(f"\n❌ Erro: {error}")

    def stop(self) -> None:
        """Para monitoramento."""
        pass

    def close(self) -> None:
        """Fecha sistema de monitoramento."""
        pass
