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
        self.last_execution_name: str = ""  # Controlar duplicação

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

        # Mapeamento de repetições para algoritmos
        self.algorithm_repetitions: Dict[str, Dict[str, Any]] = {}

    def start_task(
        self, task_type: TaskType, task_name: str, config: Dict[str, Any]
    ) -> None:
        """Inicia monitoramento de uma tarefa."""
        self.task_type = task_type
        self.task_name = task_name
        self.start_time = datetime.now()
        self.header_printed = False
        self._print_header()

    def start_item(
        self,
        item_id: str,
        item_type: str = "repetition",
        context: Optional[HierarchicalContext] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Inicia um item individual (repetição, trial, amostra)."""
        # Usar sempre o algorithm_id do contexto como chave principal
        algorithm_name = (
            context.algorithm_id if context and context.algorithm_id else item_id
        )

        # Primeiro algoritmo - adiciona à lista de ordem
        if algorithm_name not in self.current_algorithms:
            self.algorithm_order.append(algorithm_name)
            self.current_algorithms[algorithm_name] = {
                "progress": 0.0,
                "completed": False,
                "current_rep": 0,
                "total_rep": 1,
                "name_printed": False,
                "completed_reps": 0,
            }
            self.algorithm_repetitions[algorithm_name] = {}

        # Inicializar estrutura para repetição específica
        if item_id not in self.algorithm_repetitions[algorithm_name]:
            # Extrair informações de repetições se disponível no context
            current_rep = 1
            total_rep = 1
            if context and context.repetition_id:
                # Tentar extrair números de repetição do ID se for formatado como "1/10"
                rep_parts = context.repetition_id.split("/")
                if len(rep_parts) == 2:
                    try:
                        current_rep = int(rep_parts[0])
                        total_rep = int(rep_parts[1])
                    except ValueError:
                        pass  # Manter valores padrão se não conseguir parsear

            self.algorithm_repetitions[algorithm_name][item_id] = {
                "progress": 0.0,
                "completed": False,
                "rep_number": current_rep,
                "started": True,
                "start_time": datetime.now(),
            }

            # Atualizar total de repetições do algoritmo
            self.current_algorithms[algorithm_name]["total_rep"] = total_rep

            # Imprimir nome do algoritmo na primeira repetição
            if not self.current_algorithms[algorithm_name]["name_printed"]:
                print("")

                # Extrair informações do contexto ou item_id para mostrar detalhes
                dataset_info = ""
                if context and context.dataset_id:
                    dataset_info = f" | Dataset: {context.dataset_id}"
                elif "_" in item_id:
                    # Extrair dataset do item_id se formatado como algorithm_dataset
                    parts = item_id.split("_", 1)
                    if len(parts) > 1:
                        dataset_info = f" | Dataset: {parts[1]}"

                # Mostrar barra de progresso inicial
                progress_bar = self._create_progress_bar(0.0)
                rep_info = f"(0/{total_rep}) " if total_rep > 1 else ""
                status_line = f"• {algorithm_name} {progress_bar} {rep_info}0%"
                print(f"{status_line}     ", end="", flush=True)

                self.current_algorithms[algorithm_name]["name_printed"] = True

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
        print(f"=" * 50)
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
            # Nova configuração (sempre reexibir)
            if level_id != self.current_config_id:
                self.current_config_id = level_id
                if data:
                    self.current_config_index = data.get("config_index", 1)
                    self.total_configs = data.get("total_configs", 1)

                    # Extrair nome da execução se disponível
                    execution_name = data.get("execution_name", level_id)

                    # Armazenar para uso posterior no dataset, mas não exibir aqui
                    self.pending_execution_name = execution_name
                    self.pending_config_index = self.current_config_index
                    self.pending_total_configs = self.total_configs

                # Reset dataset tracking quando nova configuração
                self.current_dataset_id = None

        elif level == ExecutionLevel.DATASET:
            # Sempre reexibir informações do dataset (incluindo execução)
            if level_id != self.current_dataset_id:
                self.current_dataset_id = level_id
                if data:
                    self.current_dataset_index = data.get("dataset_index", 1)
                    self.total_datasets = data.get("total_datasets", 1)
                    total_algorithms = data.get("total_algorithms", 0)

                    # Extrair dados de execução e dataset
                    execution_name = data.get("execution_name", "Execução")
                    config_index = data.get("config_index", 1)
                    total_configs = data.get("total_configs", 1)

                    dataset_name = data.get("dataset_name", level_id)
                    algorithm_config_name = data.get(
                        "algorithm_config_name", "Algoritmos"
                    )
                    algorithm_config_index = data.get("algorithm_config_index", 1)
                    total_algorithm_configs = data.get("total_algorithm_configs", 1)

                    # Salvar dados de execução para controle
                    self.current_execution_name = execution_name
                    self.current_config_index = config_index
                    self.total_configs = total_configs

                    # Exibir linha de execução apenas se for diferente da anterior
                    # Usar dados armazenados do nível EXECUTION se disponíveis
                    display_execution_name = getattr(
                        self, "pending_execution_name", execution_name
                    )
                    display_config_index = getattr(
                        self, "pending_config_index", config_index
                    )
                    display_total_configs = getattr(
                        self, "pending_total_configs", total_configs
                    )

                    # Criar nome de execução mais descritivo
                    exec_name_display = f"{display_execution_name} - {dataset_name}"

                    if exec_name_display != self.last_execution_name:
                        print(
                            f"\n\n⚡ Execução: {exec_name_display} ({display_config_index}/{display_total_configs})"
                        )
                        self.last_execution_name = exec_name_display

                    # Limpar dados pendentes
                    if hasattr(self, "pending_execution_name"):
                        delattr(self, "pending_execution_name")
                    if hasattr(self, "pending_config_index"):
                        delattr(self, "pending_config_index")
                    if hasattr(self, "pending_total_configs"):
                        delattr(self, "pending_total_configs")

                # Reset algoritmos quando novo dataset
                self.current_algorithms = {}
                self.algorithm_order = []
                self.algorithm_repetitions = {}

                # Exibir informações do dataset
                print(f"📊 Configuração do Algoritmo: {algorithm_config_name}")
                print(f"🗂️ Dataset: {dataset_name}")

                # Obter lista de algoritmos da configuração
                algorithms_config = data.get("algorithms_config", {})
                algorithm_names = algorithms_config.get("algorithms", [])

                # Listar todos os algoritmos da configuração
                for algo_name in algorithm_names:
                    print(f"• {algo_name}")
                    print(f"[░░░░░░░░░░░░░░░░░░░░] 0%")

                # Adicionar informações do tipo de task se for otimização
                if self.task_type == TaskType.OPTIMIZATION:
                    print("💡 Tipo: Otimização de Hiperparâmetros")
                elif self.task_type == TaskType.SENSITIVITY:
                    print("💡 Tipo: Análise de Sensibilidade")
                elif self.task_type == TaskType.EXECUTION:
                    print("💡 Tipo: Execução de Algoritmos")

    def update_item(
        self,
        item_id: str,
        progress: float,
        message: str = "",
        context: Optional[HierarchicalContext] = None,
    ) -> None:
        """Atualiza um item individual (algoritmo)."""
        # Usar sempre o algorithm_id do contexto como chave principal
        algorithm_name = (
            context.algorithm_id if context and context.algorithm_id else item_id
        )

        # Extrair ID de repetição do item_id se disponível
        rep_id = item_id

        # Primeiro algoritmo - adiciona à lista de ordem
        if algorithm_name not in self.current_algorithms:
            self.algorithm_order.append(algorithm_name)
            self.current_algorithms[algorithm_name] = {
                "progress": 0.0,
                "completed": False,
                "current_rep": 0,
                "total_rep": 1,
                "name_printed": False,
                "completed_reps": 0,
            }
            self.algorithm_repetitions[algorithm_name] = {}

        # Atualizar status
        algo_status = self.current_algorithms[algorithm_name]

        # Extrair informações de repetições se disponível no context
        current_rep = 1
        total_rep = 1
        if context and context.repetition_id:
            # Tentar extrair números de repetição do ID se for formatado como "1/10"
            rep_parts = context.repetition_id.split("/")
            if len(rep_parts) == 2:
                try:
                    current_rep = int(rep_parts[0])
                    total_rep = int(rep_parts[1])
                    algo_status["current_rep"] = current_rep
                    algo_status["total_rep"] = total_rep
                except ValueError:
                    pass

        # Atualizar progresso da repetição específica
        if rep_id not in self.algorithm_repetitions[algorithm_name]:
            self.algorithm_repetitions[algorithm_name][rep_id] = {
                "progress": 0.0,
                "completed": False,
                "rep_number": current_rep,
            }

        rep_data = self.algorithm_repetitions[algorithm_name][rep_id]
        old_rep_progress = rep_data["progress"]
        rep_data["progress"] = progress
        rep_data["completed"] = progress >= 100.0
        rep_data["rep_number"] = current_rep

        # Calcular progresso total do algoritmo
        completed_reps = sum(
            1
            for rep in self.algorithm_repetitions[algorithm_name].values()
            if rep["completed"]
        )

        # Progresso da repetição atual (se não completada)
        current_rep_progress = 0.0
        if not rep_data["completed"]:
            current_rep_progress = progress / 100.0

        # Progresso total = (repetições completas + progresso atual) / total repetições
        total_progress = (completed_reps + current_rep_progress) / total_rep * 100.0

        old_total_progress = algo_status["progress"]
        algo_status["progress"] = total_progress
        algo_status["completed"] = total_progress >= 100.0
        algo_status["completed_reps"] = completed_reps

        # Decidir se deve imprimir/atualizar
        significant_change = (
            abs(total_progress - old_total_progress) >= 10.0  # Mudança de 10% ou mais
            or algo_status["completed"]
            or (
                old_rep_progress < 100.0 and rep_data["completed"]
            )  # Nova repetição completada
        )

        if significant_change:
            # Criar linha de progresso
            progress_bar = self._create_progress_bar(total_progress)

            if algo_status["completed"]:
                status_line = f"• {algorithm_name} {progress_bar} ✅ 100%"
            else:
                # Mostrar progresso no formato: (finalizadas+em_progresso/total) porcentagem%
                if total_rep > 1:
                    if current_rep_progress > 0:
                        rep_info = f"({completed_reps}+{current_rep}/{total_rep}) "
                    else:
                        rep_info = f"({completed_reps}/{total_rep}) "
                else:
                    rep_info = ""

                # Adicionar informações do message se disponível (trials, etc.)
                message_info = ""
                if message and "Trial" in message:
                    message_info = f" | {message}"

                status_line = f"• {algorithm_name} {progress_bar} {rep_info}{total_progress:.0f}%{message_info}"

            # Limpar linha completamente e imprimir nova
            print(f"\r{' ' * 80}\r{status_line}", end="", flush=True)

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
        # Extrair o nome do algoritmo do item_id
        # item_id formato: "Baseline_teste_pequeno_1"
        # Queremos extrair "Baseline"
        algorithm_name = item_id.split("_")[0]

        # Marcar repetição como concluída se bem sucedida
        if success and algorithm_name in self.current_algorithms:
            # Marcar a repetição específica como concluída
            if algorithm_name in self.algorithm_repetitions:
                if item_id in self.algorithm_repetitions[algorithm_name]:
                    self.algorithm_repetitions[algorithm_name][item_id][
                        "completed"
                    ] = True
                    self.algorithm_repetitions[algorithm_name][item_id][
                        "progress"
                    ] = 100.0

                    # Recalcular progresso total do algoritmo
                    algo_status = self.current_algorithms[algorithm_name]
                    total_rep = algo_status["total_rep"]

                    completed_reps = sum(
                        1
                        for rep in self.algorithm_repetitions[algorithm_name].values()
                        if rep["completed"]
                    )

                    total_progress = (completed_reps / total_rep) * 100.0
                    algo_status["progress"] = total_progress
                    algo_status["completed"] = total_progress >= 100.0
                    algo_status["completed_reps"] = completed_reps

                    # Imprimir atualização se houve mudança significativa
                    progress_bar = self._create_progress_bar(total_progress)

                    if algo_status["completed"]:
                        status_line = f"• {algorithm_name} {progress_bar} ✅ 100%"
                        print(f"\r{' ' * 80}\r{status_line}")
                    else:
                        rep_info = (
                            f"({completed_reps}/{total_rep}) " if total_rep > 1 else ""
                        )
                        status_line = f"• {algorithm_name} {progress_bar} {rep_info}{total_progress:.0f}%"
                        print(f"\r{' ' * 80}\r{status_line}", end="", flush=True)

    def finish_task(
        self,
        success: bool = True,
        final_results: Optional[Dict[str, Any]] = None,
        error_message: str = "",
    ) -> None:
        """Finaliza monitoramento da tarefa."""
        if success:
            print(f"\n\n✅ Monitoramento concluído!")
        else:
            print(f"❌ Erro: {error_message}")

        if self.start_time:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            print(f"⏰ Tempo total: {elapsed/60:.0f}m {elapsed%60:.0f}s")
        print("=" * 50)

    def _create_progress_bar(self, progress: float) -> str:
        """Cria uma barra de progresso ASCII."""
        width = 20
        filled = int(width * progress / 100)
        empty = width - filled

        return f"[{'█' * filled}{'░' * empty}]"

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
