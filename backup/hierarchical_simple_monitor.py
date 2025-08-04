#!/usr/bin/env python3
"""SimpleMonitor Hierárquico - ETAPA 4"""

from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional
from threading import RLock

from src.presentation.monitoring.interfaces import (
    ExecutionLevel,
    HierarchicalContext,
    MonitoringInterface,
    TaskType,
    ExecutionHierarchy,
    ActiveRun,
    CallbackEntry,
    ItemStatus
)                print("=" * 60)

    def _extract_number_from_context(self, context_id: str, default: int = 1) -> int:erarchicalSimpleMonitor(MonitoringInterface):
    """Monitor hierárquico que exibe progresso estruturado no terminal."""

    def __init__(self):
        # Estado da tarefa
        self.task_type: Optional[TaskType] = None
        self.task_name: str = ""
        self.start_time: Optional[datetime] = None
        self.is_active: bool = False
        
        # Estrutura hierárquica
        self.hierarchy: Optional[ExecutionHierarchy] = None
        
        # Runs ativas
        self.active_runs: Dict[str, ActiveRun] = {}
        
        # Callbacks recentes
        self.callbacks: List[CallbackEntry] = []
        self.max_callbacks: int = 50
        
        # Controle de exibição
        self.header_printed: bool = False
        self.last_display_time: datetime = datetime.now()
        self.display_interval: float = 1.0  # Atualizar exibição a cada 1 segundo
        
        # Controle de estado para evitar repetições
        self.last_displayed_context: Optional[tuple] = None
        
        # Thread safety
        self._lock: RLock = RLock()

    def initialize_hierarchy(self, hierarchy: ExecutionHierarchy) -> None:
        """Inicializa estrutura hierárquica."""
        with self._lock:
            self.hierarchy = hierarchy
            
    def start_task(
        self, task_type: TaskType, task_name: str, config: Dict[str, Any]
    ) -> None:
        """Inicia monitoramento de uma tarefa."""
        with self._lock:
            self.task_type = task_type
            self.task_name = task_name
            self.start_time = datetime.now()
            self.is_active = True
            self.header_printed = False
            
            # Tentar extrair e inicializar hierarquia da configuração
            if not self.hierarchy:
                try:
                    from src.presentation.monitoring.hierarchy_parser import HierarchyParser
                    hierarchy = HierarchyParser.parse_yaml_hierarchy(config)
                    self.hierarchy = hierarchy
                except ImportError:
                    # Se não conseguir importar o parser, criar hierarquia básica
                    self._create_basic_hierarchy(config)
                except Exception:
                    # Se não conseguir fazer parse, criar hierarquia básica
                    self._create_basic_hierarchy(config)
            
        self._print_header()

    def update_hierarchy_level(
        self, 
        level: ExecutionLevel, 
        current: int, 
        total: int, 
        name: str = "",
        metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """Atualiza nível específico da hierarquia."""
        with self._lock:
            if self.hierarchy:
                self.hierarchy.update_level(level, current, total, name)
                self._maybe_update_display()

    def start_run(self, run: ActiveRun) -> None:
        """Inicia uma run individual."""
        with self._lock:
            run.start_time = datetime.now()
            run.status = ItemStatus.RUNNING
            self.active_runs[run.run_id] = run
            self._maybe_update_display()

    def update_run_progress(
        self, 
        run_id: str, 
        progress: float, 
        message: str = "",
        generation: Optional[int] = None
    ) -> None:
        """Atualiza progresso de uma run específica."""
        with self._lock:
            if run_id in self.active_runs:
                run = self.active_runs[run_id]
                run.update_progress(progress, message)
                if generation is not None:
                    run.current_generation = generation
                
                self._maybe_update_display()

    def finish_run(
        self, 
        run_id: str, 
        success: bool = True, 
        result: Optional[Dict[str, Any]] = None,
        error: Optional[str] = None
    ) -> None:
        """Finaliza uma run específica."""
        with self._lock:
            if run_id in self.active_runs:
                run = self.active_runs[run_id]
                run.status = ItemStatus.COMPLETED if success else ItemStatus.FAILED
                run.progress = 100.0 if success else run.progress
                print(f"DEBUG: Finished run {run_id}: {run.algorithm_name}, success: {success}")
                
                # Atualizar exibição uma última vez
                self._maybe_update_display()
                
                # Remover run finalizada após exibição para evitar acúmulo
                # Aguardar um pouco para garantir que foi exibida
                import threading
                def remove_run():
                    import time
                    time.sleep(2)  # Aguardar 2 segundos
                    with self._lock:
                        if run_id in self.active_runs:
                            print(f"DEBUG: Removing completed run {run_id}")
                            del self.active_runs[run_id]
                
                threading.Thread(target=remove_run, daemon=True).start()
                
            else:
                print(f"DEBUG: Tried to finish run {run_id} but it's not in active_runs")

    def algorithm_callback(
        self,
        algorithm_name: str,
        progress: float,
        message: str,
        run_id: Optional[str] = None,
        generation: Optional[int] = None,
        context: Optional[HierarchicalContext] = None
    ) -> None:
        """Callback direto do algoritmo durante execução."""
        with self._lock:
            # Adicionar callback à lista
            callback = CallbackEntry(
                timestamp=datetime.now(),
                algorithm_name=algorithm_name,
                run_id=run_id,
                message=message,
                progress=progress,
                generation=generation
            )
            
            self.callbacks.append(callback)
            
            # Manter apenas os mais recentes
            if len(self.callbacks) > self.max_callbacks:
                self.callbacks = self.callbacks[-self.max_callbacks:]
            
            # Atualizar run se especificada
            if run_id and run_id in self.active_runs:
                self.update_run_progress(run_id, progress, message, generation)
            else:
                self._maybe_update_display()

    def get_active_runs(self) -> List[ActiveRun]:
        """Retorna lista de runs ativas."""
        with self._lock:
            return list(self.active_runs.values())

    def get_recent_callbacks(self, limit: int = 10) -> List[CallbackEntry]:
        """Retorna callbacks recentes."""
        with self._lock:
            return self.callbacks[-limit:] if self.callbacks else []

    def get_hierarchy_status(self) -> Dict[str, Any]:
        """Retorna status atual da hierarquia."""
        with self._lock:
            if not self.hierarchy:
                return {}
            return self.hierarchy.to_dict()

    def get_full_status(self) -> Dict[str, Any]:
        """Retorna status completo."""
        with self._lock:
            return {
                "hierarchy": self.get_hierarchy_status(),
                "active_runs": [run.to_dict() for run in self.active_runs.values()],
                "recent_callbacks": [cb.to_dict() for cb in self.get_recent_callbacks()],
                "overall_progress": self.hierarchy.get_overall_progress() if self.hierarchy else 0.0,
                "task_info": {
                    "type": self.task_type.value if self.task_type else "unknown",
                    "name": self.task_name,
                    "start_time": self.start_time.isoformat() if self.start_time else None,
                    "elapsed_time": (datetime.now() - self.start_time).total_seconds() if self.start_time else 0
                }
            }

    def _maybe_update_display(self) -> None:
        """Atualiza exibição se passou tempo suficiente."""
        now = datetime.now()
        if (now - self.last_display_time).total_seconds() >= self.display_interval:
            self._update_display()
            self.last_display_time = now

    def _print_header(self) -> None:
        """Imprime cabeçalho do monitor."""
        if self.header_printed:
            return
            
        print("=" * 60)
        if self.task_type == TaskType.EXECUTION:
            print("🚀 CSPBench - Execution Monitoring")
        elif self.task_type == TaskType.OPTIMIZATION:
            print("🔍 CSPBench - Optimization Monitoring") 
        elif self.task_type == TaskType.SENSITIVITY:
            print("📊 CSPBench - Sensitivity Analysis")
        else:
            print("⚙️  CSPBench - Task Monitoring")
            
        print("=" * 60)
        print(f"📋 Batch: {self.task_name}")
        if self.start_time:
            print(f"⏰ Started: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 60)
        self.header_printed = True

    def _update_display(self) -> None:
        """Atualiza exibição completa."""
        if not self.is_active:
            return
            
        # Exibir no novo formato organizado
        self._display_hierarchical_organized()

    def _display_hierarchical_organized(self) -> None:
        """Exibe status de forma organizada por execução -> dataset -> configuração."""
        if not self.hierarchy or not self.active_runs:
            print("DEBUG: No hierarchy or active runs")
            return
        
        print(f"DEBUG: Total active runs: {len(self.active_runs)}")
        
        # Organizar runs por contexto hierárquico - VERSÃO SIMPLIFICADA
        runs_by_context = {}
        
        for run_id, run in self.active_runs.items():
            if not run.context:
                continue
                
            # Criar chave hierárquica
            exec_key = run.context.execution_id or "Unknown Execution"
            dataset_key = run.context.dataset_id or "Unknown Dataset"
            config_key = run.context.config_id or "Unknown Configuration"
            
            context_key = f"{exec_key}|{dataset_key}|{config_key}"
            
            if context_key not in runs_by_context:
                runs_by_context[context_key] = {
                    'exec_id': exec_key,
                    'dataset_id': dataset_key, 
                    'config_id': config_key,
                    'runs': []
                }
                
            runs_by_context[context_key]['runs'].append(run)
        
        # Exibir todos os contextos que têm atividade
        for context_key, context_data in runs_by_context.items():
            exec_id = context_data['exec_id']
            dataset_id = context_data['dataset_id']
            config_id = context_data['config_id']
            runs = context_data['runs']
            
            # Verificar se há atividade neste contexto
            has_active = any(run.status == ItemStatus.RUNNING for run in runs)
            has_completed = any(run.status == ItemStatus.COMPLETED for run in runs)
            
            if has_active or has_completed:
                # Calcular numeração
                exec_num = 1 if exec_id == "Test Execution" else 2
                total_execs = 2
                
                dataset_name = self._get_friendly_name(dataset_id, "Dataset")
                dataset_num = 1 if dataset_id == "dataset_test" else (2 if dataset_id == "dataset_file" else 1)
                total_datasets = 2 if exec_id == "Test Execution" else 1
                
                config_name = self._get_friendly_name(config_id, "Configuration")
                config_num = 1 if config_id == "default_config" else 2
                total_configs = 2 if exec_id == "Test Execution" else 1
                
                # Cabeçalho
                print(f"Execution: {exec_id} ({exec_num}/{total_execs})")
                print(f"Datasets: {dataset_name} ({dataset_num}/{total_datasets})")
                print(f"Configuration: {config_name} ({config_num}/{total_configs})")
                print("-" * 40)
                
                # Exibir runs
                for run in sorted(runs, key=lambda r: r.algorithm_name):
                    progress_bar = self._create_progress_bar(run.progress)
                    status_icon = "✅" if run.status == ItemStatus.COMPLETED else ("🔄" if run.status == ItemStatus.RUNNING else "❌")
                    print(f"  {status_icon} {run.algorithm_name} (Run {run.run_info}): {progress_bar} {run.progress:.1f}%")
                
                print("=" * 60)
            # Calcular numeração baseada no TEMPLATE.yaml
            if exec_id == "Test Execution":
                exec_num = 1
            elif exec_id == "Complete Execution":
                exec_num = 2
            else:
                exec_num = 1
            total_execs = 2  # From TEMPLATE.yaml: "Test Execution" + "Complete Execution"
            
            dataset_name = self._get_friendly_name(dataset_id, "Dataset")
            
            # Calcular número do dataset e total baseado na execução atual
            if exec_id == "Test Execution":
                # Test Execution has 2 datasets: dataset_test, dataset_file
                if dataset_id == "dataset_test":
                    dataset_num = 1
                elif dataset_id == "dataset_file":
                    dataset_num = 2
                else:
                    dataset_num = 1
                total_datasets = 2
            elif exec_id == "Complete Execution":
                # Complete Execution has 1 dataset: dataset_ncbi
                dataset_num = 1
                total_datasets = 1
            else:
                dataset_num = 1
                total_datasets = 1
            
            config_name = self._get_friendly_name(config_id, "Configuration")
            
            # Calcular número da configuração e total baseado na execução atual
            if exec_id == "Test Execution":
                # Test Execution has 2 configs: default_config, aggressive_csc
                if config_id == "default_config":
                    config_num = 1
                elif config_id == "aggressive_csc":
                    config_num = 2
                else:
                    config_num = 1
                total_configs = 2
            elif exec_id == "Complete Execution":
                # Complete Execution has 1 config: aggressive_csc
                config_num = 1
                total_configs = 1
            else:
                config_num = 1
                total_configs = 1
            
            # Cabeçalho da seção
            print(f"Execution: {exec_id} ({exec_num}/{total_execs})")
            print(f"Datasets: {dataset_name} ({dataset_num}/{total_datasets})")
            print(f"Configuration: {config_name} ({config_num}/{total_configs})")
            print("-" * 40)
            
            # Ordenar runs por nome do algoritmo para exibição consistente
            runs_sorted = sorted(runs, key=lambda r: r.algorithm_name)
            
            # Exibir runs desta configuração
            for run in runs_sorted:
                progress_bar = self._create_progress_bar(run.progress)
                
                if run.status == ItemStatus.COMPLETED:
                    status_icon = "✅"
                elif run.status == ItemStatus.FAILED:
                    status_icon = "❌" 
                            elif run.status == ItemStatus.RUNNING:
                                status_icon = "🔄"
                            else:
                                status_icon = "⏸️"
                                
                            print(f"  {status_icon} {run.algorithm_name} (Run {run.run_info}): {progress_bar} {run.progress:.1f}%")
                        
                        print("=" * 60)

    def _extract_number_from_context(self, context_id: str, default: int = 1) -> int:
        """Extrai número do contexto (ex: 'execution_1' -> 1)."""
        try:
            if '_' in context_id:
                return int(context_id.split('_')[-1])
            return default
        except:
            return default
    
    def _get_friendly_name(self, context_id: str, default_prefix: str = "Item") -> str:
        """Converte ID em nome amigável."""
        name_mapping = {
            # Executions 
            "Test Execution": "Test Execution",
            "Complete Execution": "Complete Execution",
            
            # Datasets
            "dataset_test": "Test Dataset",
            "dataset_file": "File Dataset", 
            "dataset_ncbi": "NCBI Dataset",
            
            # Configurations
            "default_config": "Default Configuration",
            "aggressive_csc": "Aggressive CSC Configuration"
        }
        
        return name_mapping.get(context_id, context_id.replace('_', ' ').title())
    
    def _get_total_executions(self) -> int:
        """Retorna total de execuções."""
        if self.hierarchy:
            status = self.hierarchy.to_dict()
            return status.get("execution", {}).get("total", 1)
        return 1
    
    def _get_total_datasets(self) -> int:
        """Retorna total de datasets.""" 
        if self.hierarchy:
            status = self.hierarchy.to_dict()
            return status.get("dataset", {}).get("total", 1)
        return 1
    
    def _get_total_configs(self) -> int:
        """Retorna total de configurações."""
        if self.hierarchy:
            status = self.hierarchy.to_dict()
            return status.get("config", {}).get("total", 1)
        return 1

    def _display_hierarchy(self) -> None:
        """Exibe status hierárquico atual."""
        if not self.hierarchy:
            return
            
        status = self.hierarchy.to_dict()
        overall = self.hierarchy.get_overall_progress()
        
        print(f"\\n📊 Overall Progress: {overall:.1f}%")
        print("=" * 40)
        
        # Exibir cada nível da hierarquia
        exec_info = status.get("execution", {})
        if exec_info.get("total", 0) > 0:
            print(f"🔄 Executions: {exec_info['current']}/{exec_info['total']} ({exec_info['name']})")
            
        dataset_info = status.get("dataset", {})
        if dataset_info.get("total", 0) > 0:
            print(f"📁 Datasets: {dataset_info['current']}/{dataset_info['total']} ({dataset_info['name']})")
            
        config_info = status.get("config", {})
        if config_info.get("total", 0) > 0:
            print(f"⚙️  Configs: {config_info['current']}/{config_info['total']} ({config_info['name']})")
            
        algo_info = status.get("algorithm", {})
        if algo_info.get("total", 0) > 0:
            print(f"🧠 Algorithms: {algo_info['current']}/{algo_info['total']} ({algo_info['name']})")

    def _display_active_runs(self) -> None:
        """Exibe runs ativas."""
        if not self.active_runs:
            return
            
        print("\\n🏃 Active Runs:")
        print("-" * 40)
        
        for run in self.active_runs.values():
            # Criar barra de progresso
            progress_bar = self._create_progress_bar(run.progress)
            
            # Status icon
            if run.status == ItemStatus.COMPLETED:
                status_icon = "✅"
            elif run.status == ItemStatus.FAILED:
                status_icon = "❌"
            elif run.status == ItemStatus.RUNNING:
                status_icon = "🔄"
            else:
                status_icon = "⏸️"
                
            # Informações da run
            run_line = f"  {status_icon} {run.algorithm_name} (Run {run.run_info}): {progress_bar} {run.progress:.1f}%"
            
            # Adicionar geração se disponível
            if run.current_generation > 0:
                run_line += f" | Gen {run.current_generation}"
                
            # Adicionar mensagem se disponível
            if run.current_message:
                run_line += f" | {run.current_message}"
                
            print(run_line)

    def _display_recent_callbacks(self) -> None:
        """Exibe callbacks recentes."""
        recent = self.get_recent_callbacks(5)  # Últimos 5
        if not recent:
            return
            
        print("\\n💬 Recent Updates:")
        print("-" * 40)
        
        for callback in reversed(recent):  # Mais recente primeiro
            timestamp = callback.timestamp.strftime("%H:%M:%S")
            line = f"  [{timestamp}] {callback.algorithm_name}"
            
            if callback.generation and callback.generation > 0:
                line += f" (Gen {callback.generation})"
                
            line += f": {callback.message}"
            
            if callback.progress > 0:
                line += f" ({callback.progress:.1f}%)"
                
            print(line)

    def _create_progress_bar(self, progress: float, width: int = 10) -> str:
        """Cria barra de progresso ASCII."""
        filled = int(width * progress / 100)
        empty = width - filled
        return f"[{'█' * filled}{'░' * empty}]"

    # Métodos legados (compatibilidade)
    def start_item(self, item_id: str, item_type: str = "repetition", 
                   context: Optional[HierarchicalContext] = None,
                   metadata: Optional[Dict[str, Any]] = None) -> None:
        """Método legado - mapeia para start_run."""
        algorithm_name = context.algorithm_id if context and context.algorithm_id else item_id.split("_")[0]
        run_info = context.repetition_id if context and context.repetition_id else "1/1"
        
        run = ActiveRun(
            run_id=item_id,
            algorithm_name=algorithm_name,
            run_info=run_info,
            context=context,
            metadata=metadata or {}
        )
        self.start_run(run)

    def update_item(self, item_id: str, progress: float, message: str = "",
                    context: Optional[HierarchicalContext] = None) -> None:
        """Método legado - mapeia para update_run_progress."""
        self.update_run_progress(item_id, progress, message)

    def finish_item(self, item_id: str, success: bool = True,
                    result: Optional[Dict[str, Any]] = None,
                    error: Optional[str] = None) -> None:
        """Método legado - mapeia para finish_run."""
        self.finish_run(item_id, success, result, error)

    def update_hierarchy(self, level: ExecutionLevel, level_id: str, 
                        progress: float, message: str = "",
                        data: Optional[Dict[str, Any]] = None) -> None:
        """Método legado - funcionalidade simplificada."""
        # Este método é mantido para compatibilidade mas com funcionalidade limitada
        pass

    def _create_basic_hierarchy(self, config: Dict[str, Any]) -> None:
        """Cria hierarquia básica baseada na configuração do batch."""
        try:
            from src.presentation.monitoring.interfaces import ExecutionHierarchy, ExecutionLevel
            hierarchy = ExecutionHierarchy()
            
            # Analisar estrutura da configuração para determinar a hierarquia
            
            # Verificar se é batch estruturado
            if "execution" in config and "executions" in config.get("execution", {}):
                # Batch estruturado novo formato
                execution_config = config["execution"]
                executions = execution_config.get("executions", [])
                
                # Nível de execução
                hierarchy.update_level(ExecutionLevel.EXECUTION, 1, len(executions), 
                                     executions[0].get("name", "execution_1") if executions else "execution")
                
                # Analisar primeiro execution para determinar estrutura
                if executions:
                    first_exec = executions[0]
                    datasets = first_exec.get("datasets", [])
                    configs = first_exec.get("configs", [])
                    algorithms = first_exec.get("algorithms", [])
                    
                    if datasets:
                        hierarchy.update_level(ExecutionLevel.DATASET, 1, len(datasets), 
                                             datasets[0] if isinstance(datasets[0], str) else datasets[0].get("name", "dataset_1"))
                    
                    if configs:
                        hierarchy.update_level(ExecutionLevel.CONFIG, 1, len(configs), 
                                             configs[0] if isinstance(configs[0], str) else configs[0].get("name", "config_1"))
                    
                    if algorithms:
                        hierarchy.update_level(ExecutionLevel.ALGORITHM, 1, len(algorithms), 
                                             algorithms[0] if isinstance(algorithms[0], str) else algorithms[0].get("name", "algorithm_1"))
                        
            elif "experiments" in config:
                # Batch legado com experiments
                experiments = config["experiments"]
                
                # Estimar hierarquia baseada nos experimentos
                hierarchy.update_level(ExecutionLevel.EXECUTION, 1, len(experiments), "experiments")
                
                if experiments:
                    first_exp = experiments[0]
                    datasets = first_exp.get("datasets", [])
                    algorithms = first_exp.get("algorithms", [])
                    
                    if datasets:
                        hierarchy.update_level(ExecutionLevel.DATASET, 1, len(datasets), datasets[0])
                    
                    if algorithms:  
                        hierarchy.update_level(ExecutionLevel.ALGORITHM, 1, len(algorithms), algorithms[0])
                        
            else:
                # Estrutura simples ou compatibilidade
                datasets = config.get("datasets", [])
                algorithms = config.get("algorithms", [])
                
                if datasets:
                    hierarchy.update_level(ExecutionLevel.DATASET, 1, len(datasets), datasets[0])
                
                if algorithms:
                    hierarchy.update_level(ExecutionLevel.ALGORITHM, 1, len(algorithms), algorithms[0])
            
            self.hierarchy = hierarchy
            
        except Exception as e:
            # Se algo der errado, criar hierarquia mínima
            from src.presentation.monitoring.interfaces import ExecutionHierarchy, ExecutionLevel
            hierarchy = ExecutionHierarchy()
            hierarchy.update_level(ExecutionLevel.EXECUTION, 1, 1, "batch_execution")
            self.hierarchy = hierarchy

    def finish_task(self, success: bool = True, 
                   final_results: Optional[Dict[str, Any]] = None,
                   error_message: str = "") -> None:
        """Finaliza monitoramento da tarefa."""
        with self._lock:
            self.is_active = False
            
        if success:
            print("✅ Task completed successfully!")
        else:
            print(f"❌ Task failed: {error_message}")

        if self.start_time:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            print(f"⏰ Total time: {elapsed/60:.0f}m {elapsed%60:.0f}s")
            
        print("=" * 60)
        
        # Exibir resumo final limpo
        if final_results:
            self._display_final_summary(final_results)

    def _display_final_summary(self, results: Dict[str, Any]) -> None:
        """Exibe resumo final limpo."""
        successful = results.get('successful', 0)
        failed = results.get('failed', 0)
        
        print(f"✅ Batch completed: {successful} successful, {failed} failed")

    def get_summary(self) -> Dict[str, Any]:
        """Retorna resumo consolidado."""
        return self.get_full_status()

    def show_error(self, error: str) -> None:
        """Exibe erro ao usuário."""
        print(f"\\n❌ Error: {error}")

    def stop(self) -> None:
        """Para monitoramento."""
        with self._lock:
            self.is_active = False

    def close(self) -> None:
        """Fecha sistema de monitoramento."""
        self.stop()


if __name__ == "__main__":
    # Teste básico
    print("🧪 Testando HierarchicalSimpleMonitor...")
    
    monitor = HierarchicalSimpleMonitor()
    
    # Simular uso
    monitor.start_task(TaskType.EXECUTION, "Test Hierarchical Monitor", {})
    
    # Criar hierarquia
    from src.presentation.monitoring.interfaces import ExecutionHierarchy
    hierarchy = ExecutionHierarchy()
    hierarchy.update_level(ExecutionLevel.EXECUTION, 1, 2, "Test Execution")
    hierarchy.update_level(ExecutionLevel.DATASET, 1, 2, "dataset_test")
    hierarchy.update_level(ExecutionLevel.CONFIG, 1, 2, "default_config")
    hierarchy.update_level(ExecutionLevel.ALGORITHM, 1, 5, "Baseline")
    
    monitor.initialize_hierarchy(hierarchy)
    
    # Simular run
    from src.presentation.monitoring.interfaces import ActiveRun, ItemStatus, HierarchicalContext
    context = HierarchicalContext(execution_id="Test Execution", dataset_id="dataset_test", 
                                config_id="default_config", algorithm_id="Baseline", repetition_id="1/3")
    
    run = ActiveRun(run_id="test_run", algorithm_name="Baseline", run_info="1/3", context=context)
    monitor.start_run(run)
    
    # Simular progresso
    import time
    for i in range(0, 101, 20):
        monitor.update_run_progress("test_run", i, f"Processing step {i//20 + 1}/5")
        monitor.algorithm_callback("Baseline", i, f"Callback: step {i//20 + 1} completed", "test_run")
        time.sleep(0.5)
    
    monitor.finish_run("test_run", True)
    monitor.finish_task(True)
    
    print("✅ Teste básico concluído!")
