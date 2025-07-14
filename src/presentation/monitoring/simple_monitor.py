"""Monitor simples para terminal."""

import sys
import time
from datetime import datetime
from typing import Dict, Any

from .interfaces import MonitoringInterface, TaskType, MonitoringData


class SimpleMonitor(MonitoringInterface):
    """Monitor simples que exibe progresso no terminal."""
    
    def __init__(self):
        self.task_type: TaskType = None
        self.batch_name: str = ""
        self.start_time: datetime = None
        self.current_line_algorithm: str = ""
        self.header_printed: bool = False
        self.last_update_time: float = 0
        self.update_interval: float = 3.0  # 3 segundos
    
    def start_monitoring(self, task_type: TaskType, config: Dict[str, Any]) -> None:
        """Inicia o monitoramento."""
        self.task_type = task_type
        self.batch_name = config.get("batch_name", "Batch Sem Nome")
        self.start_time = datetime.now()
        self.header_printed = False
        self.current_line_algorithm = ""
        
    def update_progress(self, progress_data: MonitoringData) -> None:
        """Atualiza informações de progresso."""
        current_time = time.time()
        
        # Controle de frequência de atualização
        if current_time - self.last_update_time < self.update_interval:
            return
        
        self.last_update_time = current_time
        
        if not self.header_printed:
            self._print_header(progress_data)
            self.header_printed = True
            
        if self.task_type == TaskType.EXECUTION:
            self._update_execution_progress(progress_data.execution_data)
        elif self.task_type == TaskType.OPTIMIZATION:
            self._update_optimization_progress(progress_data.optimization_data)
        elif self.task_type == TaskType.SENSITIVITY:
            self._update_sensitivity_progress(progress_data.sensitivity_data)
    
    def _print_header(self, progress_data: MonitoringData) -> None:
        """Imprime cabeçalho do monitor."""
        if self.task_type == TaskType.EXECUTION:
            print("CSPBench - Monitoramento de Execução")
        elif self.task_type == TaskType.OPTIMIZATION:
            print("CSPBench - Monitoramento de Otimização")
        elif self.task_type == TaskType.SENSITIVITY:
            print("CSPBench - Análise de Sensibilidade")
            
        print("=" * 50)
        print(f"📋 Batch: {self.batch_name}")
        print(f"⏰ Iniciado: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print()
    
    def _update_execution_progress(self, data):
        """Atualiza progresso de execução."""
        if not data:
            return
            
        # Informações da execução atual
        exec_num = data.completed_executions + 1
        print(f"📊 Execução: {data.current_execution} ({exec_num}/{data.total_executions})")
        print(f"🗂️  Dataset: {data.current_dataset} ({data.current_dataset_index}/{data.total_datasets})")
        print(f"🧠 Algoritmos: {data.total_algorithms} total")
        print()
        
        # Barra de progresso geral
        if data.total_algorithms > 0:
            algo_progress = (data.completed_algorithms / data.total_algorithms) * 100
            progress_bar = self._create_progress_bar(algo_progress)
            print(f"Status Atual:")
            print(f"{progress_bar} {algo_progress:.0f}% ({data.completed_algorithms}/{data.total_algorithms})")
            print()
        
        # Algoritmos concluídos (linhas fixas)
        for algo, result in data.algorithm_results.items():
            if result.get('completed', False):
                distance = result.get('distance', 'N/A')
                exec_time = result.get('execution_time', 0)
                print(f"✅ {algo:<12} [{'█' * 20}] 100% - {exec_time:.2f}s - dist: {distance}")
        
        # Algoritmo atual (linha que será reescrita)
        if data.current_algorithm and not data.algorithm_results.get(data.current_algorithm, {}).get('completed', False):
            progress = data.algorithm_progress.get(data.current_algorithm, 0)
            callback_info = data.algorithm_callback_info.get(data.current_algorithm, "processando...")
            task_info = f" - {data.current_task_info}" if data.current_task_info else ""
            progress_bar = self._create_progress_bar(progress, width=20)
            
            # Limpa linha anterior se necessário
            if self.current_line_algorithm and self.current_line_algorithm != data.current_algorithm:
                print()
            
            line = f"⏳ {data.current_algorithm:<12} [{progress_bar}] {progress:.0f}%{task_info} - {callback_info}"
            print(f"\r{line}", end="", flush=True)
            self.current_line_algorithm = data.current_algorithm
        
        # Resultados parciais
        if data.best_distance is not None:
            print(f"\n\n📊 Resultados Parciais:")
            print(f"   Melhor distância: {data.best_distance}")
            elapsed = (datetime.now() - self.start_time).total_seconds()
            print(f"   Tempo total: {elapsed:.2f}s")
    
    def _update_optimization_progress(self, data):
        """Atualiza progresso de otimização."""
        if not data:
            return
            
        # Informações da otimização atual
        opt_num = data.completed_optimizations + 1
        print(f"📊 Otimização: {data.current_optimization} ({opt_num}/{data.total_optimizations})")
        print(f"🗂️  Dataset: {data.current_dataset} ({data.current_dataset_index}/{data.total_datasets})")
        print(f"🎯 Estudo: {data.study_name}")
        print()
        
        # Barra de progresso geral
        if data.n_trials > 0:
            trial_progress = (data.completed_trials / data.n_trials) * 100
            progress_bar = self._create_progress_bar(trial_progress)
            print(f"Progresso:")
            print(f"{progress_bar} {trial_progress:.0f}% ({data.completed_trials}/{data.n_trials})")
            print()
        
        # Trial atual (linha que será reescrita)
        if data.current_trial is not None:
            task_info = f" - {data.current_task_info}" if data.current_task_info else ""
            callback_info = data.trial_callback_info if data.trial_callback_info else ""
            
            # Limpa linha anterior se necessário
            if self.current_line_algorithm:
                print(f"\r", end="")
            
            line = f"⏳ Trial {data.current_trial}/{data.n_trials} executando...{task_info} {callback_info}"
            print(line, end="", flush=True)
            self.current_line_algorithm = f"trial_{data.current_trial}"
        
        # Resultados parciais
        if data.best_value is not None:
            print(f"\n\n🏆 Melhor até agora: {data.best_value}")
            elapsed = (datetime.now() - self.start_time).total_seconds()
            print(f"⏰ Tempo decorrido: {elapsed/60:.0f}m {elapsed%60:.0f}s")
    
    def _update_sensitivity_progress(self, data):
        """Atualiza progresso de análise de sensibilidade."""
        if not data:
            return
            
        # Informações da análise atual
        analysis_num = data.completed_analyses + 1
        print(f"📊 Análise: {data.current_analysis} ({analysis_num}/{data.total_analyses})")
        print(f"🗂️  Dataset: {data.current_dataset} ({data.current_dataset_index}/{data.total_datasets})")
        print(f"🔬 Método: {data.analysis_method} ({data.method_details})")
        print(f"📊 Parâmetros: {len(data.parameters)} ({', '.join(data.parameters)})")
        print()
        
        # Barra de progresso geral
        if data.n_samples > 0:
            sample_progress = (data.completed_samples / data.n_samples) * 100
            progress_bar = self._create_progress_bar(sample_progress)
            print(f"Progresso Geral:")
            print(f"{progress_bar} {sample_progress:.0f}% ({data.completed_samples}/{data.n_samples} amostras)")
            print()
        
        # Parâmetros analisados
        print("Parâmetros Analisados:")
        for param in data.parameters:
            if param in data.sensitivity_results:
                # Parâmetro concluído
                results = data.sensitivity_results[param]
                mu_star = results.get('mu_star', 0)
                sigma = results.get('sigma', 0)
                print(f"✅ {param:<12} [{'█' * 20}] 100% - μ*={mu_star:.3f}, σ={sigma:.3f}")
            elif param == data.current_parameter:
                # Parâmetro atual
                progress = data.parameter_progress.get(param, 0)
                progress_bar = self._create_progress_bar(progress, width=20)
                task_info = f" - {data.current_task_info}" if data.current_task_info else ""
                callback_info = data.callback_info if data.callback_info else "processando..."
                
                # Limpa linha anterior se necessário
                if self.current_line_algorithm and self.current_line_algorithm != param:
                    print()
                
                line = f"⏳ {param:<12} [{progress_bar}] {progress:.0f}%{task_info} - {callback_info}"
                print(f"\r{line}", end="", flush=True)
                self.current_line_algorithm = param
        
        # Resultados preliminares
        if data.sensitivity_results:
            print(f"\n\n🔍 Sensibilidade Preliminar:")
            sorted_params = sorted(data.sensitivity_results.items(), 
                                   key=lambda x: x[1].get('mu_star', 0), reverse=True)
            for i, (param, results) in enumerate(sorted_params[:3], 1):
                mu_star = results.get('mu_star', 0)
                if mu_star > 0.2:
                    level = "ALTA"
                elif mu_star > 0.1:
                    level = "MÉDIA"
                else:
                    level = "BAIXA"
                print(f"   {i}. {param}: {level} (maior impacto)")
        
        elapsed = (datetime.now() - self.start_time).total_seconds()
        print(f"\n⏰ Tempo: {elapsed/60:.0f}m {elapsed%60:.0f}s")
    
    def _create_progress_bar(self, progress: float, width: int = 30) -> str:
        """Cria uma barra de progresso ASCII."""
        filled_length = int(width * progress // 100)
        bar = '█' * filled_length + '░' * (width - filled_length)
        return f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" if width == 30 else f"[{bar}]"
    
    def finish_monitoring(self, results: Dict[str, Any]) -> None:
        """Finaliza o monitoramento."""
        if self.current_line_algorithm:
            print()  # Nova linha final
        
        print("\n✅ Monitoramento concluído!")
        elapsed = (datetime.now() - self.start_time).total_seconds()
        print(f"⏰ Tempo total: {elapsed/60:.0f}m {elapsed%60:.0f}s")
        print("=" * 50)
    
    def show_error(self, error: str) -> None:
        """Exibe erro."""
        print(f"\n❌ Erro: {error}")
    
    def close(self) -> None:
        """Fecha e limpa recursos do monitor."""
        if self.current_line_algorithm:
            print()  # Garante nova linha final
