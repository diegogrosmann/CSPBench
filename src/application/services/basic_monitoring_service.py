"""Basic monitoring service compatible with the new interface."""

from typing import Any, Dict, Optional

from src.infrastructure.logging_config import get_logger
from src.presentation.monitoring.interfaces import TaskType
from src.presentation.monitoring.monitor_factory import MonitorFactory


class BasicMonitoringService:
    """Basic monitoring service for integration with execution system."""

    def __init__(self, config: Dict[str, Any], web_session_manager=None):
        """
        Initialize monitoring service.

        Args:
            config: Batch/system configuration
            web_session_manager: Web session manager for progress updates
        """
        self.config = config
        self.logger = get_logger(__name__)
        self.monitor = None
        self.is_active = False
        self.web_session_manager = web_session_manager
        self.web_session_id = None
        
        # Debug log
        self.logger.info(f"BasicMonitoringService initialized with web_session_manager: {web_session_manager is not None}")

    def set_web_session_id(self, session_id: str) -> None:
        """Set the web session ID for logging."""
        self.web_session_id = session_id
        self.logger.info(f"Web session ID set to: {session_id}")

    def start_monitoring(
        self,
        task_type: TaskType,
        batch_name: str,
        batch_config: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Start monitoring.

        Args:
            task_type: Task type (EXECUTION, OPTIMIZATION, SENSITIVITY)
            batch_name: Batch name
            batch_config: Batch configuration (optional)
        """
        try:
            self.logger.info(f"Starting monitoring for task: {task_type.value}, batch: {batch_name}")
            
            # Create monitor based on configuration
            self.monitor = MonitorFactory.create_monitor(self.config)

            if not self.monitor:
                self.logger.info("Monitoring disabled in configuration")
                return

            self.logger.info(f"Monitor created: {type(self.monitor).__name__}")

            # Start monitoring in interface
            self.monitor.start_task(task_type, batch_name, batch_config or {})
            self.is_active = True

            self.logger.info(f"Monitoring started for task: {task_type.value}")
            
            # Enviar log inicial para a sessão web
            if self.web_session_manager and self.web_session_id:
                self._send_general_log("INFO", f"Starting {task_type.value}: {batch_name}")

        except Exception as e:
            self.logger.error(f"Error starting monitoring: {e}")
            self.monitor = None

    def show_error(self, error: str) -> None:
        """
        Display error in monitor.

        Args:
            error: Error message
        """
        if self.monitor:
            try:
                self.monitor.show_error(error)
            except Exception as e:
                self.logger.error(f"Error displaying error in monitor: {e}")
        
        # Enviar log de erro para a sessão web
        if self.web_session_manager and self.web_session_id:
            self._send_general_log("ERROR", f"Execution error: {error}")

    def finish_monitoring(self, results: Optional[Dict[str, Any]] = None) -> None:
        """
        Finish monitoring.

        Args:
            results: Final results (optional)
        """
        if not self.is_active:
            return

        try:
            if self.monitor:
                success = bool(results is not None and results.get("results"))
                self.monitor.finish_task(success=success, final_results=results)

            self.is_active = False
            self.logger.info("Monitoring finished")
            
            # Enviar log final para a sessão web
            if self.web_session_manager and self.web_session_id:
                if results and results.get("results"):
                    self._send_general_log("SUCCESS", "Batch execution completed successfully")
                else:
                    self._send_general_log("WARNING", "Batch execution finished without results")

        except Exception as e:
            self.logger.error(f"Error finishing monitoring: {e}")

    def close(self) -> None:
        """Close monitoring service."""
        if self.monitor:
            try:
                self.monitor.close()
            except Exception as e:
                self.logger.error(f"Error closing monitor: {e}")

    def set_web_session(self, session_id: str) -> None:
        """Set web session ID for progress updates."""
        self.web_session_id = session_id

    def _update_web_session_progress(self, progress: Dict[str, Any]) -> None:
        """Update web session with progress information."""
        if self.web_session_manager and self.web_session_id:
            try:
                self.web_session_manager.update_session(self.web_session_id, {
                    "progress": progress
                })
            except Exception as e:
                self.logger.error(f"Error updating web session progress: {e}")

    def _map_to_hierarchical_structure(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Map internal execution data to hierarchical display structure.
        
        Internal structure uses 'configs' as the basic unit, but UI needs to show:
        - Executions (execution blocks)
        - Datasets (within current execution)  
        - Algorithms (within current config)
        - Runs (repetitions within current dataset+algorithm combination)
        """
        # Default values
        progress_data = {
            "current_execution": 1,
            "total_executions": 1,
            "current_dataset": 1,
            "total_datasets": 1,
            "current_algorithm": 1,
            "total_algorithms": 1,
            "current_run": 1,
            "total_runs": 1,
        }
        
        # Map from internal data if available
        if data:
            # Execution level (maps from config_index/total_configs)
            progress_data["current_execution"] = data.get("execution_index", data.get("config_index", 1))
            progress_data["total_executions"] = data.get("total_executions", data.get("total_configs", 1))
            
            # Dataset level
            progress_data["current_dataset"] = data.get("dataset_index", 1)
            progress_data["total_datasets"] = data.get("total_datasets", 1)
            
            # Algorithm level (maps from algorithm_config_index)
            progress_data["current_algorithm"] = data.get("algorithm_index", data.get("algorithm_config_index", 1))
            progress_data["total_algorithms"] = data.get("total_algorithms", 1)
            
            # Run level (repetitions)
            progress_data["current_run"] = data.get("current_repetition", data.get("current_run", 1))
            progress_data["total_runs"] = data.get("total_repetitions", data.get("total_runs", 1))
        
        return progress_data

    def _calculate_overall_progress(self, details: Dict[str, Any]) -> int:
        """Calculate overall progress percentage from execution details."""
        try:
            current_exec = details.get("current_execution_index", 1)
            total_exec = details.get("total_executions", 1)
            
            current_dataset = details.get("current_dataset_index", 1) 
            total_datasets = details.get("total_datasets", 1)
            
            current_algo = details.get("current_algorithm_index", 1)
            total_algos = details.get("total_algorithms", 1)
            
            current_run = details.get("current_run", 1)
            total_runs = details.get("total_runs", 1)
            
            # Calculate progress based on total items processed
            if total_runs > 1:
                return int((current_run / total_runs) * 100)
            elif total_algos > 1:
                return int((current_algo / total_algos) * 100)  
            elif total_datasets > 1:
                return int((current_dataset / total_datasets) * 100)
            elif total_exec > 1:
                return int((current_exec / total_exec) * 100)
            else:
                return 0
        except Exception:
            return 0

    def _send_algorithm_log(self, algorithm_name: str, progress: float, message: str) -> None:
        """Send algorithm log to web session."""
        if not self.web_session_manager or not self.web_session_id:
            return
            
        try:
            # Criar mensagem formatada para o algoritmo
            formatted_message = f"[{algorithm_name}] {message} ({progress:.1f}%)" if message else f"[{algorithm_name}] Progress: {progress:.1f}%"
            
            # Adicionar log usando o novo método
            self.web_session_manager.add_log(
                self.web_session_id, 
                "INFO", 
                formatted_message,
                source=algorithm_name
            )
            
        except Exception as e:
            self.logger.error(f"Error sending algorithm log to web session: {e}")

    def _send_general_log(self, level: str, message: str) -> None:
        """Send general log to web session."""
        if not self.web_session_manager or not self.web_session_id:
            return
            
        try:
            # Adicionar log usando o novo método
            self.web_session_manager.add_log(
                self.web_session_id, 
                level, 
                message
            )
            
        except Exception as e:
            self.logger.error(f"Error sending general log to web session: {e}")

    def _update_current_execution(self, algorithm_name: str, progress: float, message: str, item_id: Optional[str] = None) -> None:
        """Update current execution details in the web session."""
        if not self.web_session_manager or not self.web_session_id:
            return
            
        try:
            # Extrair informações do item_id se disponível
            dataset = "-"
            run = "-"
            
            if item_id:
                # O item_id geralmente tem formato: dataset_algorithm_run ou similar
                parts = item_id.split("_")
                if len(parts) >= 3:
                    dataset = parts[0]
                    run = parts[-1]
                elif len(parts) >= 2:
                    dataset = parts[0]
                    run = parts[1]
            
            # Criar nome da execução atual
            execution_name = f"{algorithm_name}"
            if dataset != "-":
                execution_name += f" on {dataset}"
            
            current_execution = {
                "name": execution_name,
                "algorithm": algorithm_name,
                "dataset": dataset,
                "run": run,
                "progress": progress
            }
            
            # Atualizar sessão com dados de execução atual
            self.web_session_manager.update_session(self.web_session_id, {
                "current_execution": current_execution
            })
            
        except Exception as e:
            self.logger.error(f"Error updating current execution in web session: {e}")

    def report_progress(self, progress: float, message: str = "") -> None:
        """
        Report general progress update.

        Args:
            progress: Progress (0.0 to 100.0)
            message: Status message
        """
        if self.monitor:
            try:
                # For backward compatibility with direct progress reporting
                from src.presentation.monitoring.interfaces import ExecutionLevel
                self.monitor.update_hierarchy(ExecutionLevel.EXECUTION, "general", progress, message, None)
            except Exception as e:
                self.logger.error(f"Error reporting progress: {e}")

    def update_item(
        self, item_id: str, progress: float, message: str = "", context=None
    ) -> None:
        """
        Update an individual item.

        Args:
            item_id: Unique item ID
            progress: Progress (0.0 to 100.0)
            message: Status message
            context: Hierarchical context (optional)
        """
        if self.monitor:
            try:
                self.monitor.update_item(item_id, progress, message, context)
                
                # Update web session with algorithm progress
                if context and self.web_session_manager and self.web_session_id:
                    algorithm_name = getattr(context, 'algorithm_id', item_id)
                    if algorithm_name:
                        current_execution = {
                            "algorithm": algorithm_name,
                            "progress": progress
                        }
                        self.web_session_manager.update_session(self.web_session_id, {
                            "current_execution": current_execution
                        })
                        
            except Exception as e:
                self.logger.error(f"Error updating item {item_id}: {e}")

    def start_item(
        self,
        item_id: str,
        item_type: str = "repetition",
        context=None,
        metadata=None,
    ) -> None:
        """
        Start an individual item.

        Args:
            item_id: Unique item ID
            item_type: Item type (repetition, trial, sample, etc.)
            context: Hierarchical context (optional)
            metadata: Optional metadata
        """
        if self.monitor:
            try:
                self.monitor.start_item(item_id, item_type, context, metadata)
            except Exception as e:
                self.logger.error(f"Error starting item {item_id}: {e}")

    def update_hierarchy(
        self,
        level,
        level_id: str,
        progress: float,
        message: str = "",
        data=None,
    ) -> None:
        """
        Update hierarchical progress.

        Args:
            level: Hierarchical level (ExecutionLevel)
            level_id: Level ID
            progress: Progress (0.0 to 100.0)
            message: Status message
            data: Additional level-specific data
        """
        if self.monitor:
            try:
                self.monitor.update_hierarchy(level, level_id, progress, message, data)
                
                # Debug logging
                self.logger.info(f"Hierarchy update: {level} - {level_id} - {progress}% - data: {data}")
                
                # Update web session with progress
                if data and self.web_session_manager and self.web_session_id:
                    self.logger.info(f"Updating web session {self.web_session_id} with data: {data}")
                    
                    # Map internal data structure to hierarchical display structure
                    progress_data = self._map_to_hierarchical_structure(data)
                    overall_progress = self._calculate_overall_progress(data)
                    progress_data["overall_progress"] = overall_progress
                    
                    # Add current execution details
                    current_execution = {
                        "name": data.get("execution_name", "Execution"),
                        "dataset": data.get("dataset_name", "-"),
                        "algorithm": data.get("algorithm_config_name", "-"),
                        "run": str(progress_data.get("current_run", 1)),
                        "progress": overall_progress
                    }
                    
                    self.logger.info(f"Updating session with progress_data: {progress_data}")
                    
                    self.web_session_manager.update_session(self.web_session_id, {
                        "progress": progress_data,
                        "current_execution": current_execution
                    })
                    
            except Exception as e:
                self.logger.error(f"Erro ao atualizar hierarquia {level}: {e}")

    def finish_item(
        self,
        item_id: str,
        success: bool = True,
        result=None,
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
        if self.monitor:
            try:
                self.monitor.finish_item(item_id, success, result, error)
            except Exception as e:
                self.logger.error(f"Erro ao finalizar item {item_id}: {e}")

    def algorithm_callback(
        self,
        algorithm_name: str,
        progress: float,
        message: str = "",
        item_id: Optional[str] = None,
    ) -> None:
        """
        Callback direto do algoritmo durante execução.

        Args:
            algorithm_name: Nome do algoritmo
            progress: Progresso (0.0 a 100.0)
            message: Mensagem de status
            item_id: ID único do item (opcional)
        """
        if self.monitor:
            try:
                # Adaptar chamada para monitor hierárquico vs monitor legacy
                if hasattr(self.monitor, 'algorithm_callback'):
                    if "hierarchical" in type(self.monitor).__name__.lower():
                        # Monitor hierárquico usa parâmetros diferentes
                        self.monitor.algorithm_callback(
                            run_id=item_id,
                            algorithm_name=algorithm_name,
                            progress=progress,
                            message=message
                        )
                    else:
                        # Monitor legacy
                        self.monitor.algorithm_callback(
                            algorithm_name, progress, message, item_id
                        )
                
                # Enviar log do callback para a sessão web
                if self.web_session_manager and self.web_session_id:
                    self._send_algorithm_log(algorithm_name, progress, message)
                    
                    # Atualizar current_execution na sessão
                    self._update_current_execution(algorithm_name, progress, message, item_id)
                    
            except Exception as e:
                self.logger.error(
                    f"Erro no callback do algoritmo {algorithm_name}: {e}"
                )
