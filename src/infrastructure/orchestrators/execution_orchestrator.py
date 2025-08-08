"""
Execution Orchestrator

Centralizes all CSP algorithm execution logic,
for both single executions and batches.
"""

import json
import os
import pickle
import random
import re
import time
import types
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from pathlib import Path
from typing import Any, Dict, List, Optional

# IMPORTANT: Import algorithms first to load the global_registry
import algorithms
from src.domain import Dataset
from src.domain.algorithms import global_registry
from src.domain.errors import AlgorithmExecutionError
from src.infrastructure.logging_config import get_logger
from src.infrastructure.orchestrators.base_orchestrator import BaseOrchestrator
from src.infrastructure.resource_control import (
    ResourceController,
    create_resource_controller,
)


def _execute_algorithm_static(
    algorithm_name: str,
    dataset_data: dict,
    params: dict,
    rep_num: int,
    execution_metadata: dict,
):
    """
    Static function to execute algorithm - required for ProcessPoolExecutor.

    Args:
        algorithm_name: Name of algorithm to execute
        dataset_data: Serializable dataset data
        params: Algorithm parameters
        rep_num: Repetition number
        execution_metadata: Execution metadata

    Returns:
        tuple: (rep_num, repetition_results) or (rep_num, exception_info)
    """
    try:
        import time as time_module
        from src.domain.algorithms import global_registry
        from src.domain.dataset import Dataset

        # Get logger for this process
        logger = get_logger(__name__)

        # Reconstruct dataset
        dataset = Dataset(
            sequences=dataset_data["sequences"],
            metadata=dataset_data.get("metadata", {}),
        )

        # Set dataset name if provided (for identification purposes)
        dataset_name = dataset_data.get("name", "unknown")

        # Get algorithm instance
        algorithm_class = global_registry.get(algorithm_name)
        if not algorithm_class:
            raise AlgorithmExecutionError(f"Algorithm '{algorithm_name}' not found")

        # Create algorithm instance with correct parameters (STANDARD INTERFACE)
        # All algorithms use: __init__(strings, alphabet, **params) + run()
        algorithm_instance = algorithm_class(
            strings=dataset.sequences, alphabet=dataset.alphabet, **params
        )

        # Coletar warnings durante execução paralela
        collected_warnings = []
        collected_progress = []

        def warning_callback(message: str):
            """Callback para coletar warnings durante execução paralela."""
            warning_data = {
                "timestamp": time_module.time(),
                "message": message,
                "repetition": rep_num,
                "algorithm": algorithm_name,
            }
            collected_warnings.append(warning_data)
            logger.warning(f"Algorithm warning (rep {rep_num}): {message}")

        def progress_callback(message: str, progress: float = 0.0):
            """Callback para coletar progresso durante execução paralela."""
            progress_data = {
                "timestamp": time_module.time(),
                "message": message,
                "progress": progress,
                "repetition": rep_num,
                "algorithm": algorithm_name,
            }
            collected_progress.append(progress_data)
            logger.debug(f"Algorithm progress (rep {rep_num}): {message} - {progress}%")

        # Configurar callbacks
        algorithm_instance.set_warning_callback(warning_callback)
        algorithm_instance.set_progress_callback(progress_callback)

        # Execute algorithm using STANDARD INTERFACE (no parameters)
        logger.info(f"Starting repetition {rep_num} of {algorithm_name}")
        start_time = time_module.time()

        result = algorithm_instance.run()  # Standard interface - uses constructor data

        execution_time = time_module.time() - start_time

        # Build result data with only serializable components
        if isinstance(result, tuple) and len(result) >= 3:
            best_string, max_distance, metadata = result[:3]
        else:
            best_string, max_distance, metadata = str(result), 0, {}

        # Sanitize metadata
        clean_metadata = {}
        if isinstance(metadata, dict):
            for key, value in metadata.items():
                try:
                    # Test if value is serializable
                    json.dumps(value)
                    clean_metadata[key] = value
                except (TypeError, ValueError):
                    # Convert to string if not serializable
                    clean_metadata[key] = str(value)

        repetition_result = {
            "repetition": rep_num,
            "best_string": str(best_string),
            "distance": (
                int(max_distance) if isinstance(max_distance, (int, float)) else 0
            ),
            "execution_time": execution_time,
            "metadata": clean_metadata,
            "algorithm": algorithm_name,
            "dataset": dataset_name,
            # Adicionar dados coletados durante execução
            "warnings": collected_warnings,
            "progress_history": collected_progress,
        }

        logger.info(
            f"Completed repetition {rep_num} of {algorithm_name} in {execution_time:.2f}s"
        )
        return (rep_num, repetition_result)

    except Exception as e:
        logger = get_logger(__name__)
        logger.error(f"Error in repetition {rep_num} of {algorithm_name}: {e}")
        return (
            rep_num,
            {
                "error": str(e),
                "repetition": rep_num,
                "algorithm": algorithm_name,
                "dataset": dataset_data.get("name", "unknown"),
            },
        )


class ExecutionOrchestrator(BaseOrchestrator):
    """Orchestrator responsible for CSP algorithm execution."""

    def __init__(
        self,
        algorithm_registry,
        dataset_repository,
        monitoring_service=None,
        entrez_repository=None,
        session_manager=None,
    ):
        """
        Initialize execution orchestrator.

        Args:
            algorithm_registry: Algorithm registry
            dataset_repository: Dataset repository
            monitoring_service: Optional monitoring service
            entrez_repository: Optional Entrez dataset repository
            session_manager: Optional session manager for paths and sessions
        """
        super().__init__(monitoring_service)
        self._algorithm_registry = algorithm_registry
        self._dataset_repository = dataset_repository
        self._entrez_repository = entrez_repository
        self._session_manager = session_manager
        self._executions: Dict[str, Dict[str, Any]] = {}
        self._current_batch_config: Optional[Dict[str, Any]] = None
        self._partial_results_file: Optional[str] = None
        self._resource_controller: Optional[ResourceController] = None
        self._logger = get_logger(__name__)

    def execute(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Implements abstract method from BaseOrchestrator.

        Args:
            config: Execution configuration

        Returns:
            Dict[str, Any]: Resultado da execução
        """
        # Implementação padrão que delega para execute_batch
        results = self.execute_batch(config, self.monitoring_service)
        return {"results": results}

    def set_batch_config(self, batch_config: Dict[str, Any]) -> None:
        """Define current batch configuration."""
        self._current_batch_config = batch_config
        self._logger.debug(f"Batch configuration defined: {type(batch_config)}")

        # Initialize resource controller with batch resources config
        resources_config = batch_config.get("resources", {})
        if resources_config:
            try:
                self._resource_controller = create_resource_controller(resources_config)
                # Apply initial resource limits
                self._resource_controller.apply_cpu_limits()
                self._resource_controller.apply_memory_limits()
                self._logger.info(
                    "[RESOURCE] Resource controller initialized and limits applied"
                )
            except Exception as e:
                self._logger.warning(
                    f"[RESOURCE] Failed to initialize resource controller: {e}"
                )
                self._resource_controller = None

        # Configure partial saving if enabled
        if self._should_save_partial_results():
            self._setup_partial_results_file()

    def execute_single(
        self,
        algorithm_name: str,
        dataset: Dataset,
        params: Optional[Dict[str, Any]] = None,
        timeout: Optional[int] = None,
        monitoring_service=None,
        repetition_number: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Execute a single algorithm.

        Args:
            algorithm_name: Name of algorithm to execute
            dataset: Dataset for processing
            params: Algorithm-specific parameters
            timeout: Timeout in seconds
            monitoring_service: Optional monitoring service

        Returns:
            Dict[str, Any]: Execution result
        """
        from src.domain.algorithms import global_registry

        # Debug: Log available algorithms
        self._logger.debug(
            f"Available algorithms in registry: {list(global_registry.keys())}"
        )
        self._logger.debug(f"Looking for algorithm: {algorithm_name}")

        # Check if algorithm exists
        if algorithm_name not in global_registry:
            raise AlgorithmExecutionError(f"Algorithm '{algorithm_name}' not found")

        algorithm_class = global_registry[algorithm_name]
        params = params or {}

        # Apply infrastructure configurations if batch_config available
        if self._current_batch_config:
            infrastructure_config = self._current_batch_config.get("infrastructure", {})
            history_config = infrastructure_config.get("history", {})

            # Injetar parâmetros de histórico se habilitados
            # Support both template format (enabled/frequency) and implementation format (save_history/history_frequency)
            history_enabled = history_config.get(
                "enabled", history_config.get("save_history", False)
            )

            if history_enabled:
                params = params.copy()  # Não modificar o original
                params["save_history"] = True
                # Support both template format (frequency) and implementation format (history_frequency)
                frequency = history_config.get(
                    "frequency", history_config.get("history_frequency", 1)
                )
                params["history_frequency"] = frequency

        # Cria identificador único para execução
        execution_id = str(uuid.uuid4())

        try:
            # Registra início da execução
            start_time = time.time()
            self._executions[execution_id] = {
                "status": "running",
                "algorithm": algorithm_name,
                "start_time": start_time,
                "params": params.copy(),
            }

            # Instancia e executa algoritmo (STANDARD INTERFACE)
            # All algorithms use: __init__(strings, alphabet, **params) + run()
            algorithm = algorithm_class(
                strings=dataset.sequences, alphabet=dataset.alphabet, **params
            )

            # Configurar callback de progresso se fornecido
            if monitoring_service:
                start_time = time.time()

                def progress_callback(message: str):
                    # Tentar extrair progresso da mensagem
                    progress = 0.0

                    # Padrão 1: "Geração X/Y" (BLF-GA)
                    gen_match = re.search(r"Geração (\d+)/(\d+)", message)
                    if gen_match:
                        current_gen = int(gen_match.group(1))
                        max_gen = int(gen_match.group(2))
                        progress = (current_gen / max_gen) * 100.0

                    # Padrão 2: "iteração X" (H3-CSP)
                    elif "iteração" in message:
                        iter_match = re.search(r"iteração (\d+)", message)
                        if iter_match:
                            # Estimar progresso baseado no tempo (máximo 100% em 60s)
                            elapsed = time.time() - start_time
                            progress = min(elapsed / 60.0 * 100.0, 95.0)

                    # Padrão 3: Mensagens de fase
                    elif "Analisando blocos" in message:
                        progress = 10.0
                    elif "Fusão de blocos" in message:
                        progress = 30.0
                    elif "Refinamento" in message:
                        progress = 60.0
                    elif "ótima encontrada" in message or "Convergiu" in message:
                        progress = 100.0
                    elif "Timeout" in message:
                        progress = 100.0

                    # Usar progresso baseado no tempo se não conseguir extrair
                    if progress == 0.0:
                        elapsed = time.time() - start_time
                        progress = min(
                            elapsed / 30.0 * 100.0, 90.0
                        )  # Máximo 90% por tempo

                    # Garantir que progresso está no range 0-100
                    progress = max(0.0, min(100.0, progress))

                    # Usando report_algorithm_progress
                    monitoring_service.report_algorithm_progress(
                        algorithm_name=algorithm_name,
                        progress_percent=progress,
                        message=message,
                        item_id=execution_id,
                    )

                algorithm.set_progress_callback(progress_callback)

                # Configurar callback de warning se fornecido
                def warning_callback(message: str):
                    # Reportar warning com contexto de execução única
                    monitoring_service.report_warning(
                        algorithm_name=algorithm_name,
                        warning_message=message,
                        item_id=execution_id,
                        execution_context={
                            "execution_type": "single",
                            "dataset_size": len(dataset.sequences),
                            "dataset_length": (
                                len(dataset.sequences[0]) if dataset.sequences else 0
                            ),
                            "algorithm_params": params,
                        },
                    )

                algorithm.set_warning_callback(warning_callback)

            # Apply resource controls and execute algorithm
            try:
                if self._resource_controller:
                    # Check batch timeout before starting algorithm
                    self._resource_controller.check_batch_timeout()

                    # Execute with algorithm timeout
                    with self._resource_controller.algorithm_timeout():
                        best_string, max_distance, metadata = algorithm.run()
                else:
                    # Fallback execution without resource control
                    best_string, max_distance, metadata = algorithm.run()
            except TimeoutError as e:
                self._logger.error(
                    f"[RESOURCE] Algorithm {algorithm_name} timed out: {e}"
                )
                # Return timeout result
                end_time = time.time()
                result = {
                    "algorithm": algorithm_name,
                    "best_string": "",
                    "max_distance": -1,
                    "execution_time": end_time - start_time,
                    "execution_id": execution_id,
                    "params": params,
                    "metadata": {"error": str(e), "timeout": True},
                    "dataset": {
                        "size": len(dataset.sequences),
                        "length": len(dataset.sequences[0]) if dataset.sequences else 0,
                        "alphabet": dataset.alphabet,
                    },
                    "status": "timeout",
                }

                # Emit algorithm finished event for timeout
                if monitoring_service:
                    monitoring_service.report_algorithm_finished(
                        algorithm_name=algorithm_name,
                        success=False,
                        best_string="",
                        max_distance=-1,
                        execution_time=end_time - start_time,
                        metadata={"error": str(e), "timeout": True},
                        repetition_number=repetition_number,
                    )

                return result
            except Exception as e:
                if "Resource limit" in str(e) or "Memory limit" in str(e):
                    self._logger.error(
                        f"[RESOURCE] Algorithm {algorithm_name} exceeded resource limits: {e}"
                    )
                    # Return resource limit result
                    end_time = time.time()
                    result = {
                        "algorithm": algorithm_name,
                        "best_string": "",
                        "max_distance": -1,
                        "execution_time": end_time - start_time,
                        "execution_id": execution_id,
                        "params": params,
                        "metadata": {"error": str(e), "resource_limit_exceeded": True},
                        "dataset": {
                            "size": len(dataset.sequences),
                            "length": (
                                len(dataset.sequences[0]) if dataset.sequences else 0
                            ),
                            "alphabet": dataset.alphabet,
                        },
                        "status": "resource_limit",
                    }

                    # Emit algorithm finished event for resource limit
                    if monitoring_service:
                        monitoring_service.report_algorithm_finished(
                            algorithm_name=algorithm_name,
                            success=False,
                            best_string="",
                            max_distance=-1,
                            execution_time=end_time - start_time,
                            metadata={"error": str(e), "resource_limit_exceeded": True},
                            repetition_number=repetition_number,
                        )

                    return result
                else:
                    # Re-raise other exceptions
                    raise

            end_time = time.time()

            # Constroi resultado
            result = {
                "algorithm": algorithm_name,
                "best_string": best_string,
                "max_distance": max_distance,
                "execution_time": end_time - start_time,
                "execution_id": execution_id,
                "params": params,
                "metadata": metadata,
                "dataset": {
                    "size": len(dataset.sequences),
                    "length": len(dataset.sequences[0]) if dataset.sequences else 0,
                    "alphabet": dataset.alphabet,
                },
                "status": "completed",
            }

            # Atualiza status
            self._executions[execution_id].update(
                {"status": "completed", "result": result, "end_time": end_time}
            )

            # Emit algorithm finished event if monitoring service available
            if monitoring_service:
                monitoring_service.report_algorithm_finished(
                    algorithm_name=algorithm_name,
                    success=True,
                    best_string=best_string,
                    max_distance=max_distance,
                    execution_time=end_time - start_time,
                    metadata=metadata,
                    repetition_number=repetition_number,
                )

            return result

        except Exception as e:
            # Registra erro
            end_time = time.time()
            self._executions[execution_id].update(
                {"status": "failed", "error": str(e), "end_time": end_time}
            )

            # Emit algorithm finished event for error
            if monitoring_service:
                monitoring_service.report_algorithm_finished(
                    algorithm_name=algorithm_name,
                    success=False,
                    best_string="",
                    max_distance=-1,
                    execution_time=end_time - start_time,
                    metadata={"error": str(e)},
                    repetition_number=repetition_number,
                )

            raise AlgorithmExecutionError(
                f"Erro na execução de '{algorithm_name}': {e}"
            )

    def execute_batch(
        self, batch_config: Dict[str, Any], monitoring_service=None
    ) -> List[Dict[str, Any]]:
        """
        Executa batch de algoritmos.

        Args:
            batch_config: Configuração do batch
            monitoring_service: Serviço de monitoramento opcional

        Returns:
            List[Dict[str, Any]]: Lista de resultados da execução
        """
        self.set_batch_config(batch_config)

        # Detectar tipo de batch
        task_type_str = batch_config.get("task", {}).get("type", "execution")

        # Monitoramento será iniciado pela camada de aplicação para evitar duplicação
        
        results = []

        self._logger.debug(f"Task type detectado: {task_type_str}")

        if task_type_str == "execution" and "task" in batch_config:
            results = self._execute_structured_batch(batch_config, monitoring_service)
        elif "experiments" in batch_config:
            results = self._execute_legacy_batch(batch_config, monitoring_service)
        else:
            # Estrutura simples - compatibilidade
            algorithms = batch_config.get("algorithms", [])
            datasets = batch_config.get("datasets", [])
            default_params = batch_config.get("params", {})

            for algorithm_name in algorithms:
                for dataset in datasets:
                    try:
                        result = self.execute_single(
                            algorithm_name=algorithm_name,
                            dataset=dataset,
                            params=default_params.get(algorithm_name, {}),
                        )
                        results.append(result)

                        # Salvar resultado parcial se habilitado
                        if self._should_save_partial_results():
                            self._save_partial_result(result)

                    except Exception as e:
                        error_result = {
                            "algorithm": algorithm_name,
                            "dataset": getattr(dataset, "name", "unknown"),
                            "status": "failed",
                            "error": str(e),
                        }
                        results.append(error_result)

                        # Salvar resultado de erro parcial se habilitado
                        if self._should_save_partial_results():
                            self._save_partial_result(error_result)

        return results

    def _execute_structured_batch(
        self, batch_config: Dict[str, Any], monitoring_service=None
    ) -> List[Dict[str, Any]]:
        """Execute batch with new structure (datasets + algorithms + executions)."""
        self._logger.debug("Processando batch de execução estruturado")

        # Nova estrutura com datasets e algoritmos por ID
        datasets_config = batch_config.get("datasets", [])
        algorithm_configs = batch_config.get("algorithm_configs", [])
        tasks = batch_config["execution"]["tasks"]

        self._logger.debug(f"Datasets: {len(datasets_config)}")
        self._logger.debug(f"Algorithm configs: {len(algorithm_configs)}")
        self._logger.debug(f"Tasks: {len(tasks)}")

        # Inicializar dados de monitoramento
        if monitoring_service:
            self._setup_monitoring_data(
                tasks, datasets_config, algorithm_configs, monitoring_service
            )

        from src.infrastructure import FileDatasetRepository

        dataset_repo = FileDatasetRepository("./datasets")

        results = []
        task_index = 0

        for task in tasks:
            task_index += 1
            task_id = task.get("id", f"task_{task_index}")

            # Notificar início da execução
            if monitoring_service:
                monitoring_service.notify_execution_started(
                    execution_name=task_id,
                    metadata={
                        "index": task_index,
                        "total_tasks": len(tasks),
                        "datasets": task.get("datasets", []),
                        "algorithm_configs": task.get("algorithm_configs", []),
                    },
                )

            # Iterar sobre datasets primeiro
            dataset_ids = task["datasets"]
            algorithm_config_ids = task["algorithm_configs"]

            for dataset_idx, dataset_id in enumerate(dataset_ids, 1):

                dataset_config = next(
                    (d for d in datasets_config if d["id"] == dataset_id), None
                )
                if not dataset_config:
                    results.append(
                        {
                            "task_id": task.get("id", "unknown"),
                            "dataset_id": dataset_id,
                            "status": "error",
                            "error": f"Dataset com ID '{dataset_id}' não encontrado",
                        }
                    )
                    continue

                # Iterar sobre configurações de algoritmo para este dataset
                for algo_config_idx, algorithm_config_id in enumerate(algorithm_config_ids, 1):

                    # Resolver configuração do algoritmo
                    algorithm_config = next(
                        (a for a in algorithm_configs if a["id"] == algorithm_config_id), None
                    )
                    if not algorithm_config:
                        self._logger.error(
                            f"Algoritmo com ID '{algorithm_config_id}' não encontrado"
                        )
                        continue

                    # Atualizar informações do dataset no monitoramento
                    if monitoring_service:
                        # Contar algoritmos únicos desta configuração
                        unique_algorithms = set(algorithm_config["algorithms"])

                        # Atualizar hierarquia de dataset (que já inclui execução)
                        from src.application.monitoring.progress_events import (
                            ExecutionLevel,
                        )

                        # Obter nome do dataset
                        dataset_name = dataset_config.get("name", dataset_id)

                        # Obter nome da configuração de algoritmo
                        algorithm_config_name = algorithm_config.get(
                            "name", "Algorithms"
                        )

                        # Use monitoring service method that handles None monitor
                        if monitoring_service:
                            monitoring_service.update_hierarchy(
                                level=ExecutionLevel.DATASET,
                                level_id=f"{dataset_id}_{algorithm_config_id}",
                                progress=0.0,
                                message=f"Processando dataset {dataset_name}",
                                data={
                                    "task_id": task_id,
                                    "config_index": task_index,
                                    "total_configs": len(tasks),
                                    "dataset_name": dataset_name,
                                    "dataset_index": dataset_idx,
                                    "total_datasets": len(dataset_ids),
                                    "algorithm_config_name": algorithm_config_name,
                                    "algorithm_config_index": algo_config_idx,
                                    "total_algorithm_configs": len(algorithm_config_ids),
                                    "total_algorithms": len(unique_algorithms),
                                },
                            )

                    # Carregar dataset e executar algoritmos desta configuração
                    dataset_results = self._execute_dataset_algorithms_for_config(
                        task,
                        dataset_config,
                        dataset_id,
                        dataset_repo,
                        algorithm_config,
                        monitoring_service,
                    )
                    results.extend(dataset_results)

            # Configurações completadas são controladas pela hierarquia
            # Não precisamos mais usar update_execution_data

        return results

    def _execute_legacy_batch(
        self, batch_config: Dict[str, Any], monitoring_service=None
    ) -> List[Dict[str, Any]]:
        """Execute batch with legacy structure (experiments)."""
        from src.infrastructure import FileDatasetRepository

        dataset_repo = FileDatasetRepository("./datasets")

        results = []
        for exp in batch_config.get("experiments", []):
            try:
                print(f"[DEBUG] Processando experimento: {exp}")
                # Carrega dataset
                dataset = dataset_repo.load(exp["dataset"])
                print(f"[DEBUG] Dataset loaded: {len(dataset.sequences)} sequences")

                result = self.execute_single(
                    algorithm_name=exp["algorithm"],
                    dataset=dataset,
                    params=exp.get("params", {}),
                )
                result["dataset"] = exp["dataset"]
                results.append(result)
                print("[DEBUG] Resultado adicionado com sucesso")

                # Salvar resultado parcial se habilitado
                if self._should_save_partial_results():
                    self._save_partial_result(result)

            except Exception as e:
                print(f"[DEBUG] Erro no experimento: {e}")
                error_result = {
                    "algorithm": exp["algorithm"],
                    "dataset": exp["dataset"],
                    "status": "error",
                    "error": str(e),
                }
                results.append(error_result)

                # Salvar resultado de erro parcial se habilitado
                if self._should_save_partial_results():
                    self._save_partial_result(error_result)

        return results

    def _should_save_partial_results(self) -> bool:
        """Check if partial results should be saved."""
        if not self._current_batch_config:
            return False

        infrastructure = self._current_batch_config.get("infrastructure", {})
        result_config = infrastructure.get("result", {})
        return result_config.get("save_partial_results", False)

    def _setup_partial_results_file(self) -> None:
        """Configure file for partial results saving."""
        try:
            # Use the session manager passed in constructor if available
            if self._session_manager:
                session_folder = self._session_manager.create_session()
                results_dir = Path(self._session_manager.get_result_dir())
                print(f"📁 Sessão criada: {session_folder}")
                print(f"📁 Salvando resultados parciais em: {results_dir}")
            else:
                # Create a local SessionManager as fallback
                from src.infrastructure import SessionManager

                session_manager = SessionManager(self._current_batch_config or {})
                session_folder = session_manager.create_session()
                results_dir = Path(session_manager.get_result_dir())
                print(f"📁 Sessão criada: {session_folder}")
                print(f"📁 Salvando resultados parciais em: {results_dir}")
        except Exception as e:
            # Fallback para diretório padrão
            import time

            timestamp = time.strftime("%Y%m%d_%H%M%S")
            results_dir = Path("outputs") / timestamp
            print(f"📁 Usando diretório fallback: {results_dir}, erro: {e}")

        results_dir.mkdir(parents=True, exist_ok=True)
        self._partial_results_file = str(results_dir / "partial_results.json")

        print(f"💾 Arquivo de resultados parciais: {self._partial_results_file}")

        # Inicializar arquivo com array vazio
        with open(self._partial_results_file, "w", encoding="utf-8") as f:
            json.dump([], f)

        print("✅ Sistema de salvamento parcial inicializado")

    def _save_partial_result(self, result: Dict[str, Any]) -> None:
        """Save a partial result to file."""
        if not self._partial_results_file:
            return

        try:
            # Carregar resultados existentes
            if os.path.exists(self._partial_results_file):
                with open(self._partial_results_file, encoding="utf-8") as f:
                    existing_results = json.load(f)
            else:
                existing_results = []

            # Adicionar novo resultado
            existing_results.append(result)

            # Salvar de volta
            with open(self._partial_results_file, "w", encoding="utf-8") as f:
                json.dump(existing_results, f, indent=2, ensure_ascii=False)

            print(
                f"💾 Resultado salvo [{len(existing_results)}]: {result.get('algorithm', 'N/A')} - {result.get('status', 'N/A')}"
            )

        except Exception as e:
            print(f"⚠️  Erro ao salvar resultado parcial: {e}")

    # Métodos auxiliares para organização do código...
    def _setup_monitoring_data(
        self, tasks, datasets_config, algorithm_configs, monitoring_service
    ):
        """Configura dados iniciais de monitoramento."""
        # Calcular totais de datasets considerando todas as tasks
        total_dataset_executions = 0
        for task in tasks:
            total_dataset_executions += len(task.get("datasets", []))

        # Contar total de algoritmos únicos
        unique_algorithms = set()
        for task in tasks:
            algorithm_config_ids = task["algorithm_configs"]
            for algorithm_config_id in algorithm_config_ids:
                algorithm_config = next(
                    (a for a in algorithm_configs if a["id"] == algorithm_config_id), None
                )
                if algorithm_config:
                    algorithms_list = algorithm_config["algorithms"]
                    unique_algorithms.update(algorithms_list)

        total_algorithms = len(unique_algorithms)

        # Dados iniciais configurados, mas não usamos mais update_execution_data
        # O monitoramento agora é feito através de update_hierarchy

        # Forçar primeira atualização para mostrar interface
        time.sleep(0.1)

    def _update_execution_monitoring(
        self,
        execution,
        execution_index,
        task_id,
        datasets_config,
        algorithms_config,
        monitoring_service,
    ):
        """Atualiza dados de monitoramento para execução atual."""
        # Contar algoritmos únicos nesta execução
        execution_algorithms = set()
        for algo_id in execution["algorithms"]:
            algo_config = next(
                (a for a in algorithms_config if a["id"] == algo_id), None
            )
            if algo_config:
                execution_algorithms.update(algo_config["algorithms"])

        # Atualizar hierarquia de execução
        from src.application.monitoring.progress_events import ExecutionLevel

        # Obter total de execuções do contexto (será passado pelo orchestrator)
        # Por enquanto, usando o index como fallback
        total_executions = execution_index  # Isso será melhorado no orchestrator

        # Use monitoring service method that handles None monitor
        if monitoring_service:
            monitoring_service.update_hierarchy(
                level=ExecutionLevel.EXECUTION,
                level_id=task_id,
                progress=0.0,
                message=f"Iniciando execução {task_id}",
                data={
                    "task_id": task_id,
                    "config_index": execution_index,
                    "total_configs": total_executions,
                },
            )

    def _execute_dataset_algorithms_for_config(
        self,
        task,
        dataset_config,
        dataset_id,
        dataset_repo,
        algorithm_config,
        monitoring_service,
    ):
        """Executa algoritmos de uma configuração específica para um dataset."""
        results = []

        try:
            # Validar algorithm_config
            if algorithm_config is None:
                error_msg = f"Algorithm config is None for dataset {dataset_id}"
                self._logger.error(error_msg)
                return [
                    {
                        "task_id": task.get("id", "unknown"),
                        "dataset_id": dataset_id,
                        "status": "error",
                        "error": error_msg,
                        "execution_time": 0.0,
                    }
                ]

            # Carregar dataset baseado no tipo
            dataset_type = dataset_config["type"]

            if dataset_type == "file":
                filename = dataset_config["parameters"]["filename"]
                dataset = dataset_repo.load(filename)
            elif dataset_type == "synthetic":
                # For synthetic datasets, create using generator
                dataset = self._create_dataset_from_config(dataset_config)
            elif dataset_type == "entrez":
                # For entrez datasets, use the entrez repository
                if not self._entrez_repository:
                    raise ValueError(
                        "Entrez dataset repository not available. "
                        "Check NCBI_EMAIL environment variable and Biopython installation."
                    )

                params = dataset_config.get("parameters", {})
                query = params.get("query")
                if not query:
                    raise ValueError(
                        "Field 'query' is required for entrez dataset type"
                    )

                db = params.get("db", "nucleotide")
                retmax = params.get("retmax", 20)

                # Fetch dataset from NCBI
                # Extract only additional parameters, avoiding duplicates
                additional_params = {
                    k: v
                    for k, v in params.items()
                    if k not in ["query", "db", "retmax"]
                }

                sequences, metadata = self._entrez_repository.fetch_dataset(
                    query=query, db=db, retmax=retmax, **additional_params
                )

                # Create Dataset object
                from src.domain import Dataset

                dataset = Dataset(
                    sequences=sequences,
                    metadata={
                        "type": "entrez",
                        "query": query,
                        "db": db,
                        "retmax": retmax,
                        "n_obtained": len(sequences),
                        "L": len(sequences[0]) if sequences else 0,
                        **metadata,
                    },
                )

                self._logger.info(
                    "Created Entrez dataset: n=%d, L=%d, query='%s'",
                    len(sequences),
                    len(sequences[0]) if sequences else 0,
                    query,
                )
            else:
                raise ValueError(f"Unsupported dataset type: {dataset_type}")

            self._logger.info(
                f"Dataset {dataset_id} loaded: {len(dataset.sequences)} sequences"
            )

            # Execute algorithms for this configuration
            algorithm_names = algorithm_config["algorithms"]
            algorithm_params = algorithm_config.get("algorithm_params", {})
            repetitions = task.get("repetitions", 1)

            for algorithm_name in algorithm_names:
                # Obter parâmetros específicos do algoritmo
                params = algorithm_params.get(algorithm_name, {})

                # Executar repetições com paralelismo
                algorithm_results = self._execute_algorithm_repetitions_parallel(
                    algorithm_name=algorithm_name,
                    dataset=dataset,
                    params=params,
                    repetitions=repetitions,
                    execution_context={
                        "task_id": task.get("id", "unknown"),
                        "dataset_id": dataset_id,
                        "algorithm_config_id": algorithm_config["id"],
                        "config_id": algorithm_config["id"],
                    },
                    monitoring_service=monitoring_service,
                )

                results.extend(algorithm_results)

        except Exception as e:
            self._logger.error(
                f"Erro no carregamento/processamento do dataset {dataset_id}: {e}"
            )
            results.append(
                {
                    "task_id": task.get("id", "unknown"),
                    "dataset_id": dataset_id,
                    "status": "error",
                    "error": str(e),
                    "execution_time": 0.0,
                }
            )

        return results

    def _execute_single_repetition(
        self,
        algorithm_name: str,
        dataset_data: Dict[str, Any],  # Dados serializáveis do dataset
        params: Dict[str, Any],
        execution_context: Dict[str, Any],
        rep_number: int,
        total_repetitions: int,
    ) -> Dict[str, Any]:
        """
        Executa uma única repetição de um algoritmo.

        Este método é projetado para ser usado com ProcessPoolExecutor,
        portanto deve ser independente de estado do orchestrator.

        Args:
            algorithm_name: Nome do algoritmo a executar
            dataset_data: Dados serializáveis do dataset (sequences, alphabet)
            params: Parâmetros do algoritmo
            execution_context: Contexto da execução (nomes, IDs, etc.)
            rep_number: Número da repetição (1-based)
            total_repetitions: Total de repetições

        Returns:
            Dict[str, Any]: Resultado da execução com contexto
        """
        try:
            # Recriar objeto Dataset a partir dos dados serializáveis
            from src.domain import Dataset

            dataset = Dataset(dataset_data["sequences"], dataset_data["alphabet"])

            # Executar algoritmo (sem monitoring_service pois não é thread-safe)
            result = self.execute_single(
                algorithm_name, dataset, params, monitoring_service=None
            )

            # Sanitizar metadados para serialização
            if "metadata" in result:
                result["metadata"] = self._sanitize_metadata_for_process(
                    result["metadata"]
                )

            # Adicionar informações de contexto
            result.update(
                {
                    "task_id": execution_context.get("task_id", "unknown"),
                    "dataset_id": execution_context.get("dataset_id", "unknown"),
                    "algorithm_id": execution_context.get("algorithm_id", "unknown"),
                    "algorithm_name": algorithm_name,
                    "repetition": rep_number,
                    "total_repetitions": total_repetitions,
                    "status": "success",
                }
            )

            # Sanitizar todo o resultado para garantir serialização
            result = self._sanitize_metadata_for_process(result)

            # Teste final de serialização para debug
            try:
                import pickle

                pickle.dumps(result)
            except Exception as pickle_error:
                self._logger.error(
                    f"Resultado ainda não é serializável após sanitização: {pickle_error}"
                )
                # Fallback: converter tudo para strings
                result = {
                    "algorithm": algorithm_name,
                    "best_string": str(result.get("best_string", "")),
                    "max_distance": result.get("max_distance", -1),
                    "execution_time": result.get("execution_time", 0.0),
                    "execution_id": str(result.get("execution_id", "")),
                    "params": str(result.get("params", {})),
                    "metadata": str(result.get("metadata", {})),
                    "dataset": {
                        "size": len(dataset_data["sequences"]),
                        "length": (
                            len(dataset_data["sequences"][0])
                            if dataset_data["sequences"]
                            else 0
                        ),
                        "alphabet": dataset_data["alphabet"],
                    },
                    "status": "completed",
                    "task_id": execution_context.get("task_id", "unknown"),
                    "dataset_id": execution_context.get("dataset_id", "unknown"),
                    "algorithm_id": execution_context.get("algorithm_id", "unknown"),
                    "algorithm_name": algorithm_name,
                    "repetition": rep_number,
                    "total_repetitions": total_repetitions,
                    "serialization_fallback": True,
                }

            return result

        except Exception as e:
            error_result = {
                "task_id": execution_context.get("task_id", "unknown"),
                "dataset_id": execution_context.get("dataset_id", "unknown"),
                "algorithm_id": execution_context.get("algorithm_id", "unknown"),
                "algorithm_name": algorithm_name,
                "repetition": rep_number,
                "total_repetitions": total_repetitions,
                "status": "error",
                "error": str(e),
                "execution_time": 0.0,
            }

            return error_result

    def _sanitize_metadata_for_process(self, metadata: Any) -> Any:
        """
        Sanitiza metadados para serialização em ProcessPoolExecutor.

        Remove ou converte objetos não serializáveis para garantir que
        os resultados possam ser transferidos entre processos.
        """
        import pickle
        import types

        if isinstance(metadata, dict):
            sanitized = {}
            for key, value in metadata.items():
                try:
                    # Testa se o valor pode ser serializado
                    pickle.dumps(value)
                    sanitized[key] = self._sanitize_metadata_for_process(value)
                except (TypeError, AttributeError):
                    # Converte objetos não serializáveis
                    if isinstance(value, types.ModuleType):
                        sanitized[key] = f"<module '{value.__name__}'>"
                    elif callable(value):
                        sanitized[key] = (
                            f"<function '{getattr(value, '__name__', str(value))}'>"
                        )
                    elif hasattr(value, "__dict__") and hasattr(value, "__class__"):
                        # Para objetos complexos, extrai atributos serializáveis
                        try:
                            class_name = value.__class__.__name__
                            repr_str = str(value)
                            sanitized[key] = {
                                "__class__": class_name,
                                "__repr__": repr_str,
                            }
                        except:
                            sanitized[key] = str(value)
                    else:
                        sanitized[key] = str(value)
            return sanitized
        elif isinstance(metadata, (list, tuple)):
            sanitized_items = []
            for item in metadata:
                try:
                    pickle.dumps(item)
                    sanitized_items.append(self._sanitize_metadata_for_process(item))
                except (TypeError, AttributeError):
                    sanitized_items.append(str(item))
            return type(metadata)(sanitized_items)
        else:
            # Para tipos primitivos, retorna como está
            return metadata

    def _get_max_workers(self) -> int:
        """
        Obtém o número máximo de workers para paralelização.

        Returns:
            int: Número de workers a usar
        """
        if self._current_batch_config and self._current_batch_config is not None:
            resources = self._current_batch_config.get("resources", {})
            parallel_config = resources.get("parallel", {})

            # Handle different formats for parallel config
            if isinstance(parallel_config, int):
                # Direct integer value for max_workers
                return parallel_config if parallel_config > 0 else 1
            elif isinstance(parallel_config, dict) and parallel_config is not None:
                # Dictionary format with max_workers key
                max_workers = parallel_config.get("max_workers")
                if max_workers is not None and max_workers > 0:
                    return max_workers

        # Fallback para número de CPUs
        return cpu_count() or 1

    def _execute_dataset_algorithms(
        self,
        execution,
        dataset_config,
        dataset_id,
        dataset_repo,
        algorithms_config,
        monitoring_service,
    ):
        """Executa algoritmos para um dataset específico."""
        results = []

        try:
            # Load dataset
            if dataset_config["type"] == "file":
                filename = dataset_config["parameters"]["filename"]
                dataset = dataset_repo.load(filename)
            else:
                # For synthetic datasets, create using generator
                dataset = self._create_dataset_from_config(dataset_config)

            self._logger.info(
                f"Dataset {dataset_id} loaded: {len(dataset.sequences)} sequences"
            )

            # Get algorithm configurations from execution
            algorithm_ids = execution["algorithms"]
            repetitions = execution.get("repetitions", 1)

            for algorithm_id in algorithm_ids:
                # Resolver configuração do algoritmo
                algorithm_config = next(
                    (a for a in algorithms_config if a["id"] == algorithm_id), None
                )
                if not algorithm_config:
                    self._logger.error(
                        f"Algoritmo com ID '{algorithm_id}' não encontrado"
                    )
                    results.append(
                        {
                            "task_id": execution.get("id", "unknown"),
                            "dataset_id": dataset_id,
                            "algorithm_id": algorithm_id,
                            "status": "error",
                            "error": f"Algoritmo com ID '{algorithm_id}' não encontrado",
                        }
                    )
                    continue

                # Executar cada algoritmo da configuração
                algorithm_names = algorithm_config["algorithms"]
                algorithm_params = algorithm_config.get("algorithm_params", {})

                for algorithm_name in algorithm_names:
                    # Obter parâmetros específicos do algoritmo
                    params = algorithm_params.get(algorithm_name, {})

                    # Executar repetições com paralelismo
                    algorithm_results = self._execute_algorithm_repetitions_parallel(
                        algorithm_name=algorithm_name,
                        dataset=dataset,
                        params=params,
                        repetitions=repetitions,
                        execution_context={
                            "task_id": execution.get("id", "unknown"),
                            "dataset_id": dataset_id,
                            "algorithm_id": algorithm_id,
                            "config_id": algorithm_id,
                        },
                        monitoring_service=monitoring_service,
                    )

                    results.extend(algorithm_results)

        except Exception as e:
            self._logger.error(
                f"Erro no carregamento/processamento do dataset {dataset_id}: {e}"
            )
            results.append(
                {
                    "task_id": execution.get("id", "unknown"),
                    "dataset_id": dataset_id,
                    "status": "error",
                    "error": f"Erro no dataset: {e}",
                }
            )

        return results

    def _create_dataset_from_config(self, dataset_config: Dict[str, Any]):
        """Cria dataset a partir da configuração."""
        from src.domain.dataset import SyntheticDatasetGenerator

        dataset_type = dataset_config["type"]
        params = dataset_config.get("parameters", {})

        if dataset_type == "synthetic":
            generator = SyntheticDatasetGenerator()

            # Se há parâmetros de noise, usar generate_from_center
            if "noise" in params and params["noise"] > 0:
                # Gerar string central primeiro
                n = params.get("n", 10)
                L = params.get("L", 20)
                alphabet = params.get("alphabet", "ACTG")
                seed = params.get("seed")

                # Criar string central aleatória
                import random

                rng = random.Random(seed)
                center = "".join(rng.choice(alphabet) for _ in range(L))

                return generator.generate_from_center(
                    center=center,
                    n=n,
                    noise_rate=params.get("noise", 0.0),
                    alphabet=alphabet,
                    seed=seed,
                )
            else:
                # Usar generate_random para datasets sem ruído
                return generator.generate_random(
                    n=params.get("n", 10),
                    length=params.get("L", 20),
                    alphabet=params.get("alphabet", "ACTG"),
                    seed=params.get("seed"),
                )
        else:
            raise ValueError(
                f"Dataset type '{dataset_type}' not supported in synthetic dataset generation"
            )

    def get_execution_status(self, execution_id: str) -> str:
        """Obtém status de uma execução específica."""
        if execution_id not in self._executions:
            return "not_found"
        return self._executions[execution_id]["status"]

    def cancel_execution(self, execution_id: str) -> bool:
        """Cancela uma execução em andamento."""
        if execution_id not in self._executions:
            return False

        execution = self._executions[execution_id]
        if execution["status"] == "running":
            execution["status"] = "cancelled"
            execution["end_time"] = time.time()
            return True

        return False

    def _execute_algorithm_repetitions_parallel(
        self,
        algorithm_name: str,
        dataset: Dataset,
        params: Dict[str, Any],
        repetitions: int,
        execution_context: Dict[str, Any],
        monitoring_service=None,
    ) -> List[Dict[str, Any]]:
        """
        Executa repetições de um algoritmo em paralelo usando ProcessPoolExecutor.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para processamento
            params: Parâmetros do algoritmo
            repetitions: Número de repetições
            execution_context: Contexto da execução
            monitoring_service: Serviço de monitoramento

        Returns:
            List[Dict[str, Any]]: Lista de resultados das repetições
        """
        max_workers = self._get_max_workers()

        # Se max_workers = 1, usar execução sequencial
        if max_workers == 1:
            return self._execute_algorithm_repetitions_sequential(
                algorithm_name,
                dataset,
                params,
                repetitions,
                execution_context,
                monitoring_service,
            )

        # Execução paralela
        self._logger.debug(
            f"Executando {repetitions} repetições de {algorithm_name} com {max_workers} workers"
        )

        results = []

        # Preparar argumentos para ProcessPoolExecutor
        args_list = []
        for rep in range(repetitions):
            # Preparar dados serializáveis do dataset
            dataset_data = {
                "name": execution_context.get(
                    "dataset_id", "unknown"
                ),  # Use dataset_id from context
                "sequences": dataset.sequences,
                "metadata": getattr(dataset, "metadata", {}),  # Safe access to metadata
            }

            # Preparar metadados de execução serializáveis
            execution_metadata = {
                "session_id": str(execution_context.get("session_id", "")),
                "timestamp": execution_context.get("timestamp", ""),
                "repetitions": repetitions,
            }

            args_list.append(
                (
                    algorithm_name,
                    dataset_data,
                    params,
                    rep + 1,  # rep_num
                    execution_metadata,
                )
            )

        # Executar em paralelo
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submeter todas as tarefas usando função estática
            future_to_rep = {}
            for i, args in enumerate(args_list):
                future = executor.submit(_execute_algorithm_static, *args)
                future_to_rep[future] = i + 1

            # Coletar resultados conforme completam
            for future in as_completed(future_to_rep):
                rep_number = future_to_rep[future]
                rep_id = f"{algorithm_name}_{execution_context.get('dataset_id', 'unknown')}_{rep_number}"

                try:
                    # Inicializar monitoramento se disponível
                    if monitoring_service:
                        # Iniciar item antes da execução
                        monitoring_service.start_item(rep_id, "repetition")

                    # Obter resultado da função estática
                    returned_rep_num, result_data = future.result()

                    # Processar warnings coletados durante execução paralela
                    if monitoring_service and isinstance(result_data, dict):
                        # Emitir warnings individuais
                        warnings = result_data.get("warnings", [])
                        for warning_data in warnings:
                            monitoring_service.report_warning(
                                algorithm_name=warning_data.get(
                                    "algorithm", algorithm_name
                                ),
                                warning_message=warning_data.get("message", ""),
                                item_id=rep_id,
                                repetition_number=warning_data.get(
                                    "repetition", returned_rep_num
                                ),
                                execution_context={
                                    "execution_type": "batch_parallel",
                                    "timestamp": warning_data.get("timestamp"),
                                    "dataset_id": execution_context.get(
                                        "dataset_id", "unknown"
                                    ),
                                    "algorithm_params": execution_context.get(
                                        "algorithm_params", {}
                                    ),
                                },
                            )

                        # Log progresso histórico se disponível
                        progress_history = result_data.get("progress_history", [])
                        if progress_history:
                            self._logger.debug(
                                f"Collected {len(progress_history)} progress entries for rep {returned_rep_num}"
                            )
                            # Os dados de progresso são mantidos no resultado para análise posterior

                    # Verificar se houve erro (result_data contém 'error')
                    if isinstance(result_data, dict) and "error" in result_data:
                        self._logger.error(
                            f"Erro na execução do algoritmo {algorithm_name} (rep {returned_rep_num}): {result_data.get('error')}"
                        )

                        # Criar resultado de erro compatível
                        error_result = {
                            "task_id": execution_context.get("task_id", "unknown"),
                            "dataset_id": execution_context.get(
                                "dataset_id", "unknown"
                            ),
                            "algorithm_id": execution_context.get(
                                "algorithm_id", "unknown"
                            ),
                            "algorithm_name": algorithm_name,
                            "repetition": returned_rep_num,
                            "total_repetitions": repetitions,
                            "status": "error",
                            "error": result_data.get("error"),
                            "execution_time": 0.0,
                        }

                        # Notificar monitoramento de erro
                        if monitoring_service:
                            monitoring_service.finish_item(
                                rep_id,
                                False,
                                error_result,
                                result_data.get("error", "Unknown error"),
                            )

                            # Emit algorithm finished event for error in parallel execution
                            monitoring_service.report_algorithm_finished(
                                algorithm_name=algorithm_name,
                                success=False,
                                best_string="",
                                max_distance=-1,
                                execution_time=0.0,
                                metadata={
                                    "error": result_data.get("error", "Unknown error")
                                },
                                repetition_number=returned_rep_num,
                            )

                        results.append(error_result)
                    else:
                        self._logger.debug(
                            f"Algoritmo {algorithm_name} executado com sucesso (rep {returned_rep_num}/{repetitions})"
                        )

                        # Criar resultado de sucesso compatível
                        success_result = {
                            "task_id": execution_context.get("task_id", "unknown"),
                            "dataset_id": execution_context.get(
                                "dataset_id", "unknown"
                            ),
                            "algorithm_id": execution_context.get(
                                "algorithm_id", "unknown"
                            ),
                            "algorithm_name": algorithm_name,
                            "repetition": returned_rep_num,
                            "total_repetitions": repetitions,
                            "status": "success",
                            "best_string": result_data.get("best_string", ""),
                            "distance": result_data.get("distance", 0),
                            "execution_time": result_data.get("execution_time", 0.0),
                            "metadata": result_data.get("metadata", {}),
                            # Adicionar dados coletados de monitoramento
                            "warnings": result_data.get("warnings", []),
                            "progress_history": result_data.get("progress_history", []),
                        }

                        # Notificar monitoramento de conclusão com progresso 100%
                        if monitoring_service:
                            # Enviar progresso 100% para o algoritmo
                            monitoring_service.algorithm_callback(
                                algorithm_name=algorithm_name,
                                progress=100.0,
                                message=f"Completed run {returned_rep_num}/{repetitions}",
                                item_id=rep_id,
                            )
                            monitoring_service.finish_item(rep_id, True, success_result)

                            # Emit algorithm finished event for success in parallel execution
                            monitoring_service.report_algorithm_finished(
                                algorithm_name=algorithm_name,
                                success=True,
                                best_string=result_data.get("best_string", ""),
                                max_distance=result_data.get("distance", 0),
                                execution_time=result_data.get("execution_time", 0.0),
                                metadata=result_data.get("metadata", {}),
                                repetition_number=returned_rep_num,
                            )

                        results.append(success_result)

                except Exception as e:
                    self._logger.error(
                        f"Erro ao processar resultado da repetição {rep_number} de {algorithm_name}: {e}"
                    )

                    error_result = {
                        "task_id": execution_context.get("task_id", "unknown"),
                        "dataset_id": execution_context.get("dataset_id", "unknown"),
                        "algorithm_id": execution_context.get(
                            "algorithm_id", "unknown"
                        ),
                        "algorithm_name": algorithm_name,
                        "repetition": rep_number,
                        "total_repetitions": repetitions,
                        "status": "error",
                        "error": str(e),
                        "execution_time": 0.0,
                    }

                    results.append(error_result)

                    # Notificar monitoramento de erro
                    if monitoring_service:
                        monitoring_service.finish_item(
                            rep_id, False, error_result, str(e)
                        )

        return results

    def _execute_algorithm_repetitions_sequential(
        self,
        algorithm_name: str,
        dataset: Dataset,
        params: Dict[str, Any],
        repetitions: int,
        execution_context: Dict[str, Any],
        monitoring_service=None,
    ) -> List[Dict[str, Any]]:
        """
        Executa repetições de um algoritmo sequencialmente (fallback).

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para processamento
            params: Parâmetros do algoritmo
            repetitions: Número de repetições
            execution_context: Contexto da execução
            monitoring_service: Serviço de monitoramento

        Returns:
            List[Dict[str, Any]]: Lista de resultados das repetições
        """
        results = []

        for rep in range(repetitions):
            rep_id = f"{algorithm_name}_{execution_context.get('dataset_id', 'unknown')}_{rep+1}"

            try:
                # Notificar monitoramento de novo item
                if monitoring_service:
                    # Iniciar item antes da execução
                    monitoring_service.start_item(rep_id, "repetition")
                    monitoring_service.update_item(rep_id, 0.0, "Iniciando")

                # Executar algoritmo
                result = self.execute_single(
                    algorithm_name,
                    dataset,
                    params,
                    monitoring_service,
                    repetition_number=rep + 1,
                )

                # Adicionar informações de contexto
                result.update(
                    {
                        "task_id": execution_context.get("task_id", "unknown"),
                        "dataset_id": execution_context.get("dataset_id", "unknown"),
                        "algorithm_id": execution_context.get(
                            "algorithm_id", "unknown"
                        ),
                        "algorithm_name": algorithm_name,
                        "repetition": rep + 1,
                        "total_repetitions": repetitions,
                        "status": "success",
                    }
                )

                results.append(result)

                # Notificar monitoramento de conclusão
                if monitoring_service:
                    monitoring_service.finish_item(rep_id, True, result)

                self._logger.debug(
                    f"Algoritmo {algorithm_name} executado com sucesso (rep {rep+1}/{repetitions})"
                )

            except Exception as e:
                self._logger.error(
                    f"Erro na execução do algoritmo {algorithm_name} (rep {rep+1}): {e}"
                )

                error_result = {
                    "task_id": execution_context.get("task_id", "unknown"),
                    "dataset_id": execution_context.get("dataset_id", "unknown"),
                    "algorithm_id": execution_context.get("algorithm_id", "unknown"),
                    "algorithm_name": algorithm_name,
                    "repetition": rep + 1,
                    "total_repetitions": repetitions,
                    "status": "error",
                    "error": str(e),
                    "execution_time": 0.0,
                }

                results.append(error_result)

                # Notificar monitoramento de erro
                if monitoring_service:
                    monitoring_service.finish_item(rep_id, False, error_result, str(e))

        return results

    def cleanup(self) -> None:
        """Cleanup orchestrator resources."""
        if self._resource_controller:
            self._resource_controller.cleanup()
            self._resource_controller = None
            self._logger.info("[RESOURCE] Orchestrator cleanup completed")

    def __del__(self):
        """Destructor to ensure cleanup."""
        try:
            self.cleanup()
        except Exception:
            pass  # Ignore cleanup errors during destruction
