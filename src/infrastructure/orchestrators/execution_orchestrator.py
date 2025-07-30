"""
Execution Orchestrator

Centralizes all CSP algorithm execution logic,
for both single executions and batches.
"""

import json
import os
import time
import uuid
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# IMPORTANT: Import algorithms first to load the global_registry
import algorithms

from src.domain import Dataset
from src.domain.errors import AlgorithmExecutionError
from src.infrastructure.logging_config import get_logger
from src.infrastructure.orchestrators.base_orchestrator import BaseOrchestrator
from src.infrastructure.resource_control import create_resource_controller, ResourceController


class ExecutionOrchestrator(BaseOrchestrator):
    """Orchestrator responsible for CSP algorithm execution."""

    def __init__(self, algorithm_registry, dataset_repository, monitoring_service=None, entrez_repository=None):
        """
        Initialize execution orchestrator.

        Args:
            algorithm_registry: Algorithm registry
            dataset_repository: Dataset repository
            monitoring_service: Optional monitoring service
            entrez_repository: Optional Entrez dataset repository
        """
        super().__init__(monitoring_service)
        self._algorithm_registry = algorithm_registry
        self._dataset_repository = dataset_repository
        self._entrez_repository = entrez_repository
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
                self._logger.info("[RESOURCE] Resource controller initialized and limits applied")
            except Exception as e:
                self._logger.warning(f"[RESOURCE] Failed to initialize resource controller: {e}")
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
        self._logger.debug(f"Available algorithms in registry: {list(global_registry.keys())}")
        self._logger.debug(f"Looking for algorithm: {algorithm_name}")

        # Check if algorithm exists
        if algorithm_name not in global_registry:
            raise AlgorithmExecutionError(
                f"Algorithm '{algorithm_name}' not found"
            )

        algorithm_class = global_registry[algorithm_name]
        params = params or {}

        # Apply infrastructure configurations if batch_config available
        if self._current_batch_config:
            infrastructure_config = self._current_batch_config.get("infrastructure", {})
            history_config = infrastructure_config.get("history", {})

            # Injetar parâmetros de histórico se habilitados
            # Support both template format (enabled/frequency) and implementation format (save_history/history_frequency)
            history_enabled = history_config.get("enabled", history_config.get("save_history", False))
            
            if history_enabled:
                params = params.copy()  # Não modificar o original
                params["save_history"] = True
                # Support both template format (frequency) and implementation format (history_frequency)
                frequency = history_config.get("frequency", history_config.get("history_frequency", 1))
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

            # Instancia e executa algoritmo
            algorithm = algorithm_class(
                strings=dataset.sequences, alphabet=dataset.alphabet, **params
            )

            # Configurar callback de progresso se fornecido
            if monitoring_service:

                def progress_callback(message: str):
                    # Usando algorithm_callback da MonitoringInterface
                    monitoring_service.algorithm_callback(
                        algorithm_name=algorithm_name,
                        progress=0.5,  # Progresso genérico, algoritmo pode não informar progresso específico
                        message=message,
                        item_id=execution_id,
                    )

                algorithm.set_progress_callback(progress_callback)

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
                self._logger.error(f"[RESOURCE] Algorithm {algorithm_name} timed out: {e}")
                # Return timeout result
                end_time = time.time()
                return {
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
            except Exception as e:
                if "Resource limit" in str(e) or "Memory limit" in str(e):
                    self._logger.error(f"[RESOURCE] Algorithm {algorithm_name} exceeded resource limits: {e}")
                    # Return resource limit result
                    end_time = time.time()
                    return {
                        "algorithm": algorithm_name,
                        "best_string": "",
                        "max_distance": -1,
                        "execution_time": end_time - start_time,
                        "execution_id": execution_id,
                        "params": params,
                        "metadata": {"error": str(e), "resource_limit_exceeded": True},
                        "dataset": {
                            "size": len(dataset.sequences),
                            "length": len(dataset.sequences[0]) if dataset.sequences else 0,
                            "alphabet": dataset.alphabet,
                        },
                        "status": "resource_limit",
                    }
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

            return result

        except Exception as e:
            # Registra erro
            self._executions[execution_id].update(
                {"status": "failed", "error": str(e), "end_time": time.time()}
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

        # Determinar tipo de monitoramento
        if monitoring_service:
            from src.presentation.monitoring.interfaces import TaskType

            task_type = getattr(TaskType, task_type_str.upper(), TaskType.EXECUTION)
            monitoring_service.start_monitoring(task_type, batch_config)

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
        algorithms_config = batch_config.get("algorithms", [])
        executions = batch_config["execution"]["executions"]

        self._logger.debug(f"Datasets: {len(datasets_config)}")
        self._logger.debug(f"Algorithms: {len(algorithms_config)}")
        self._logger.debug(f"Executions: {len(executions)}")

        # Inicializar dados de monitoramento
        if monitoring_service:
            self._setup_monitoring_data(
                executions, datasets_config, algorithms_config, monitoring_service
            )

        from src.infrastructure import FileDatasetRepository

        dataset_repo = FileDatasetRepository("./datasets")

        results = []
        execution_index = 0

        for execution in executions:
            execution_index += 1
            execution_name = execution.get("name", f"Execution {execution_index}")

            # Iterar sobre configurações de algoritmo
            algorithm_ids = execution["algorithms"]
            for algo_config_idx, algorithm_id in enumerate(algorithm_ids, 1):

                # Resolver configuração do algoritmo
                algorithm_config = next(
                    (a for a in algorithms_config if a["id"] == algorithm_id), None
                )
                if not algorithm_config:
                    self._logger.error(
                        f"Algoritmo com ID '{algorithm_id}' não encontrado"
                    )
                    continue

                # Resolver datasets para esta configuração
                dataset_ids = execution["datasets"]
                for dataset_idx, dataset_id in enumerate(dataset_ids, 1):

                    dataset_config = next(
                        (d for d in datasets_config if d["id"] == dataset_id), None
                    )
                    if not dataset_config:
                        results.append(
                            {
                                "execution_name": execution.get("name", "unknown"),
                                "dataset_id": dataset_id,
                                "status": "error",
                                "error": f"Dataset com ID '{dataset_id}' não encontrado",
                            }
                        )
                        continue

                    # Atualizar informações do dataset no monitoramento
                    if monitoring_service:
                        # Contar algoritmos únicos desta configuração
                        unique_algorithms = set(algorithm_config["algorithms"])

                        # Atualizar hierarquia de dataset (que já inclui execução)
                        from src.presentation.monitoring.interfaces import (
                            ExecutionLevel,
                        )

                        # Obter nome do dataset
                        dataset_name = dataset_config.get("name", dataset_id)

                        # Obter nome da configuração de algoritmo
                        algorithm_config_name = algorithm_config.get(
                            "name", "Algorithms"
                        )

                        monitoring_service.monitor.update_hierarchy(
                            level=ExecutionLevel.DATASET,
                            level_id=f"{dataset_id}_{algorithm_id}",
                            progress=0.0,
                            message=f"Processando dataset {dataset_name}",
                            data={
                                "execution_name": execution_name,
                                "config_index": execution_index,
                                "total_configs": len(executions),
                                "dataset_name": dataset_name,
                                "dataset_index": dataset_idx,
                                "total_datasets": len(dataset_ids),
                                "algorithm_config_name": algorithm_config_name,
                                "algorithm_config_index": algo_config_idx,
                                "total_algorithm_configs": len(algorithm_ids),
                                "total_algorithms": len(unique_algorithms),
                            },
                        )

                    # Carregar dataset e executar algoritmos desta configuração
                    dataset_results = self._execute_dataset_algorithms_for_config(
                        execution,
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
                print(f"[DEBUG] Resultado adicionado com sucesso")

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
        from src.infrastructure import SessionManager

        try:
            session_manager = SessionManager(self._current_batch_config or {})
            session_folder = session_manager.create_session()
            results_dir = Path(session_manager.get_result_dir())
            print(f"📁 Sessão criada: {session_folder}")
            print(f"📁 Salvando resultados parciais em: {results_dir}")
        except Exception as e:
            # Fallback para diretório padrão
            base_dir = Path("./outputs/results")
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            results_dir = base_dir / timestamp
            print(f"📁 Usando diretório fallback: {results_dir}, erro: {e}")

        results_dir.mkdir(parents=True, exist_ok=True)
        self._partial_results_file = str(results_dir / "partial_results.json")

        print(f"💾 Arquivo de resultados parciais: {self._partial_results_file}")

        # Inicializar arquivo com array vazio
        with open(self._partial_results_file, "w", encoding="utf-8") as f:
            json.dump([], f)

        print(f"✅ Sistema de salvamento parcial inicializado")

    def _save_partial_result(self, result: Dict[str, Any]) -> None:
        """Save a partial result to file."""
        if not self._partial_results_file:
            return

        try:
            # Carregar resultados existentes
            if os.path.exists(self._partial_results_file):
                with open(self._partial_results_file, "r", encoding="utf-8") as f:
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
        self, executions, datasets_config, algorithms_config, monitoring_service
    ):
        """Configura dados iniciais de monitoramento."""
        # Calcular totais de datasets considerando todas as execuções
        total_dataset_executions = 0
        for execution in executions:
            total_dataset_executions += len(execution.get("datasets", []))

        # Contar total de algoritmos únicos
        unique_algorithms = set()
        for execution in executions:
            algorithm_ids = execution["algorithms"]
            for algorithm_id in algorithm_ids:
                algorithm_config = next(
                    (a for a in algorithms_config if a["id"] == algorithm_id), None
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
        execution_name,
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
        from src.presentation.monitoring.interfaces import ExecutionLevel

        # Obter total de execuções do contexto (será passado pelo orchestrator)
        # Por enquanto, usando o index como fallback
        total_executions = execution_index  # Isso será melhorado no orchestrator

        monitoring_service.monitor.update_hierarchy(
            level=ExecutionLevel.EXECUTION,
            level_id=execution_name,
            progress=0.0,
            message=f"Iniciando execução {execution_name}",
            data={
                "execution_name": execution_name,
                "config_index": execution_index,
                "total_configs": total_executions,
            },
        )

    def _execute_dataset_algorithms_for_config(
        self,
        execution,
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
                return [{
                    "execution_name": execution.get("name", "unknown"),
                    "dataset_id": dataset_id,
                    "status": "error",
                    "error": error_msg,
                    "execution_time": 0.0,
                }]
            
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
                    raise ValueError("Field 'query' is required for entrez dataset type")
                
                db = params.get("db", "nucleotide")
                retmax = params.get("retmax", 20)
                
                # Fetch dataset from NCBI
                # Extract only additional parameters, avoiding duplicates
                additional_params = {k: v for k, v in params.items() 
                                   if k not in ['query', 'db', 'retmax']}
                
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
                        **metadata
                    }
                )
                
                self._logger.info(
                    "Created Entrez dataset: n=%d, L=%d, query='%s'",
                    len(sequences), len(sequences[0]) if sequences else 0, query
                )
            else:
                raise ValueError(f"Unsupported dataset type: {dataset_type}")

            self._logger.info(
                f"Dataset {dataset_id} loaded: {len(dataset.sequences)} sequences"
            )

            # Execute algorithms for this configuration
            algorithm_names = algorithm_config["algorithms"]
            algorithm_params = algorithm_config.get("algorithm_params", {})
            repetitions = execution.get("repetitions", 1)

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
                        "execution_name": execution.get("name", "unknown"),
                        "dataset_id": dataset_id,
                        "algorithm_id": algorithm_config["id"],
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
                    "execution_name": execution.get("name", "unknown"),
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
        dataset: Dataset,
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
            dataset: Dataset para processamento
            params: Parâmetros do algoritmo
            execution_context: Contexto da execução (nomes, IDs, etc.)
            rep_number: Número da repetição (1-based)
            total_repetitions: Total de repetições

        Returns:
            Dict[str, Any]: Resultado da execução com contexto
        """
        try:
            # Executar algoritmo (sem monitoring_service pois não é thread-safe)
            result = self.execute_single(
                algorithm_name, dataset, params, monitoring_service=None
            )

            # Sanitizar metadados para serialização
            if "metadata" in result:
                result["metadata"] = self._sanitize_metadata_for_process(result["metadata"])

            # Adicionar informações de contexto
            result.update(
                {
                    "execution_name": execution_context.get(
                        "execution_name", "unknown"
                    ),
                    "dataset_id": execution_context.get("dataset_id", "unknown"),
                    "algorithm_id": execution_context.get("algorithm_id", "unknown"),
                    "algorithm_name": algorithm_name,
                    "repetition": rep_number,
                    "total_repetitions": total_repetitions,
                    "status": "success",
                }
            )

            return result

        except Exception as e:
            error_result = {
                "execution_name": execution_context.get("execution_name", "unknown"),
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
        import types
        import pickle
        
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
                        sanitized[key] = f"<function '{getattr(value, '__name__', str(value))}'>"
                    elif hasattr(value, '__dict__') and hasattr(value, '__class__'):
                        # Para objetos complexos, extrai atributos serializáveis
                        try:
                            class_name = value.__class__.__name__
                            repr_str = str(value)
                            sanitized[key] = {
                                '__class__': class_name,
                                '__repr__': repr_str
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
                            "execution_name": execution.get("name", "unknown"),
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
                            "execution_name": execution.get("name", "unknown"),
                            "dataset_id": dataset_id,
                            "algorithm_id": algorithm_id,
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
                    "execution_name": execution.get("name", "unknown"),
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
            args_list.append(
                (
                    algorithm_name,
                    dataset,
                    params,
                    execution_context,
                    rep + 1,  # 1-based
                    repetitions,
                )
            )

        # Executar em paralelo
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submeter todas as tarefas
            future_to_rep = {}
            for i, args in enumerate(args_list):
                future = executor.submit(self._execute_single_repetition, *args)
                future_to_rep[future] = i + 1

            # Coletar resultados conforme completam
            for future in as_completed(future_to_rep):
                rep_number = future_to_rep[future]
                rep_id = f"{algorithm_name}_{execution_context.get('dataset_id', 'unknown')}_{rep_number}"

                try:
                    # Inicializar monitoramento se disponível
                    if monitoring_service:
                        from src.presentation.monitoring.interfaces import (
                            HierarchicalContext,
                        )

                        context = HierarchicalContext(
                            dataset_id=execution_context.get("dataset_id", "unknown"),
                            algorithm_id=algorithm_name,
                            repetition_id=f"{rep_number}/{repetitions}",
                        )
                        # Iniciar item antes da execução
                        monitoring_service.start_item(rep_id, "repetition", context)

                    # Obter resultado
                    result = future.result()

                    # Verificar se houve erro
                    if result.get("status") == "error":
                        self._logger.error(
                            f"Erro na execução do algoritmo {algorithm_name} (rep {rep_number}): {result.get('error')}"
                        )

                        # Notificar monitoramento de erro
                        if monitoring_service:
                            monitoring_service.finish_item(
                                rep_id,
                                False,
                                result,
                                result.get("error", "Unknown error"),
                            )
                    else:
                        self._logger.debug(
                            f"Algoritmo {algorithm_name} executado com sucesso (rep {rep_number}/{repetitions})"
                        )

                        # Notificar monitoramento de conclusão
                        if monitoring_service:
                            monitoring_service.finish_item(rep_id, True, result)

                    results.append(result)

                except Exception as e:
                    self._logger.error(
                        f"Erro ao processar resultado da repetição {rep_number} de {algorithm_name}: {e}"
                    )

                    error_result = {
                        "execution_name": execution_context.get(
                            "execution_name", "unknown"
                        ),
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
                    from src.presentation.monitoring.interfaces import (
                        HierarchicalContext,
                    )

                    context = HierarchicalContext(
                        dataset_id=execution_context.get("dataset_id", "unknown"),
                        algorithm_id=algorithm_name,
                        repetition_id=f"{rep+1}/{repetitions}",
                    )
                    # Iniciar item antes da execução
                    monitoring_service.start_item(rep_id, "repetition", context)
                    monitoring_service.update_item(rep_id, 0.0, "Iniciando", context)

                # Executar algoritmo
                result = self.execute_single(
                    algorithm_name, dataset, params, monitoring_service
                )

                # Adicionar informações de contexto
                result.update(
                    {
                        "execution_name": execution_context.get(
                            "execution_name", "unknown"
                        ),
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
                    "execution_name": execution_context.get(
                        "execution_name", "unknown"
                    ),
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
