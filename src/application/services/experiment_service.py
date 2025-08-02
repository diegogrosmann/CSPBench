"""
CSPBench Experiment Service

Orchestrates the execution of experiments, optimizations, and sensitivity analyses.
Implements application use cases without infrastructure dependency.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

# IMPORTANT: Import algorithms first to load the global_registry
from src.application.ports import (
    AlgorithmRegistry,
    DatasetRepository,
    EntrezDatasetRepository,
    ExecutorPort,
    ExportPort,
)
from src.application.services.config_parser import (
    BatchConfig,
    ConfigParser,
    SystemConfig,
)
from src.domain import (
    AlgorithmExecutionError,
    AlgorithmNotFoundError,
    BatchConfigurationError,
    BatchExecutionError,
    Dataset,
    DatasetNotFoundError,
    OptimizationConfigurationError,
    OptimizationExecutionError,
    SensitivityConfigurationError,
    SensitivityExecutionError,
)
from src.domain.dataset import SyntheticDatasetGenerator
from src.infrastructure.logging_config import LoggerConfig, get_logger
from src.infrastructure.session_manager import SessionManager
from src.presentation.monitoring.interfaces import TaskType


class ExperimentService:
    """
    Main service for CSP experiment execution.

    Coordinates execution of different types of experiments through defined ports,
    keeping application logic separate from infrastructure.
    """

    def __init__(
        self,
        dataset_repo: DatasetRepository,
        exporter: ExportPort,
        executor: ExecutorPort,
        algo_registry: AlgorithmRegistry,
        entrez_repo: Optional[EntrezDatasetRepository] = None,
        monitoring_service: Optional[Any] = None,
        session_manager: Optional[Any] = None,
    ):
        """
        Initialize the experiments service.

        Args:
            dataset_repo: Dataset repository
            exporter: Port for results exportation
            executor: Port for algorithm execution
            algo_registry: Algorithm registry
            entrez_repo: Entrez dataset repository (optional)
            monitoring_service: Monitoring service (optional)
            session_manager: Session manager for paths and sessions (optional)
        """
        self._dataset_repo = dataset_repo
        self._exporter = exporter
        self._executor = executor
        self._algo_registry = algo_registry
        self._entrez_repo = entrez_repo
        self._monitoring_service = monitoring_service
        self._session_manager = session_manager
        self._logger = get_logger(__name__)

    def _update_batch_logging_with_config(self, logging_config) -> None:
        """Update logging configuration based on parsed logging config."""
        try:
            if logging_config.level:
                LoggerConfig.set_level(logging_config.level)
                self._logger.info(f"Log level updated to: {logging_config.level}")

            # TODO: Implement additional logging configurations like file output, formatters, etc.
            # This can be extended as needed when more detailed logging control is required

        except Exception as e:
            self._logger.warning(f"Error updating batch logging configuration: {e}")

    def _update_batch_logging(self, batch_config: Dict[str, Any]) -> None:
        """Update logging configuration based on specific batch."""
        try:
            # Check logging configuration in current format
            if "logging" in batch_config:
                log_config = batch_config["logging"]

                # Update level if specified
                if "level" in log_config:
                    new_level = log_config["level"]
                    LoggerConfig.set_level(new_level)
                    self._logger.info(f"Log level updated to: {new_level}")

            # Support legacy format as well
            elif "advanced" in batch_config and "logs" in batch_config["advanced"]:
                log_config = batch_config["advanced"]["logs"]

                # Check if logs are enabled
                if not log_config.get("enabled", True):
                    return

                # Update level if specified (environment variable takes precedence)
                if "level" in log_config:
                    import os

                    new_level = os.getenv("LOG_LEVEL", log_config["level"])
                    LoggerConfig.set_level(new_level)
                    self._logger.info(f"Log level updated to: {new_level}")
        except Exception as e:
            self._logger.warning(f"Error updating batch logging configuration: {e}")

    def _convert_batch_config_to_dict(self, batch_config) -> Dict[str, Any]:
        """Convert BatchConfig object back to dict for compatibility."""
        result = {}

        # Metadata section
        result["metadata"] = {
            "name": batch_config.metadata.name,
            "description": batch_config.metadata.description,
            "author": batch_config.metadata.author,
            "version": batch_config.metadata.version,
            "creation_date": batch_config.metadata.creation_date,
            "tags": batch_config.metadata.tags,
        }

        # Task section
        result["task"] = {"type": batch_config.task.type}

        # Datasets section
        result["datasets"] = []
        for dataset in batch_config.datasets:
            result["datasets"].append(
                {
                    "id": dataset.id,
                    "name": dataset.name,
                    "type": dataset.type,
                    "parameters": dataset.parameters,
                }
            )

        # Algorithms section
        result["algorithms"] = []
        for algorithm in batch_config.algorithms:
            result["algorithms"].append(
                {
                    "id": algorithm.id,
                    "name": algorithm.name,
                    "description": algorithm.description,
                    "algorithms": algorithm.algorithms,
                    "algorithm_params": algorithm.algorithm_params,
                }
            )

        # Task-specific sections
        if batch_config.execution:
            result["execution"] = {"executions": []}
            for exec_config in batch_config.execution.get("executions", []):
                if hasattr(exec_config, "name"):  # ExecutionConfig object
                    result["execution"]["executions"].append(
                        {
                            "name": exec_config.name,
                            "datasets": exec_config.datasets,
                            "algorithms": exec_config.algorithms,
                            "repetitions": exec_config.repetitions,
                        }
                    )
                else:  # Already a dict
                    result["execution"]["executions"].append(exec_config)

        if batch_config.optimization:
            result["optimization"] = batch_config.optimization
        if batch_config.sensitivity:
            result["sensitivity"] = batch_config.sensitivity

        return result

    def run_batch(self, batch_cfg: str) -> Dict[str, Any]:
        """
        Execute batch experiments from configuration.

        Args:
            batch_cfg: Path or content of batch configuration

        Returns:
            Dict[str, Any]: Consolidated batch results

        Raises:
            BatchConfigurationError: If invalid configuration
            BatchExecutionError: If execution error
            DatasetNotFoundError: If dataset not found
            AlgorithmNotFoundError: If algorithm not found
        """
        try:
            self._logger.info(f"Starting batch execution: {batch_cfg}")

            # Parse batch configuration using new parser
            parsed_config = ConfigParser.parse_config(batch_cfg)

            # Extract metadata
            metadata = parsed_config.metadata
            self._logger.info(f"Batch loaded: {metadata.name}")

            # Extract all configuration sections
            infrastructure_config = parsed_config.infrastructure
            export_config = parsed_config.export
            plots_config = parsed_config.plots
            monitoring_config = parsed_config.monitoring
            logging_config = parsed_config.logging
            system_config = parsed_config.system
            resources_config = parsed_config.resources

            # Update existing SessionManager with batch-specific configuration if needed
            # Don't create a new one to preserve the current session
            if not hasattr(self, "_session_manager") or not self._session_manager:
                # Only create new SessionManager if none exists
                batch_dict = self._convert_batch_config_to_dict(parsed_config)
                self._session_manager = SessionManager(batch_dict)
            # If SessionManager already exists, keep using it (preserves current session)

            # Apply global system configurations
            self._apply_system_config(system_config, parsed_config)

            # Update logging configuration if specified in batch
            self._update_batch_logging_with_config(logging_config)

            # Determine task type from parsed config
            task_type = parsed_config.task.type
            self._logger.info(f"Task type detected: {task_type}")

            # Configure monitoring if enabled
            if (
                monitoring_config
                and monitoring_config.enabled
                and self._monitoring_service
            ):
                self._configure_monitoring(monitoring_config)

            # Process according to task type
            if task_type == "execution":
                # Convert BatchConfig back to dict for compatibility
                batch_dict = self._convert_batch_config_to_dict(parsed_config)
                results = self._process_execution_batch(
                    batch_dict, resources_config.__dict__ if resources_config else {}
                )
            elif task_type == "optimization":
                batch_dict = self._convert_batch_config_to_dict(parsed_config)
                results = self._process_optimization_batch(
                    batch_dict, resources_config.__dict__ if resources_config else {}
                )
            elif task_type == "sensitivity":
                batch_dict = self._convert_batch_config_to_dict(parsed_config)
                results = self._process_sensitivity_batch(
                    batch_dict, resources_config.__dict__ if resources_config else {}
                )
            else:
                raise BatchConfigurationError(f"Unsupported task type: {task_type}")

            # Consolidate results
            consolidated_results = self._consolidate_batch_results(results["results"])

            total_experiments = len(results["results"])
            successful = consolidated_results.get("summary", {}).get("successful", 0)
            failed = consolidated_results.get("summary", {}).get("failed", 0)
            self._logger.info(
                f"Batch completed: {total_experiments} experiments, {successful} successful, {failed} failed"
            )

            # Export results using new export configuration
            if export_config.enabled:
                export_data = {
                    "batch_summary": consolidated_results,
                    "detailed_results": results.get("results", []),
                }

                # Add task-specific data
                if task_type == "sensitivity":
                    export_data["sensitivity_summaries"] = results.get(
                        "sensitivity_summaries", []
                    )

                export_path = self._export_batch_results_with_config(
                    export_data, export_config, task_type
                )
                self._logger.info(f"Results exported to: {export_path}")

            # Generate plots if enabled
            if plots_config.enabled:
                self._generate_plots_with_config(results, plots_config, task_type)

            return consolidated_results

        except Exception as e:
            self._logger.error(f"Erro durante execução de batch: {e}")
            if isinstance(
                e,
                (BatchConfigurationError, DatasetNotFoundError, AlgorithmNotFoundError),
            ):
                raise
            else:
                raise BatchExecutionError(
                    f"Erro inesperado durante execução: {e}"
                ) from e

    def optimize(self, opt_cfg: str) -> Dict[str, Any]:
        """
        Execute hyperparameter optimization.

        Args:
            opt_cfg: Path or content of optimization configuration

        Returns:
            Dict[str, Any]: Optimization results

        Raises:
            OptimizationConfigurationError: If configuration is invalid
            OptimizationExecutionError: If execution error
            DatasetNotFoundError: If dataset not found
            AlgorithmNotFoundError: If algorithm not found
        """
        try:
            # Parse optimization configuration
            opt_config = self._parse_optimization_config(opt_cfg)

            # Validate configuration
            self._validate_optimization_config(opt_config)

            # Load dataset
            dataset = self._load_dataset(opt_config["dataset"])

            # Check if algorithm exists
            algorithm_name = opt_config["algorithm"]
            if not self._algo_registry.algorithm_exists(algorithm_name):
                raise AlgorithmNotFoundError(f"Algorithm not found: {algorithm_name}")

            # Execute optimization
            results = self._executor.execute_optimization(
                algorithm_name,
                dataset,
                opt_config,
                monitoring_service=self._monitoring_service,
            )

            # Exportar resultados se configurado
            if opt_config.get("export", {}).get("enabled", False):
                export_config = opt_config["export"]

                # Use session manager to determine export path
                if hasattr(self, "_session_manager") and self._session_manager:
                    destination = self._session_manager.get_session_result_path()
                else:
                    import os

                    base_output_dir = os.getenv("OUTPUT_BASE_DIRECTORY", "outputs")
                    destination = export_config.get(
                        "destination", f"{base_output_dir}/optimization_results"
                    )

                self._exporter.export_optimization_results(
                    results,
                    destination,
                )

            return results

        except Exception as e:
            if isinstance(
                e,
                (
                    OptimizationConfigurationError,
                    DatasetNotFoundError,
                    AlgorithmNotFoundError,
                ),
            ):
                raise
            raise OptimizationExecutionError(
                f"Erro na execução da otimização: {e}"
            ) from e

    def sensitivity(self, sens_cfg: str) -> Dict[str, Any]:
        """
        Execute parameter sensitivity analysis.

        Args:
            sens_cfg: Path or content of analysis configuration

        Returns:
            Dict[str, Any]: Sensitivity analysis results

        Raises:
            SensitivityConfigurationError: If configuration is invalid
            SensitivityExecutionError: If execution error
            DatasetNotFoundError: If dataset not found
            AlgorithmNotFoundError: If algorithm not found
        """
        try:
            # Parse sensitivity configuration
            sens_config = self._parse_sensitivity_config(sens_cfg)

            # Validate configuration
            self._validate_sensitivity_config(sens_config)

            # Load dataset
            dataset = self._load_dataset(sens_config["dataset"])

            # Check if algorithm exists
            algorithm_name = sens_config["algorithm"]
            if not self._algo_registry.algorithm_exists(algorithm_name):
                raise AlgorithmNotFoundError(f"Algorithm not found: {algorithm_name}")

            # Execute sensitivity analysis
            results = self._executor.execute_sensitivity_analysis(
                algorithm_name, dataset, sens_config
            )

            # Exportar resultados se configurado
            if sens_config.get("export", {}).get("enabled", False):
                import os

                base_output_dir = os.getenv("OUTPUT_BASE_DIRECTORY", "outputs")
                export_config = sens_config["export"]
                self._exporter.export_results(
                    results,
                    export_config.get("format", "json"),
                    export_config.get(
                        "destination", f"{base_output_dir}/sensitivity_results"
                    ),
                )

            return results

        except Exception as e:
            if isinstance(
                e,
                (
                    SensitivityConfigurationError,
                    DatasetNotFoundError,
                    AlgorithmNotFoundError,
                ),
            ):
                raise
            raise SensitivityExecutionError(
                f"Erro na execução da análise de sensibilidade: {e}"
            ) from e

    def run_single_experiment(
        self,
        algorithm_name: str,
        dataset_id: str,
        params: Optional[Dict[str, Any]] = None,
        timeout: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Execute single experiment.

        Args:
            algorithm_name: Algorithm name
            dataset_id: Dataset identifier
            params: Algorithm parameters
            timeout: Timeout in seconds

        Returns:
            Dict[str, Any]: Experiment results
        """
        # Check if algorithm exists
        if not self._algo_registry.algorithm_exists(algorithm_name):
            raise AlgorithmNotFoundError(f"Algorithm not found: {algorithm_name}")

        # Load dataset
        dataset = self._load_dataset(dataset_id)

        # Create batch configuration for single execution (legacy format)
        batch_config = {
            "experiments": [
                {
                    "algorithm": algorithm_name,
                    "dataset": dataset_id,
                    "params": params or {},
                }
            ],
        }

        # Add timeout if provided
        if timeout is not None:
            batch_config["experiments"][0]["timeout"] = timeout

        # Execute as batch and return first result
        results = self._executor.execute_batch(batch_config)
        if results:
            result = results[0]
            # Garantir que timeout está nos metadados se foi fornecido
            if timeout is not None and "metadata" in result:
                result["metadata"]["timeout"] = timeout
            return result
        else:
            raise AlgorithmExecutionError(f"Erro na execução de {algorithm_name}")

    def list_available_algorithms(self) -> List[str]:
        """List available algorithms."""
        return self._algo_registry.list_algorithms()

    def list_available_datasets(self) -> List[str]:
        """List available datasets."""
        return self._dataset_repo.list_available()

    def get_algorithm_info(self, algorithm_name: str) -> Dict[str, Any]:
        """Get algorithm information."""
        if not self._algo_registry.algorithm_exists(algorithm_name):
            raise AlgorithmNotFoundError(f"Algorithm not found: {algorithm_name}")

        return self._algo_registry.get_algorithm_metadata(algorithm_name)

    # Métodos privados auxiliares

    def _process_execution_batch(
        self, batch_config: Dict[str, Any], resources_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Process execution batch using existing executor.

        Args:
            batch_config: Batch configuration
            resources_config: Resource configuration

        Returns:
            Dict[str, Any]: Execution results
        """
        self._logger.info("Processando batch de execução")

        # Iniciar monitoramento se disponível
        if self._monitoring_service:
            batch_name = batch_config.get("metadata", {}).get("name", "Execution")
            self._monitoring_service.start_monitoring(TaskType.EXECUTION, batch_name)

        try:
            # Include resources configuration in batch config
            batch_config_with_resources = batch_config.copy()
            batch_config_with_resources["resources"] = resources_config

            results = self._executor.execute_batch(
                batch_config_with_resources, self._monitoring_service
            )

            # Finalizar monitoramento
            if self._monitoring_service:
                self._monitoring_service.finish_monitoring({"results": results})

            return {"results": results}
        except Exception as e:
            # Mostrar erro no monitoramento
            if self._monitoring_service:
                self._monitoring_service.show_error(str(e))
                self._monitoring_service.finish_monitoring({})
            raise

    def _process_optimization_batch(
        self, batch_config: Dict[str, Any], resources_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Process optimization batch with support for multiple configurations.

        Args:
            batch_config: Batch configuration

        Returns:
            Dict[str, Any]: Consolidated optimization results
        """
        self._logger.info("Processando batch de otimização")

        # Iniciar monitoramento se disponível
        if self._monitoring_service:
            batch_name = batch_config.get("metadata", {}).get("name", "Optimization")
            self._monitoring_service.start_monitoring(TaskType.OPTIMIZATION, batch_name)

        try:
            # Use o método de compatibilidade do backup
            from .config_parser_backup import ConfigurationParser

            optimization_configs = ConfigurationParser.parse_optimization_configs(
                batch_config
            )

            all_results = []
            optimization_summaries = []

            # Calcular totais para monitoramento
            total_optimizations = len(optimization_configs)

            # Calculate total executions (datasets × algorithms per configuration)
            total_executions = 0
            for opt_config in optimization_configs:
                # Resolve algorithms to calculate total
                algorithm_names, _ = self._resolve_algorithm_configuration(
                    getattr(
                        opt_config,
                        "target_algorithm",
                        getattr(opt_config, "algorithm", None),
                    ),
                    batch_config,
                )
                target_datasets = getattr(
                    opt_config, "target_datasets", getattr(opt_config, "datasets", [])
                )
                total_executions += len(target_datasets) * len(algorithm_names)

            total_datasets = sum(
                len(getattr(config, "target_datasets", getattr(config, "datasets", [])))
                for config in optimization_configs
            )

            current_optimization_index = 0
            current_execution_index = 0

            for opt_config in optimization_configs:
                current_optimization_index += 1
                self._logger.info(f"Executing optimization: {opt_config.name}")

                # Resolve algorithms and base parameters by ID
                algorithm_names, all_algorithm_params = (
                    self._resolve_algorithm_configuration(
                        getattr(
                            opt_config,
                            "target_algorithm",
                            getattr(opt_config, "algorithm", None),
                        ),
                        batch_config,
                    )
                )

                # Process each algorithm individually
                for algorithm_name in algorithm_names:
                    current_dataset_index = 0

                    # Get base parameters specific for this algorithm
                    base_params = all_algorithm_params.get(algorithm_name, {})

                    # Process each dataset in configuration
                    for dataset_id in getattr(
                        opt_config,
                        "target_datasets",
                        getattr(opt_config, "datasets", []),
                    ):
                        current_dataset_index += 1
                        current_execution_index += 1
                        self._logger.info(f"Processando dataset: {dataset_id}")

                        try:
                            # Resolve dataset configuration and load it
                            dataset_config = self._resolve_dataset_config(
                                dataset_id, batch_config.get("datasets", [])
                            )

                            # Load dataset based on configuration
                            if dataset_config["type"] == "file":
                                filename = dataset_config["parameters"]["filename"]
                                dataset = self._load_dataset(filename)
                            else:
                                # For synthetic datasets, create using existing method
                                dataset = self._create_dataset_from_config(
                                    dataset_config
                                )

                            self._logger.info(
                                f"Dataset {dataset_id} loaded: {len(dataset.sequences)} sequences"
                            )

                            # Extrair configurações de recursos
                            from .config_parser_backup import ConfigurationParser

                            resources_config = (
                                ConfigurationParser.parse_resources_config(batch_config)
                            )

                            # Processar parâmetros de otimização - nova estrutura
                            optimization_params = {}
                            if algorithm_name in opt_config.parameters:
                                optimization_params = opt_config.parameters[
                                    algorithm_name
                                ]
                            else:
                                # Fallback para estrutura antiga (compatibilidade)
                                optimization_params = opt_config.parameters

                            # Prepare configuration for executor
                            executor_config = {
                                "study_name": f"{opt_config.study_name}_{algorithm_name}_{dataset_id}",
                                "direction": opt_config.direction,
                                "n_trials": opt_config.trials,
                                "timeout_per_trial": opt_config.timeout_per_trial,
                                "parameters": optimization_params,
                                "base_params": base_params,  # Base parameters from configuration
                                "optuna_config": opt_config.optuna_config or {},
                                "resources": resources_config,  # Include resource configurations
                                "internal_jobs": resources_config.get(
                                    "internal_jobs", 4
                                ),  # Internal parallelism
                            }

                            # Execute optimization
                            optimization_results = self._executor.execute_optimization(
                                algorithm_name,  # Usar nome do algoritmo resolvido
                                dataset,
                                executor_config,
                                self._monitoring_service,
                                config_index=current_execution_index,
                                total_configs=total_executions,
                                dataset_index=current_dataset_index,
                                total_datasets=len(
                                    getattr(
                                        opt_config,
                                        "target_datasets",
                                        getattr(opt_config, "datasets", []),
                                    )
                                ),
                                dataset_name=dataset_id,  # Passar nome original do dataset
                            )

                            all_results.append(optimization_results)

                            # Criar sumário para esta otimização
                            optimization_summaries.append(
                                {
                                    "optimization_name": opt_config.name,
                                    "algorithm": algorithm_name,  # Usar nome do algoritmo resolvido
                                    "dataset": dataset_id,
                                    "best_value": optimization_results.get(
                                        "best_value"
                                    ),
                                    "best_params": optimization_results.get(
                                        "best_params"
                                    ),
                                    "n_trials": optimization_results.get("n_trials"),
                                    "total_time": optimization_results.get(
                                        "total_time"
                                    ),
                                }
                            )

                            self._logger.info(f"Otimização concluída para {dataset_id}")

                        except Exception as e:
                            self._logger.error(
                                f"Erro na otimização de {dataset_id}: {e}"
                            )
                            # Adicionar resultado de erro
                            error_result = {
                                "optimization_name": opt_config.name,
                                "algorithm": algorithm_name,  # Usar nome do algoritmo atual
                                "dataset": dataset_id,
                                "error": str(e),
                                "status": "failed",
                            }
                            all_results.append(error_result)
                            optimization_summaries.append(error_result)

            return {
                "results": all_results,
                "optimization_summaries": optimization_summaries,
            }
        except Exception as e:
            # Mostrar erro no monitoramento
            if self._monitoring_service:
                self._monitoring_service.show_error(str(e))
                self._monitoring_service.finish_monitoring({})
            raise

    def _process_sensitivity_batch(
        self, batch_config: Dict[str, Any], resources_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Processa batch de análise de sensibilidade com suporte a múltiplas configurações.

        Args:
            batch_config: Configuração do batch

        Returns:
            Dict[str, Any]: Resultados consolidados da análise de sensibilidade
        """
        self._logger.info("Processando batch de análise de sensibilidade")

        # Iniciar monitoramento se disponível
        if self._monitoring_service:
            batch_name = batch_config.get("metadata", {}).get(
                "name", "Sensitivity Analysis"
            )
            self._monitoring_service.start_monitoring(TaskType.SENSITIVITY, batch_name)

        try:
            # Use o método de compatibilidade do backup
            from .config_parser_backup import ConfigurationParser

            sensitivity_configs = ConfigurationParser.parse_sensitivity_configs(
                batch_config
            )

            all_results = []
            sensitivity_summaries = []

            for sens_config in sensitivity_configs:
                self._logger.info(f"Executing sensitivity analysis: {sens_config.name}")

                # Processar cada dataset na configuração
                for dataset_id in getattr(
                    sens_config, "target_datasets", getattr(sens_config, "datasets", [])
                ):
                    self._logger.info(f"Processando dataset: {dataset_id}")

                    try:
                        # Resolver configuração do dataset e carregá-lo
                        dataset_config = self._resolve_dataset_config(
                            dataset_id, batch_config.get("datasets", [])
                        )

                        # Load dataset based on configuration
                        if dataset_config["type"] == "file":
                            filename = dataset_config["parameters"]["filename"]
                            dataset = self._load_dataset(filename)
                        else:
                            # For synthetic datasets, create using existing method
                            dataset = self._create_dataset_from_config(dataset_config)

                        self._logger.info(
                            f"Dataset {dataset_id} loaded: {len(dataset.sequences)} sequences"
                        )

                        # Extrair configurações de recursos
                        from .config_parser_backup import ConfigurationParser

                        resources_config = ConfigurationParser.parse_resources_config(
                            batch_config
                        )

                        # Preparar configuração para o executor
                        # Resolver algoritmo e parâmetros base por ID
                        algorithm_names, base_params = (
                            self._resolve_algorithm_configuration(
                                getattr(
                                    sens_config,
                                    "target_algorithm",
                                    getattr(sens_config, "algorithm", None),
                                ),
                                batch_config,
                            )
                        )

                        # Para sensitivity analysis, usar apenas o primeiro algoritmo
                        algorithm_name = (
                            algorithm_names[0]
                            if isinstance(algorithm_names, list)
                            else algorithm_names
                        )

                        print(f"DEBUG: Resolved algorithm_names = {algorithm_names}")
                        print(f"DEBUG: Using algorithm_name = {algorithm_name}")
                        print(f"DEBUG: base_params = {base_params}")
                        print("DEBUG: sens_config attributes:")
                        print(
                            f"  - method: {getattr(sens_config, 'method', 'NOT_FOUND')}"
                        )
                        print(
                            f"  - analysis_method: {getattr(sens_config, 'analysis_method', 'NOT_FOUND')}"
                        )
                        print(
                            f"  - samples: {getattr(sens_config, 'samples', 'NOT_FOUND')}"
                        )
                        print(
                            f"  - n_samples: {getattr(sens_config, 'n_samples', 'NOT_FOUND')}"
                        )
                        print(
                            f"  - repetitions: {getattr(sens_config, 'repetitions', 'NOT_FOUND')}"
                        )
                        print(
                            f"  - repetitions_per_sample: {getattr(sens_config, 'repetitions_per_sample', 'NOT_FOUND')}"
                        )

                        # Processar parâmetros de sensibilidade - nova estrutura
                        sensitivity_params = {}
                        if algorithm_name in sens_config.parameters:
                            sensitivity_params = sens_config.parameters[algorithm_name]
                        else:
                            # Fallback para estrutura antiga (compatibilidade)
                            sensitivity_params = sens_config.parameters

                        print(f"DEBUG: sensitivity_params = {sensitivity_params}")

                        executor_config = {
                            "analysis_method": getattr(
                                sens_config,
                                "method",
                                getattr(sens_config, "analysis_method", "morris"),
                            ),
                            "n_samples": getattr(
                                sens_config,
                                "samples",
                                getattr(sens_config, "n_samples", 1000),
                            ),
                            "repetitions_per_sample": getattr(
                                sens_config,
                                "repetitions",
                                getattr(sens_config, "repetitions_per_sample", 3),
                            ),
                            "parameters": sensitivity_params,
                            "base_params": base_params,  # Parâmetros base da configuração
                            "output_metrics": sens_config.output_metrics,
                            "method_config": getattr(
                                sens_config, "method_config", {}
                            ),  # Use getattr with default
                            "resources": resources_config,  # Incluir configurações de recursos
                            "internal_jobs": resources_config.get(
                                "internal_jobs", 4
                            ),  # Paralelismo interno
                        }

                        print(f"DEBUG: executor_config = {executor_config}")

                        # Executar análise de sensibilidade
                        sensitivity_results = (
                            self._executor.execute_sensitivity_analysis(
                                algorithm_name, dataset, executor_config
                            )
                        )

                        all_results.append(sensitivity_results)

                        # Criar sumário para esta análise
                        sensitivity_summaries.append(
                            {
                                "analysis_name": sens_config.name,
                                "algorithm": algorithm_name,  # Usar nome do algoritmo resolvido
                                "dataset": dataset_id,
                                "n_samples": sensitivity_results.get("n_samples"),
                                "parameters_analyzed": list(sensitivity_params.keys()),
                                "total_time": sensitivity_results.get("total_time"),
                            }
                        )

                        self._logger.info(
                            f"Análise de sensibilidade concluída para {dataset_id}"
                        )

                    except Exception as e:
                        self._logger.error(
                            f"Erro na análise de sensibilidade de {dataset_id}: {e}"
                        )
                        # Adicionar resultado de erro
                        error_result = {
                            "analysis_name": sens_config.name,
                            "algorithm": getattr(
                                sens_config,
                                "target_algorithm",
                                getattr(sens_config, "algorithm", None),
                            ),  # Manter ID original no erro
                            "dataset": dataset_id,
                            "error": str(e),
                            "status": "failed",
                        }
                        all_results.append(error_result)
                        sensitivity_summaries.append(error_result)

        except Exception as e:
            # Mostrar erro no monitoramento
            if self._monitoring_service:
                self._monitoring_service.show_error(str(e))
                self._monitoring_service.finish_monitoring({})
            raise

        # Finalizar monitoramento
        if self._monitoring_service:
            self._monitoring_service.finish_monitoring({})

        # Retornar resultados consolidados
        return {
            "results": all_results,
            "sensitivity_summaries": sensitivity_summaries,
            "summary": {
                "total_analyses": len(sensitivity_summaries),
                "successful": len(
                    [r for r in all_results if r.get("status") != "failed"]
                ),
                "failed": len([r for r in all_results if r.get("status") == "failed"]),
            },
        }

    def _parse_batch_config(self, batch_cfg: str) -> Dict[str, Any]:
        """Parseia configuração de batch."""
        import json
        from pathlib import Path

        import yaml

        # Se não é string, retorna como está (já é dict)
        if not isinstance(batch_cfg, str):
            return batch_cfg

        # Se é caminho de arquivo, carrega
        file_path = Path(batch_cfg)
        if file_path.exists():
            with open(file_path, encoding="utf-8") as f:
                if file_path.suffix.lower() in [".yaml", ".yml"]:
                    config = yaml.safe_load(f)
                elif file_path.suffix.lower() == ".json":
                    config = json.load(f)
                else:
                    raise BatchConfigurationError(
                        f"Formato não suportado: {file_path.suffix}"
                    )
        else:
            # Tenta parsear como JSON direto
            try:
                config = json.loads(batch_cfg)
            except json.JSONDecodeError:
                raise BatchConfigurationError(
                    f"Erro ao parsear configuração: {batch_cfg}"
                )

        # Validação básica
        if not isinstance(config, dict):
            raise BatchConfigurationError("Configuração deve ser um objeto/dicionário")

        # Detectar tipo de task se não especificado
        if "task" not in config:
            # Compatibilidade com formato antigo
            if "experiments" in config:
                config["task"] = {"type": "execution"}
            elif "optimization" in config:
                config["task"] = {"type": "optimization"}
            elif "sensitivity" in config:
                config["task"] = {"type": "sensitivity"}
            else:
                # Formato novo - task é obrigatório
                raise BatchConfigurationError(
                    "Campo 'task' é obrigatório na nova estrutura"
                )

        return config

    def _parse_optimization_config(self, opt_cfg: str) -> Dict[str, Any]:
        """Parseia configuração de otimização."""
        import json
        from pathlib import Path

        import yaml

        # Se não é string, retorna como está (já é dict)
        if not isinstance(opt_cfg, str):
            return opt_cfg

        # Se é caminho de arquivo, carrega
        file_path = Path(opt_cfg)
        if file_path.exists():
            with open(file_path, encoding="utf-8") as f:
                if file_path.suffix.lower() in [".yaml", ".yml"]:
                    config = yaml.safe_load(f)
                elif file_path.suffix.lower() == ".json":
                    config = json.load(f)
                else:
                    raise OptimizationConfigurationError(
                        f"Formato não suportado: {file_path.suffix}"
                    )
        else:
            # Tenta parsear como JSON direto
            try:
                config = json.loads(opt_cfg)
            except json.JSONDecodeError:
                raise OptimizationConfigurationError(
                    f"Erro ao parsear configuração: {opt_cfg}"
                )

        # Validação básica
        if not isinstance(config, dict):
            raise OptimizationConfigurationError(
                "Configuração deve ser um objeto/dicionário"
            )

        return config

    def _parse_sensitivity_config(self, sens_cfg: str) -> Dict[str, Any]:
        """Parseia configuração de sensibilidade."""
        # TODO: Implementar parsing
        if isinstance(sens_cfg, str):
            raise SensitivityConfigurationError(
                "Parsing de arquivo não implementado ainda"
            )
        return sens_cfg

    def _validate_batch_config(self, config: Dict[str, Any]) -> None:
        """Valida configuração de batch."""
        # Detectar se é estrutura nova ou legada
        is_new_structure = (
            "task" in config and "datasets" in config and "algorithms" in config
        )
        is_legacy_structure = "experiments" in config
        is_optimization_structure = (
            "task" in config
            and config.get("task", {}).get("type") == "optimization"
            and "algorithm" in config
            and "dataset" in config
        )
        is_sensitivity_structure = (
            "task" in config
            and config.get("task", {}).get("type") == "sensitivity"
            and "algorithm" in config
            and "dataset" in config
        )

        if is_new_structure:
            # Validar estrutura nova
            required_fields = ["task", "datasets", "algorithms"]
            for field in required_fields:
                if field not in config:
                    raise BatchConfigurationError(f"Campo obrigatório ausente: {field}")

            # Validar task
            task = config["task"]
            if "type" not in task:
                raise BatchConfigurationError("Campo 'type' obrigatório em task")

            task_type = task["type"]

            # Validar estrutura específica por tipo de task
            if task_type == "execution":
                if "execution" not in task:
                    raise BatchConfigurationError(
                        "Campo 'execution' obrigatório para task do tipo 'execution'"
                    )

                execution_config = task["execution"]
                if "executions" not in execution_config:
                    raise BatchConfigurationError(
                        "Campo 'executions' obrigatório em execution"
                    )

                # Validar cada execução
                for i, exec_item in enumerate(execution_config["executions"]):
                    if "dataset" not in exec_item:
                        raise BatchConfigurationError(
                            f"Campo 'dataset' obrigatório na execução {i}"
                        )
                    if "algorithm" not in exec_item:
                        raise BatchConfigurationError(
                            f"Campo 'algorithm' obrigatório na execução {i}"
                        )

            elif task_type == "optimization":
                if "optimization" not in task:
                    raise BatchConfigurationError(
                        "Campo 'optimization' obrigatório para task do tipo 'optimization'"
                    )

            elif task_type == "sensitivity":
                if "sensitivity" not in task:
                    raise BatchConfigurationError(
                        "Campo 'sensitivity' obrigatório para task do tipo 'sensitivity'"
                    )

            # Validar datasets
            if not isinstance(config["datasets"], list):
                raise BatchConfigurationError("Campo 'datasets' deve ser uma lista")

            for dataset in config["datasets"]:
                if "id" not in dataset:
                    raise BatchConfigurationError("Campo 'id' obrigatório em dataset")
                if "type" not in dataset:
                    raise BatchConfigurationError("Campo 'type' obrigatório em dataset")

            # Validar algorithms
            if not isinstance(config["algorithms"], list):
                raise BatchConfigurationError("Campo 'algorithms' deve ser uma lista")

            for algorithm in config["algorithms"]:
                if "id" not in algorithm:
                    raise BatchConfigurationError("Campo 'id' obrigatório em algorithm")
                if "algorithms" not in algorithm:
                    raise BatchConfigurationError(
                        "Campo 'algorithms' obrigatório em algorithm"
                    )

        elif is_optimization_structure:
            # Validar estrutura de otimização
            required_fields = ["task", "algorithm", "dataset", "optimization"]
            for field in required_fields:
                if field not in config:
                    raise BatchConfigurationError(f"Campo obrigatório ausente: {field}")

        elif is_sensitivity_structure:
            # Validar estrutura de sensibilidade
            required_fields = ["task", "algorithm", "dataset", "sensitivity"]
            for field in required_fields:
                if field not in config:
                    raise BatchConfigurationError(f"Campo obrigatório ausente: {field}")

        elif is_legacy_structure:
            # Validar estrutura legada
            required_fields = ["experiments"]
            for field in required_fields:
                if field not in config:
                    raise BatchConfigurationError(f"Campo obrigatório ausente: {field}")

            # Validar cada experimento
            for exp in config["experiments"]:
                if "algorithm" not in exp:
                    raise BatchConfigurationError(
                        "Campo 'algorithm' obrigatório em experimento"
                    )
                if "dataset" not in exp:
                    raise BatchConfigurationError(
                        "Campo 'dataset' obrigatório em experimento"
                    )
        else:
            raise BatchConfigurationError(
                "Estrutura de batch não reconhecida. Use 'task' + 'datasets' + 'algorithms' (nova), 'experiments' (legada), ou estrutura de otimização/sensibilidade"
            )

    def _validate_optimization_config(self, config: Dict[str, Any]) -> None:
        """Valida configuração de otimização."""
        required_fields = ["algorithm", "dataset", "optimization"]
        for field in required_fields:
            if field not in config:
                raise OptimizationConfigurationError(
                    f"Campo obrigatório ausente: {field}"
                )

    def _validate_sensitivity_config(self, config: Dict[str, Any]) -> None:
        """Valida configuração de sensibilidade."""
        required_fields = ["algorithm", "dataset", "parameters"]
        for field in required_fields:
            if field not in config:
                raise SensitivityConfigurationError(
                    f"Campo obrigatório ausente: {field}"
                )

    def _load_dataset(self, dataset_id: str) -> Dataset:
        """Carrega dataset por ID."""
        if not self._dataset_repo.exists(dataset_id):
            raise DatasetNotFoundError(f"Dataset não encontrado: {dataset_id}")

        return self._dataset_repo.load(dataset_id)

    def _resolve_dataset_config(
        self, dataset_id: str, datasets_config: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Resolve configuração de dataset por ID."""
        for dataset in datasets_config:
            if dataset["id"] == dataset_id:
                return dataset
        raise BatchConfigurationError(f"Dataset com ID '{dataset_id}' não encontrado")

    def _resolve_algorithm_config(
        self, algorithm_id: str, algorithms_config: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Resolve configuração de algoritmo por ID."""
        for algorithm in algorithms_config:
            if algorithm["id"] == algorithm_id:
                return algorithm
        raise BatchConfigurationError(
            f"Algoritmo com ID '{algorithm_id}' não encontrado"
        )

    def _resolve_algorithm_configuration(
        self, algorithm_id: str, config: Dict[str, Any]
    ) -> tuple[List[str], Dict[str, Dict[str, Any]]]:
        """
        Resolve configuração de algoritmo por ID.

        Args:
            algorithm_id: ID da configuração do algoritmo
            config: Configuração completa do batch

        Returns:
            tuple: (list_of_algorithm_names, all_algorithm_params)
        """
        algorithms_config = config.get("algorithms", [])

        # Encontrar configuração por ID
        algorithm_config = next(
            (a for a in algorithms_config if a["id"] == algorithm_id), None
        )

        if not algorithm_config:
            raise ValueError(
                f"Configuração de algoritmo com ID '{algorithm_id}' não encontrada"
            )

        # Extrair algoritmos da configuração
        algorithms = algorithm_config.get("algorithms", [])
        if not algorithms:
            raise ValueError(
                f"Nenhum algoritmo definido na configuração '{algorithm_id}'"
            )

        # Retornar todos os algoritmos
        algorithm_names = algorithms

        # Extrair parâmetros específicos de todos os algoritmos
        algorithm_params = algorithm_config.get("algorithm_params", {})

        return algorithm_names, algorithm_params

    def _create_dataset_from_config(self, dataset_config: Dict[str, Any]) -> Dataset:
        """Cria dataset a partir da configuração."""
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

        elif dataset_type == "file":
            filename = params.get("filename")
            if not filename:
                raise BatchConfigurationError(
                    "Campo 'filename' obrigatório para dataset do tipo 'file'"
                )
            return self._dataset_repo.load(filename)
        elif dataset_type == "entrez":
            # Check if Entrez repository is available
            if not self._entrez_repo:
                raise BatchConfigurationError(
                    "Entrez dataset repository not configured. "
                    "Check NCBI_EMAIL environment variable and Biopython installation."
                )

            # Extract Entrez parameters
            query = params.get("query")
            if not query:
                raise BatchConfigurationError(
                    "Field 'query' is required for entrez dataset type"
                )

            db = params.get("db", "nucleotide")
            retmax = params.get("retmax", 20)

            # Fetch dataset from NCBI
            try:
                sequences, metadata = self._entrez_repo.fetch_dataset(
                    query=query, db=db, retmax=retmax, **params
                )

                # Create Dataset object
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

                return dataset

            except Exception as e:
                self._logger.error("Failed to create Entrez dataset: %s", str(e))
                raise BatchConfigurationError(
                    f"Failed to fetch Entrez dataset: {str(e)}"
                ) from e
        else:
            raise BatchConfigurationError(
                f"Tipo de dataset '{dataset_type}' não suportado"
            )

    def _consolidate_batch_results(
        self, results: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Consolida resultados de batch."""
        return {
            "total_experiments": len(results),
            "results": results,
            "summary": {
                "successful": sum(1 for r in results if r.get("status") == "success"),
                "failed": sum(1 for r in results if r.get("status") == "error"),
                "algorithms_used": list(
                    set(
                        r.get("algorithm_name", r.get("algorithm"))
                        for r in results
                        if r.get("algorithm_name") or r.get("algorithm")
                    )
                ),
                "datasets_used": list(
                    set(
                        r.get("dataset_id", r.get("dataset"))
                        for r in results
                        if r.get("dataset_id") or r.get("dataset")
                    )
                ),
            },
        }

    def _export_batch_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """Exporta resultados de batch."""
        # Se results já é uma lista de resultados, usar diretamente
        if isinstance(results, list):
            results_list = results
        elif "results" in results:
            results_list = results["results"]
        else:
            # Se não tem structure esperada, exportar como está
            results_list = [results]

        return self._exporter.export_batch_results(
            results_list, format_type, destination
        )

    def _configure_monitoring(self, monitoring_config) -> None:
        """Configure monitoring service based on configuration."""
        if self._monitoring_service:
            # Configure monitoring interface and update interval
            try:
                # Update monitoring interface if supported
                if hasattr(self._monitoring_service, "set_interface"):
                    self._monitoring_service.set_interface(monitoring_config.interface)

                # Update update interval if supported
                if hasattr(self._monitoring_service, "set_update_interval"):
                    self._monitoring_service.set_update_interval(
                        monitoring_config.update_interval
                    )

                self._logger.info(
                    f"Monitoring configured: interface={monitoring_config.interface}, interval={monitoring_config.update_interval}"
                )
            except Exception as e:
                self._logger.warning(f"Error configuring monitoring: {e}")

    def _export_batch_results_with_config(
        self, results: Dict[str, Any], export_config, task_type: str
    ) -> str:
        """Export batch results using new export configuration."""
        try:
            # Determine formats based on configuration (support multiple formats)
            formats_to_export = []
            if export_config.formats:
                if export_config.formats.get("json", False):
                    formats_to_export.append("json")
                if export_config.formats.get("csv", False):
                    formats_to_export.append("csv")
                if export_config.formats.get("parquet", False):
                    formats_to_export.append("parquet")
                if export_config.formats.get("pickle", False):
                    formats_to_export.append("pickle")

            # Default to JSON if no formats specified
            if not formats_to_export:
                formats_to_export = ["json"]

            # Process destination template
            destination = export_config.destination
            destination = destination.replace("{task_type}", task_type)

            # Handle session replacement differently depending on SessionManager availability
            if hasattr(self, "_session_manager") and self._session_manager:
                # When using SessionManager, FileExporter already has the full session path
                # We should only use relative paths to avoid duplication
                if "{session}" in destination:
                    # For template like "outputs_test/{session}", extract just the relative part
                    # Since FileExporter already has the session base, use empty or just filename
                    destination = ""  # Let FileExporter generate filename automatically
                else:
                    # If no {session} template, use destination as-is (should be relative)
                    pass
            else:
                # Fallback: use current timestamp if no session manager
                from datetime import datetime

                session_fallback = datetime.now().strftime("%Y%m%d_%H%M%S")
                destination = destination.replace("{session}", session_fallback)

            # Filter results based on include configuration
            if export_config.include:
                filtered_results = {}
                for item in export_config.include:
                    if item == "summary" and "batch_summary" in results:
                        filtered_results["summary"] = results["batch_summary"]
                    elif item == "detailed_results" and "detailed_results" in results:
                        filtered_results["detailed_results"] = results[
                            "detailed_results"
                        ]
                    elif item == "plots":
                        # TODO: Add plot data when plots are implemented
                        pass
                    elif item == "logs":
                        # TODO: Add log data when advanced logging is implemented
                        pass
                results = filtered_results

            # Export in all specified formats
            exported_files = []
            for format_type in formats_to_export:
                try:
                    # Try to export with format options if supported
                    try:
                        exported_file = self._exporter.export_results(
                            results,
                            format_type,
                            destination,
                            export_config.format_options,
                        )
                    except TypeError:
                        # Fallback for exporters that don't support options parameter
                        exported_file = self._exporter.export_results(
                            results, format_type, destination
                        )
                    exported_files.append(exported_file)
                except Exception as format_error:
                    self._logger.warning(
                        f"Failed to export in {format_type} format: {format_error}"
                    )
                    # Continue with other formats

            # Return primary exported file (first successful export)
            return exported_files[0] if exported_files else None
        except Exception as e:
            self._logger.error(f"Error exporting results: {e}")
            # Fallback to simple export
            return self._export_batch_results(
                results, "json", f"{task_type}_results.json"
            )

    def _generate_plots_with_config(
        self, results: Dict[str, Any], plots_config, task_type: str
    ) -> None:
        """Generate plots based on configuration."""
        if not plots_config.enabled:
            return

        try:
            self._logger.info("Generating plots...")

            # Debug: show results structure
            self._logger.debug(
                f"Results structure: {list(results.keys()) if isinstance(results, dict) else type(results)}"
            )

            # Get current session path using standard pattern
            # Use session manager to determine base path
            if hasattr(self, "_session_manager") and self._session_manager:
                session_base_path = Path(
                    self._session_manager.get_session_result_path()
                )
                plots_dir = session_base_path / "plots"
            else:
                # Fallback to old pattern using infrastructure settings from .env
                import os
                from datetime import datetime

                # Get infrastructure configuration from environment variables
                session_format = os.getenv("OUTPUT_SESSION_FOLDER_FORMAT", "%Y%m%d_%H%M%S")
                current_session = datetime.now().strftime(session_format)
                base_output_dir = os.getenv("OUTPUT_BASE_DIRECTORY", "outputs")
                session_base_path = Path(base_output_dir) / "results" / current_session
                plots_dir = session_base_path / "plots"

            plots_dir.mkdir(parents=True, exist_ok=True)

            self._logger.info(f"Plots directory created: {plots_dir}")

            # Generate plots based on task type
            if task_type == "execution":
                self._generate_execution_plots(results, plots_config, plots_dir)
            elif task_type == "optimization":
                self._generate_optimization_plots(results, plots_config, plots_dir)
            elif task_type == "sensitivity":
                self._generate_sensitivity_plots(results, plots_config, plots_dir)
            else:
                self._logger.warning(
                    f"Unknown task type for plot generation: {task_type}"
                )

            self._logger.info("Plot generation completed successfully")

        except Exception as e:
            self._logger.error(f"Error generating plots: {e}")

    def _generate_execution_plots(
        self, results: Dict[str, Any], plots_config, plots_dir: Path
    ) -> None:
        """Generate plots for execution tasks."""
        try:
            from src.infrastructure.io.report_generators.execution_report_generator import (
                ExecutionReportGenerator,
            )

            # Create report generator with plots configuration
            config_dict = {
                "plots": {
                    "enabled": plots_config.enabled,
                    "plot_convergence": plots_config.plot_convergence,
                    "plot_comparison": plots_config.plot_comparison,
                    "plot_boxplots": plots_config.plot_boxplots,
                    "plot_scatter": plots_config.plot_scatter,
                    "plot_heatmap": plots_config.plot_heatmap,
                    "plot_runtime": plots_config.plot_runtime,
                    "plot_success_rate": plots_config.plot_success_rate,
                    "formats": plots_config.formats,
                }
            }

            session_path = plots_dir.parent
            report_generator = ExecutionReportGenerator(config_dict, session_path)

            # Prepare batch data for report generator
            batch_data = {
                "batch_results": results.get("results", []),  # Fixed key
                "summary": results.get("summary", {}),
            }

            # Debug log
            self._logger.debug(f"Batch data structure: {batch_data}")

            # Generate plots by calling the report generator
            report_generator.generate_report(batch_data)

            # Copy plots from report/plots to main plots directory for easy access
            report_plots_dir = plots_dir.parent / "report" / "plots"
            if report_plots_dir.exists():
                import shutil

                for plot_file in report_plots_dir.glob("*.*"):
                    shutil.copy2(plot_file, plots_dir)
                self._logger.info(f"Plots copied to: {plots_dir}")

            self._logger.info("Execution plots generated successfully")

        except ImportError as e:
            self._logger.error(f"Cannot import execution report generator: {e}")
        except Exception as e:
            self._logger.error(f"Error generating execution plots: {e}")

    def _generate_optimization_plots(
        self, results: Dict[str, Any], plots_config, plots_dir: Path
    ) -> None:
        """Generate plots for optimization tasks."""
        try:
            from src.infrastructure.orchestrators.optimization_report_generator import (
                OptimizationReportGenerator,
            )

            # Create configuration for optimization report generator
            config_dict = {
                "plots": {
                    "enabled": plots_config.enabled,
                    "plot_convergence": plots_config.plot_convergence,
                    "plot_comparison": plots_config.plot_comparison,
                    "plot_boxplots": plots_config.plot_boxplots,
                    "plot_scatter": plots_config.plot_scatter,
                    "plot_heatmap": plots_config.plot_heatmap,
                    "plot_runtime": plots_config.plot_runtime,
                    "plot_success_rate": plots_config.plot_success_rate,
                    "plot_optimization_history": plots_config.plot_optimization_history,
                    "plot_parameter_importance": plots_config.plot_parameter_importance,
                    "plot_parallel_coordinate": plots_config.plot_parallel_coordinate,
                    "formats": plots_config.formats,
                }
            }

            session_path = plots_dir.parent
            report_generator = OptimizationReportGenerator(config_dict, session_path)

            # Note: This would require optuna study object
            # For now, log that optimization plots need study object
            self._logger.info(
                "Optimization plots require optuna study object (not available in current context)"
            )

        except ImportError as e:
            self._logger.error(f"Cannot import optimization report generator: {e}")
        except Exception as e:
            self._logger.error(f"Error generating optimization plots: {e}")

    def _generate_sensitivity_plots(
        self, results: Dict[str, Any], plots_config, plots_dir: Path
    ) -> None:
        """Generate plots for sensitivity analysis tasks."""
        try:
            from src.infrastructure.io.report_generators.sensitivity_report_generator import (
                SensitivityReportGenerator,
            )

            session_path = plots_dir.parent
            report_generator = SensitivityReportGenerator(session_path)

            # Extract sensitivity analysis results
            analysis_results = results.get("sensitivity_results", [])
            if analysis_results:
                report_dir = session_path / "reports"
                report_dir.mkdir(exist_ok=True)
                report_generator._generate_sensitivity_plots(
                    analysis_results, report_dir
                )
                self._logger.info("Sensitivity plots generated successfully")
            else:
                self._logger.warning(
                    "No sensitivity analysis results found for plot generation"
                )

        except ImportError as e:
            self._logger.error(f"Cannot import sensitivity report generator: {e}")

    def _apply_system_config(
        self, system_config: SystemConfig, parsed_config: BatchConfig
    ) -> None:
        """Apply global system configurations like global_seed substitution."""
        if not system_config:
            return

        # Apply global seed substitution if configured
        if system_config.global_seed is not None:
            self._logger.info(f"Applying global_seed: {system_config.global_seed}")
            self._substitute_global_seed(parsed_config, system_config.global_seed)

        # Apply other system configurations
        # Check force_cleanup from .env (infrastructure setting)
        import os
        force_cleanup = os.getenv("FORCE_CLEANUP", "false").lower() == "true"
        if force_cleanup:
            self._logger.info("Force cleanup enabled (from .env)")

        if system_config.checkpointing:
            self._logger.info(
                f"Checkpointing configured: {system_config.checkpointing}"
            )

        if system_config.error_handling:
            self._logger.info(
                f"Error handling configured: {system_config.error_handling}"
            )

        if system_config.progress_tracking:
            self._logger.info(
                f"Progress tracking configured: {system_config.progress_tracking}"
            )

        if system_config.environment:
            self._logger.info(f"Environment configured: {system_config.environment}")

    def _substitute_global_seed(
        self, parsed_config: BatchConfig, global_seed: int
    ) -> None:
        """Substitute all local seeds with global_seed throughout the configuration."""
        self._logger.debug(f"Substituting seeds with global_seed: {global_seed}")

        # Substitute in algorithm parameters
        for algorithm_config in parsed_config.algorithms:
            if algorithm_config.algorithm_params:
                for alg_name, params in algorithm_config.algorithm_params.items():
                    if isinstance(params, dict) and "seed" in params:
                        old_seed = params.get("seed")
                        params["seed"] = global_seed
                        self._logger.debug(
                            f"Algorithm {alg_name}: seed {old_seed} -> {global_seed}"
                        )

        # Substitute in dataset parameters
        for dataset_config in parsed_config.datasets:
            if dataset_config.parameters and "seed" in dataset_config.parameters:
                old_seed = dataset_config.parameters.get("seed")
                dataset_config.parameters["seed"] = global_seed
                self._logger.debug(
                    f"Dataset {dataset_config.id}: seed {old_seed} -> {global_seed}"
                )

        # Substitute in optimization/sensitivity trial parameters if present
        if hasattr(parsed_config, "optimization") and parsed_config.optimization:
            if (
                parsed_config.optimization.parameters
                and "seed" in parsed_config.optimization.parameters
            ):
                old_seed = parsed_config.optimization.parameters.get("seed")
                parsed_config.optimization.parameters["seed"] = global_seed
                self._logger.debug(f"Optimization: seed {old_seed} -> {global_seed}")

        if hasattr(parsed_config, "sensitivity") and parsed_config.sensitivity:
            if (
                parsed_config.sensitivity.parameters
                and "seed" in parsed_config.sensitivity.parameters
            ):
                old_seed = parsed_config.sensitivity.parameters.get("seed")
                parsed_config.sensitivity.parameters["seed"] = global_seed
                self._logger.debug(f"Sensitivity: seed {old_seed} -> {global_seed}")
