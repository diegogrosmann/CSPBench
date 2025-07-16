"""
Serviço de Experimentos CSPBench

Orquestra a execução de experimentos, otimizações e análises de sensibilidade.
Implementa casos de uso da aplicação sem dependência de infraestrutura.
"""

from typing import Any, Dict, List, Optional

from src.application.ports import (
    AlgorithmRegistry,
    DatasetRepository,
    ExecutorPort,
    ExportPort,
)
from src.application.services.config_parser import (
    ConfigurationParser,
    ConfigurationValidator,
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
from src.presentation.monitoring.interfaces import TaskType


class ExperimentService:
    """
    Serviço principal para execução de experimentos CSP.

    Coordena a execução de diferentes tipos de experimentos através das portas
    definidas, mantendo a lógica de aplicação separada da infraestrutura.
    """

    def __init__(
        self,
        dataset_repo: DatasetRepository,
        exporter: ExportPort,
        executor: ExecutorPort,
        algo_registry: AlgorithmRegistry,
        monitoring_service: Optional[Any] = None,
    ):
        """
        Inicializa o serviço de experimentos.

        Args:
            dataset_repo: Repositório de datasets
            exporter: Port para exportação de resultados
            executor: Port para execução de algoritmos
            algo_registry: Registry de algoritmos
            monitoring_service: Serviço de monitoramento (opcional)
        """
        self._dataset_repo = dataset_repo
        self._exporter = exporter
        self._executor = executor
        self._algo_registry = algo_registry
        self._monitoring_service = monitoring_service
        self._logger = get_logger(__name__)

    def _update_batch_logging(self, batch_config: Dict[str, Any]) -> None:
        """Atualiza configuração de logging baseado no batch específico."""
        try:
            # Verificar configuração de logging no formato atual
            if "logging" in batch_config:
                log_config = batch_config["logging"]

                # Atualizar nível se especificado
                if "level" in log_config:
                    new_level = log_config["level"]
                    LoggerConfig.set_level(new_level)
                    self._logger.info(f"Nível de log atualizado para: {new_level}")

            # Suportar formato legado também
            elif "advanced" in batch_config and "logs" in batch_config["advanced"]:
                log_config = batch_config["advanced"]["logs"]

                # Verificar se logs estão habilitados
                if not log_config.get("enable", True):
                    return

                # Atualizar nível se especificado
                if "log_level" in log_config:
                    new_level = log_config["log_level"]
                    LoggerConfig.set_level(new_level)
                    self._logger.info(f"Nível de log atualizado para: {new_level}")
        except Exception as e:
            self._logger.warning(
                f"Erro ao atualizar configuração de logging do batch: {e}"
            )

    def run_batch(self, batch_cfg: str) -> Dict[str, Any]:
        """
        Executa experimentos em lote a partir de configuração.

        Args:
            batch_cfg: Caminho ou conteúdo da configuração de batch

        Returns:
            Dict[str, Any]: Resultados consolidados do batch

        Raises:
            BatchConfigurationError: Se configuração inválida
            BatchExecutionError: Se erro na execução
            DatasetNotFoundError: Se dataset não encontrado
            AlgorithmNotFoundError: Se algoritmo não encontrado
        """
        try:
            self._logger.info(f"Iniciando execução de batch: {batch_cfg}")

            # Parsear configuração de batch usando novo parser modular
            batch_config = self._parse_batch_config(batch_cfg)

            # Extrair metadados usando novo parser
            metadata = ConfigurationParser.parse_metadata(batch_config)
            self._logger.info(f"Batch carregado: {metadata.nome}")

            # Atualizar configuração de logging se especificada no batch
            self._update_batch_logging(batch_config)

            # Validar estrutura e determinar tipo
            task_type = ConfigurationValidator.validate_batch_structure(batch_config)
            self._logger.info(f"Tipo de task detectado: {task_type}")

            # Processar de acordo com o tipo de task
            if task_type == "execution":
                results = self._process_execution_batch(batch_config)
            elif task_type == "optimization":
                results = self._process_optimization_batch(batch_config)
            elif task_type == "sensitivity":
                results = self._process_sensitivity_batch(batch_config)
            else:
                raise BatchConfigurationError(
                    f"Tipo de task não suportado: {task_type}"
                )

            # Consolidar resultados
            consolidated_results = self._consolidate_batch_results(results["results"])

            total_experiments = len(results["results"])
            successful = consolidated_results.get("summary", {}).get("successful", 0)
            failed = consolidated_results.get("summary", {}).get("failed", 0)
            self._logger.info(
                f"Batch concluído: {total_experiments} experimentos, {successful} sucessos, {failed} falhas"
            )

            # Exportar se configurado
            output_config = batch_config.get("output", {})
            export_config = batch_config.get("export", {})

            # Suportar ambas as configurações: output.save_results e export.enabled
            should_export = output_config.get(
                "save_results", False
            ) or export_config.get("enabled", False)

            if should_export:
                # Determinar nome do arquivo baseado no tipo de task
                default_filename = f"{task_type}_results.json"

                # Usar configurações de export se disponível, senão usar padrões
                format_type = (
                    export_config.get("formats", {}).get("json", True)
                    and "json"
                    or "txt"
                )
                destination = export_config.get("destination", default_filename)

                # Estruturar dados de export baseado no tipo de task
                export_data = {
                    "batch_summary": consolidated_results,
                    "detailed_results": results.get("results", []),
                }

                # Adicionar campos específicos por tipo
                if task_type == "sensitivity":
                    export_data["sensitivity_summaries"] = results.get(
                        "sensitivity_summaries", []
                    )

                export_path = self._export_batch_results(
                    export_data,
                    format_type,
                    destination,
                )
                self._logger.info(f"Resultados exportados para: {export_path}")

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
        Executa otimização de hiperparâmetros.

        Args:
            opt_cfg: Caminho ou conteúdo da configuração de otimização

        Returns:
            Dict[str, Any]: Resultados da otimização

        Raises:
            OptimizationConfigurationError: Se configuração inválida
            OptimizationExecutionError: Se erro na execução
            DatasetNotFoundError: Se dataset não encontrado
            AlgorithmNotFoundError: Se algoritmo não encontrado
        """
        try:
            # Parsear configuração de otimização
            opt_config = self._parse_optimization_config(opt_cfg)

            # Validar configuração
            self._validate_optimization_config(opt_config)

            # Carregar dataset
            dataset = self._load_dataset(opt_config["dataset"])

            # Verificar se algoritmo existe
            algorithm_name = opt_config["algorithm"]
            if not self._algo_registry.algorithm_exists(algorithm_name):
                raise AlgorithmNotFoundError(
                    f"Algoritmo não encontrado: {algorithm_name}"
                )

            # Executar otimização
            results = self._executor.execute_optimization(
                algorithm_name,
                dataset,
                opt_config,
                monitoring_service=self._monitoring_service,
            )

            # Exportar resultados se configurado
            if opt_config.get("export", {}).get("enabled", False):
                export_config = opt_config["export"]
                self._exporter.export_optimization_results(
                    results,
                    export_config.get("destination", "outputs/optimization_results"),
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
        Executa análise de sensibilidade de parâmetros.

        Args:
            sens_cfg: Caminho ou conteúdo da configuração de análise

        Returns:
            Dict[str, Any]: Resultados da análise de sensibilidade

        Raises:
            SensitivityConfigurationError: Se configuração inválida
            SensitivityExecutionError: Se erro na execução
            DatasetNotFoundError: Se dataset não encontrado
            AlgorithmNotFoundError: Se algoritmo não encontrado
        """
        try:
            # Parsear configuração de sensibilidade
            sens_config = self._parse_sensitivity_config(sens_cfg)

            # Validar configuração
            self._validate_sensitivity_config(sens_config)

            # Carregar dataset
            dataset = self._load_dataset(sens_config["dataset"])

            # Verificar se algoritmo existe
            algorithm_name = sens_config["algorithm"]
            if not self._algo_registry.algorithm_exists(algorithm_name):
                raise AlgorithmNotFoundError(
                    f"Algoritmo não encontrado: {algorithm_name}"
                )

            # Executar análise de sensibilidade
            results = self._executor.execute_sensitivity_analysis(
                algorithm_name, dataset, sens_config
            )

            # Exportar resultados se configurado
            if sens_config.get("export", {}).get("enabled", False):
                export_config = sens_config["export"]
                self._exporter.export_results(
                    results,
                    export_config.get("format", "json"),
                    export_config.get("destination", "outputs/sensitivity_results"),
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
        Executa experimento único.

        Args:
            algorithm_name: Nome do algoritmo
            dataset_id: Identificador do dataset
            params: Parâmetros do algoritmo
            timeout: Timeout em segundos

        Returns:
            Dict[str, Any]: Resultados do experimento
        """
        # Verificar se algoritmo existe
        if not self._algo_registry.algorithm_exists(algorithm_name):
            raise AlgorithmNotFoundError(f"Algoritmo não encontrado: {algorithm_name}")

        # Carregar dataset
        dataset = self._load_dataset(dataset_id)

        # Criar configuração de batch para execução única
        batch_config = {
            "task": {"type": "execution"},
            "experiments": [
                {
                    "algorithm": algorithm_name,
                    "dataset": dataset_id,
                    "params": params or {},
                }
            ],
        }

        # Adicionar timeout se fornecido
        if timeout is not None:
            batch_config["experiments"][0]["timeout"] = timeout

        # Executar como batch e retornar primeiro resultado
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
        """Lista algoritmos disponíveis."""
        return self._algo_registry.list_algorithms()

    def list_available_datasets(self) -> List[str]:
        """Lista datasets disponíveis."""
        return self._dataset_repo.list_available()

    def get_algorithm_info(self, algorithm_name: str) -> Dict[str, Any]:
        """Obtém informações de algoritmo."""
        if not self._algo_registry.algorithm_exists(algorithm_name):
            raise AlgorithmNotFoundError(f"Algoritmo não encontrado: {algorithm_name}")

        return self._algo_registry.get_algorithm_metadata(algorithm_name)

    # Métodos privados auxiliares

    def _process_execution_batch(self, batch_config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Processa batch de execução usando o executor existente.

        Args:
            batch_config: Configuração do batch

        Returns:
            Dict[str, Any]: Resultados da execução
        """
        self._logger.info("Processando batch de execução")

        # Iniciar monitoramento se disponível
        if self._monitoring_service:
            batch_name = batch_config.get("metadados", {}).get("nome", "Execução")
            self._monitoring_service.start_monitoring(TaskType.EXECUTION, batch_name)

        try:
            results = self._executor.execute_batch(
                batch_config, self._monitoring_service
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
        self, batch_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Processa batch de otimização com suporte a múltiplas configurações.

        Args:
            batch_config: Configuração do batch

        Returns:
            Dict[str, Any]: Resultados consolidados da otimização
        """
        self._logger.info("Processando batch de otimização")

        # Iniciar monitoramento se disponível
        if self._monitoring_service:
            batch_name = batch_config.get("metadados", {}).get("nome", "Otimização")
            self._monitoring_service.start_monitoring(TaskType.OPTIMIZATION, batch_name)

        try:
            # Parsear configurações de otimização usando novo parser
            optimization_configs = ConfigurationParser.parse_optimization_configs(
                batch_config
            )

            all_results = []
            optimization_summaries = []

            # Calcular totais para monitoramento
            total_optimizations = len(optimization_configs)

            # Calcular total de execuções (datasets × algoritmos por configuração)
            total_executions = 0
            for opt_config in optimization_configs:
                # Resolver algoritmos para calcular total
                algorithm_names, _ = self._resolve_algorithm_configuration(
                    opt_config.target_algorithm, batch_config
                )
                total_executions += len(opt_config.target_datasets) * len(
                    algorithm_names
                )

            total_datasets = sum(
                len(config.target_datasets) for config in optimization_configs
            )

            current_optimization_index = 0
            current_execution_index = 0

            for opt_config in optimization_configs:
                current_optimization_index += 1
                self._logger.info(f"Executando otimização: {opt_config.nome}")

                # Resolver algoritmos e parâmetros base por ID
                algorithm_names, all_algorithm_params = (
                    self._resolve_algorithm_configuration(
                        opt_config.target_algorithm, batch_config
                    )
                )

                # Processar cada algoritmo individualmente
                for algorithm_name in algorithm_names:
                    current_dataset_index = 0

                    # Obter parâmetros base específicos para este algoritmo
                    base_params = all_algorithm_params.get(algorithm_name, {})

                    # Processar cada dataset na configuração
                    for dataset_id in opt_config.target_datasets:
                        current_dataset_index += 1
                        current_execution_index += 1
                        self._logger.info(f"Processando dataset: {dataset_id}")

                        try:
                            # Resolver configuração do dataset e carregá-lo
                            dataset_config = self._resolve_dataset_config(
                                dataset_id, batch_config.get("datasets", [])
                            )

                            # Carregar dataset baseado na configuração
                            if dataset_config["tipo"] == "file":
                                filename = dataset_config["parametros"]["filename"]
                                dataset = self._load_dataset(filename)
                            else:
                                # Para datasets sintéticos, criar usando o método existente
                                dataset = self._create_dataset_from_config(
                                    dataset_config
                                )

                            self._logger.info(
                                f"Dataset {dataset_id} carregado: {len(dataset.sequences)} sequências"
                            )

                            # Extrair configurações de recursos
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

                            # Preparar configuração para o executor
                            executor_config = {
                                "study_name": f"{opt_config.study_name}_{algorithm_name}_{dataset_id}",
                                "direction": opt_config.direction,
                                "n_trials": opt_config.n_trials,
                                "timeout_per_trial": opt_config.timeout_per_trial,
                                "parameters": optimization_params,
                                "base_params": base_params,  # Parâmetros base da configuração
                                "optuna_config": opt_config.optuna_config or {},
                                "resources": resources_config,  # Incluir configurações de recursos
                                "internal_jobs": resources_config.get(
                                    "internal_jobs", 4
                                ),  # Paralelismo interno
                            }

                            # Executar otimização
                            optimization_results = self._executor.execute_optimization(
                                algorithm_name,  # Usar nome do algoritmo resolvido
                                dataset,
                                executor_config,
                                self._monitoring_service,
                                config_index=current_execution_index,
                                total_configs=total_executions,
                                dataset_index=current_dataset_index,
                                total_datasets=len(opt_config.target_datasets),
                                dataset_name=dataset_id,  # Passar nome original do dataset
                            )

                            all_results.append(optimization_results)

                            # Criar sumário para esta otimização
                            optimization_summaries.append(
                                {
                                    "optimization_name": opt_config.nome,
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
                                "optimization_name": opt_config.nome,
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
        self, batch_config: Dict[str, Any]
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
            batch_name = batch_config.get("metadados", {}).get(
                "nome", "Análise de Sensibilidade"
            )
            self._monitoring_service.start_monitoring(TaskType.SENSITIVITY, batch_name)

        try:
            # Parsear configurações de sensibilidade usando novo parser
            sensitivity_configs = ConfigurationParser.parse_sensitivity_configs(
                batch_config
            )

            all_results = []
            sensitivity_summaries = []

            for sens_config in sensitivity_configs:
                self._logger.info(
                    f"Executando análise de sensibilidade: {sens_config.nome}"
                )

                # Processar cada dataset na configuração
                for dataset_id in sens_config.target_datasets:
                    self._logger.info(f"Processando dataset: {dataset_id}")

                    try:
                        # Resolver configuração do dataset e carregá-lo
                        dataset_config = self._resolve_dataset_config(
                            dataset_id, batch_config.get("datasets", [])
                        )

                        # Carregar dataset baseado na configuração
                        if dataset_config["tipo"] == "file":
                            filename = dataset_config["parametros"]["filename"]
                            dataset = self._load_dataset(filename)
                        else:
                            # Para datasets sintéticos, criar usando o método existente
                            dataset = self._create_dataset_from_config(dataset_config)

                        self._logger.info(
                            f"Dataset {dataset_id} carregado: {len(dataset.sequences)} sequências"
                        )

                        # Extrair configurações de recursos
                        resources_config = ConfigurationParser.parse_resources_config(
                            batch_config
                        )

                        # Preparar configuração para o executor
                        # Resolver algoritmo e parâmetros base por ID
                        algorithm_name, base_params = (
                            self._resolve_algorithm_configuration(
                                sens_config.target_algorithm, batch_config
                            )
                        )

                        # Processar parâmetros de sensibilidade - nova estrutura
                        sensitivity_params = {}
                        if algorithm_name in sens_config.parameters:
                            sensitivity_params = sens_config.parameters[algorithm_name]
                        else:
                            # Fallback para estrutura antiga (compatibilidade)
                            sensitivity_params = sens_config.parameters

                        executor_config = {
                            "analysis_method": sens_config.analysis_method,
                            "n_samples": sens_config.n_samples,
                            "repetitions_per_sample": sens_config.repetitions_per_sample,
                            "parameters": sensitivity_params,
                            "base_params": base_params,  # Parâmetros base da configuração
                            "output_metrics": sens_config.output_metrics,
                            "method_config": sens_config.method_config or {},
                            "resources": resources_config,  # Incluir configurações de recursos
                            "internal_jobs": resources_config.get(
                                "internal_jobs", 4
                            ),  # Paralelismo interno
                        }

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
                                "analysis_name": sens_config.nome,
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
                            "analysis_name": sens_config.nome,
                            "algorithm": sens_config.target_algorithm,  # Manter ID original no erro
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
            with open(file_path, "r", encoding="utf-8") as f:
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
            with open(file_path, "r", encoding="utf-8") as f:
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
                if "tipo" not in dataset:
                    raise BatchConfigurationError("Campo 'tipo' obrigatório em dataset")

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
        dataset_type = dataset_config["tipo"]
        params = dataset_config.get("parametros", {})

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
            # TODO: Implementar suporte a Entrez
            raise BatchConfigurationError(
                "Tipo de dataset 'entrez' ainda não implementado"
            )
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
