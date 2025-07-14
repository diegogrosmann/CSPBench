"""
Orquestrador de Execução

Concentra toda a lógica de execução de algoritmos CSP,
tanto para execuções únicas quanto batches.
"""

import json
import os
import time
import uuid
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from src.domain import Dataset
from src.domain.errors import AlgorithmExecutionError
from src.infrastructure.logging_config import get_logger
from src.infrastructure.orchestrators.base_orchestrator import BaseOrchestrator


class ExecutionOrchestrator(BaseOrchestrator):
    """Orquestrador responsável pela execução de algoritmos CSP."""

    def __init__(self, algorithm_registry, dataset_repository, monitoring_service=None):
        """
        Inicializa orquestrador de execução.

        Args:
            algorithm_registry: Registry de algoritmos
            dataset_repository: Repositório de datasets
            monitoring_service: Serviço de monitoramento opcional
        """
        super().__init__(monitoring_service)
        self._algorithm_registry = algorithm_registry
        self._dataset_repository = dataset_repository
        self._executions: Dict[str, Dict[str, Any]] = {}
        self._current_batch_config: Optional[Dict[str, Any]] = None
        self._partial_results_file: Optional[str] = None
        self._logger = get_logger(__name__)

    def execute(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Implementa método abstrato do BaseOrchestrator.

        Args:
            config: Configuração da execução

        Returns:
            Dict[str, Any]: Resultado da execução
        """
        # Implementação padrão que delega para execute_batch
        results = self.execute_batch(config, self.monitoring_service)
        return {"results": results}

    def set_batch_config(self, batch_config: Dict[str, Any]) -> None:
        """Define configuração do batch atual."""
        self._current_batch_config = batch_config
        self._logger.debug(f"Configuração de batch definida: {type(batch_config)}")

        # Configurar salvamento parcial se habilitado
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
        Executa um algoritmo único.

        Args:
            algorithm_name: Nome do algoritmo a executar
            dataset: Dataset para processamento
            params: Parâmetros específicos do algoritmo
            timeout: Timeout em segundos
            monitoring_service: Serviço de monitoramento opcional

        Returns:
            Dict[str, Any]: Resultado da execução
        """
        from algorithms import global_registry

        # Verifica se algoritmo existe
        if algorithm_name not in global_registry:
            raise AlgorithmExecutionError(
                f"Algoritmo '{algorithm_name}' não encontrado"
            )

        algorithm_class = global_registry[algorithm_name]
        params = params or {}

        # Aplicar configurações de infraestrutura se batch_config disponível
        if self._current_batch_config:
            infrastructure_config = self._current_batch_config.get("infrastructure", {})
            history_config = infrastructure_config.get("history", {})

            # Injetar parâmetros de histórico se habilitados
            if history_config.get("save_history", False):
                params = params.copy()  # Não modificar o original
                params["save_history"] = True
                params["history_frequency"] = history_config.get("history_frequency", 1)

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

            # Executa algoritmo
            best_string, max_distance, metadata = algorithm.run()
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
        """Executa batch com estrutura nova (datasets + algorithms + executions)."""
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
        dataset_execution_index = 0

        for execution in executions:
            execution_index += 1
            execution_name = execution.get("nome", f"Execução {execution_index}")

            # Atualizar configuração atual
            if monitoring_service:
                self._update_execution_monitoring(
                    execution,
                    execution_index,
                    execution_name,
                    datasets_config,
                    algorithms_config,
                    monitoring_service,
                )

            # Resolver datasets
            dataset_ids = execution["datasets"]
            for dataset_id in dataset_ids:
                dataset_execution_index += 1

                dataset_config = next(
                    (d for d in datasets_config if d["id"] == dataset_id), None
                )
                if not dataset_config:
                    results.append(
                        {
                            "execution_name": execution.get("nome", "unknown"),
                            "dataset_id": dataset_id,
                            "status": "error",
                            "error": f"Dataset com ID '{dataset_id}' não encontrado",
                        }
                    )
                    continue

                # Carregar dataset e executar algoritmos
                dataset_results = self._execute_dataset_algorithms(
                    execution,
                    dataset_config,
                    dataset_id,
                    dataset_repo,
                    algorithms_config,
                    monitoring_service,
                )
                results.extend(dataset_results)

            # Atualizar configurações completadas
            if monitoring_service:
                monitoring_service.update_execution_data(
                    completed_configs=execution_index,
                    completed_executions=execution_index,
                )

        return results

    def _execute_legacy_batch(
        self, batch_config: Dict[str, Any], monitoring_service=None
    ) -> List[Dict[str, Any]]:
        """Executa batch com estrutura legada (experiments)."""
        from src.infrastructure import FileDatasetRepository

        dataset_repo = FileDatasetRepository("./datasets")

        results = []
        for exp in batch_config.get("experiments", []):
            try:
                print(f"[DEBUG] Processando experimento: {exp}")
                # Carrega dataset
                dataset = dataset_repo.load(exp["dataset"])
                print(f"[DEBUG] Dataset carregado: {len(dataset.sequences)} sequências")

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
        """Verifica se deve salvar resultados parciais."""
        if not self._current_batch_config:
            return False

        infrastructure = self._current_batch_config.get("infrastructure", {})
        result_config = infrastructure.get("result", {})
        return result_config.get("save_partial_results", False)

    def _setup_partial_results_file(self) -> None:
        """Configura arquivo para salvamento de resultados parciais."""
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
        """Salva um resultado parcial no arquivo."""
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

        # Atualizar dados iniciais
        monitoring_service.update_execution_data(
            current_execution=(
                executions[0].get("nome", "Execução 1") if executions else "Execução"
            ),
            total_executions=len(executions),
            completed_executions=0,
            total_algorithms=total_algorithms,
            completed_algorithms=0,
            total_datasets=total_dataset_executions,
            current_dataset_index=0,
            current_config_name=(
                executions[0].get("nome", "Execução 1") if executions else "Execução"
            ),
            total_configs=len(executions),
            completed_configs=0,
        )

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

        # Criar descrição da execução
        execution_desc = f"{execution_name} ({len(execution['datasets'])} datasets, {len(execution_algorithms)} algoritmos)"

        monitoring_service.update_execution_data(
            current_execution=execution_desc,
            current_config_name=execution_name,
            completed_configs=execution_index - 1,
            current_task_info=f"Iniciando execução {execution_name}",
        )

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
            # Carregar dataset
            if dataset_config["tipo"] == "file":
                filename = dataset_config["parametros"]["filename"]
                dataset = dataset_repo.load(filename)
            else:
                # Para datasets sintéticos, criar usando gerador
                dataset = self._create_dataset_from_config(dataset_config)

            self._logger.info(
                f"Dataset {dataset_id} carregado: {len(dataset.sequences)} sequências"
            )

            # Obter configurações de algoritmos da execução
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
                            "execution_name": execution.get("nome", "unknown"),
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

                    # Executar repetições
                    for rep in range(repetitions):
                        rep_id = f"{algorithm_name}_{dataset_id}_{rep+1}"

                        try:
                            # Notificar monitoramento de novo item
                            if monitoring_service:
                                from src.presentation.monitoring.interfaces import (
                                    HierarchicalContext,
                                )

                                context = HierarchicalContext(
                                    dataset_id=dataset_id,
                                    algorithm_id=algorithm_name,
                                    repetition_id=f"{rep+1}/{repetitions}",
                                )
                                monitoring_service.update_item(
                                    rep_id, 0.0, "Iniciando", context
                                )

                            # Executar algoritmo
                            result = self.execute_single(
                                algorithm_name, dataset, params, monitoring_service
                            )

                            # Adicionar informações de contexto
                            result.update(
                                {
                                    "execution_name": execution.get("nome", "unknown"),
                                    "dataset_id": dataset_id,
                                    "algorithm_id": algorithm_id,
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
                                "execution_name": execution.get("nome", "unknown"),
                                "dataset_id": dataset_id,
                                "algorithm_id": algorithm_id,
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
                                monitoring_service.finish_item(
                                    rep_id, False, error_result, str(e)
                                )

        except Exception as e:
            self._logger.error(
                f"Erro no carregamento/processamento do dataset {dataset_id}: {e}"
            )
            results.append(
                {
                    "execution_name": execution.get("nome", "unknown"),
                    "dataset_id": dataset_id,
                    "status": "error",
                    "error": f"Erro no dataset: {e}",
                }
            )

        return results

    def _create_dataset_from_config(self, dataset_config: Dict[str, Any]):
        """Cria dataset a partir da configuração."""
        from src.domain.dataset import SyntheticDatasetGenerator

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
        else:
            raise ValueError(
                f"Tipo de dataset '{dataset_type}' não suportado em dataset sintético"
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
