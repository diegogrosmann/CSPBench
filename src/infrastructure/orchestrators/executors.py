"""
Executors de Algoritmos

Implementa estratégias de execução para algoritmos CSP.
"""

import json
import os
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.domain import Dataset
from src.domain.errors import AlgorithmExecutionError, AlgorithmNotFoundError


class SequentialExecutor:
    """Executor sequencial simples para algoritmos CSP."""

    def __init__(self):
        """Inicializa executor sequencial."""
        self._executions: Dict[str, Dict[str, Any]] = {}
        self._current_batch_config: Optional[Dict[str, Any]] = None
        self._partial_results_file: Optional[str] = None

    def set_batch_config(self, batch_config: Dict[str, Any]) -> None:
        """Define configuração do batch atual para aplicar configurações de infraestrutura."""
        self._current_batch_config = batch_config

        # Configurar salvamento parcial se habilitado
        if self._should_save_partial_results():
            self._setup_partial_results_file()

    def _should_save_partial_results(self) -> bool:
        """Verifica se deve salvar resultados parciais."""
        if not self._current_batch_config:
            return False

        # Verificar configuração no batch
        infrastructure = self._current_batch_config.get("infrastructure", {})
        result_config = infrastructure.get("result", {})

        return result_config.get("save_partial_results", False)

    def _setup_partial_results_file(self) -> None:
        """Configura arquivo para salvamento de resultados parciais."""

        # Criar diretório de resultados se não existe
        from src.infrastructure import SessionManager

        try:
            # Criar SessionManager e forçar criação de nova sessão
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

    def execute_single(
        self,
        algorithm_name: str,
        dataset: Dataset,
        params: Optional[Dict[str, Any]] = None,
        timeout: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Executa um único algoritmo.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para execução
            params: Parâmetros do algoritmo
            timeout: Timeout em segundos

        Returns:
            Dict: Resultado da execução
        """
        # Garante que os algoritmos estão carregados
        import algorithms  # Force algorithm loading
        from src.domain import global_registry

        # Verifica se algoritmo existe
        if algorithm_name not in global_registry:
            raise AlgorithmNotFoundError(f"Algoritmo '{algorithm_name}' não encontrado")

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

            result_string, max_distance, metadata = algorithm.run()
            end_time = time.time()

            # Calcula estatísticas
            execution_time = end_time - start_time

            result = {
                "execution_id": execution_id,
                "algorithm": algorithm_name,
                "result_string": result_string,
                "max_distance": max_distance,
                "execution_time": execution_time,
                "params": params,
                "metadata": metadata,
                "dataset_info": {
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

    def execute_batch(self, batch_config: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Executa batch de experimentos.

        Args:
            batch_config: Configuração do batch

        Returns:
            List: Resultados das execuções
        """
        # Configurar batch config para aplicar configurações de infraestrutura
        self.set_batch_config(batch_config)

        results = []

        # Detectar estrutura do batch
        task_type = batch_config.get("task", {}).get("type", "execution")

        if task_type == "execution" and "task" in batch_config:
            # Nova estrutura com datasets e algoritmos por ID
            datasets_config = batch_config.get("datasets", [])
            algorithms_config = batch_config.get("algorithms", [])
            executions = batch_config["task"]["execution"]["executions"]

            from src.infrastructure import FileDatasetRepository

            dataset_repo = FileDatasetRepository("./datasets")

            for execution in executions:
                # Resolver dataset
                dataset_id = execution["dataset"]
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

                # Resolver algoritmo
                algorithm_id = execution["algorithm"]
                algorithm_config = next(
                    (a for a in algorithms_config if a["id"] == algorithm_id), None
                )
                if not algorithm_config:
                    results.append(
                        {
                            "execution_name": execution.get("nome", "unknown"),
                            "algorithm_id": algorithm_id,
                            "status": "error",
                            "error": f"Algoritmo com ID '{algorithm_id}' não encontrado",
                        }
                    )
                    continue

                try:
                    # Criar dataset
                    dataset = self._create_dataset_from_config(
                        dataset_config, dataset_repo
                    )

                    # Configurações da execução
                    runs_per_algorithm_per_base = execution.get(
                        "runs_per_algorithm_per_base", 1
                    )
                    num_bases = execution.get("num_bases", 1)
                    timeout = execution.get("timeout", None)

                    # Para cada algoritmo na configuração
                    for algorithm_name in algorithm_config["algorithms"]:
                        # Obter parâmetros específicos do algoritmo
                        algorithm_params = algorithm_config.get(
                            "algorithm_params", {}
                        ).get(algorithm_name, {})

                        # Executar múltiplas vezes
                        for base_idx in range(num_bases):
                            for run_idx in range(runs_per_algorithm_per_base):
                                try:
                                    result = self.execute_single(
                                        algorithm_name,
                                        dataset,
                                        algorithm_params,
                                        timeout,
                                    )

                                    result.update(
                                        {
                                            "execution_name": execution.get(
                                                "nome", "unknown"
                                            ),
                                            "dataset_id": dataset_id,
                                            "algorithm_id": algorithm_id,
                                            "algorithm_name": algorithm_name,
                                            "base_index": base_idx,
                                            "run_index": run_idx,
                                            "status": "success",
                                        }
                                    )

                                except Exception as e:
                                    result = {
                                        "execution_name": execution.get(
                                            "nome", "unknown"
                                        ),
                                        "dataset_id": dataset_id,
                                        "algorithm_id": algorithm_id,
                                        "algorithm_name": algorithm_name,
                                        "base_index": base_idx,
                                        "run_index": run_idx,
                                        "status": "error",
                                        "error": str(e),
                                    }

                                results.append(result)

                                # Salvar resultado parcial se habilitado
                                if self._should_save_partial_results():
                                    self._save_partial_result(result)

                except Exception as e:
                    results.append(
                        {
                            "execution_name": execution.get("nome", "unknown"),
                            "dataset_id": dataset_id,
                            "algorithm_id": algorithm_id,
                            "status": "error",
                            "error": str(e),
                        }
                    )

        elif "experiments" in batch_config:
            # Estrutura legada 'experiments'
            from src.infrastructure import FileDatasetRepository

            dataset_repo = FileDatasetRepository("./datasets")

            for exp in batch_config.get("experiments", []):
                try:
                    print(f"[DEBUG] Processando experimento: {exp}")
                    # Carrega dataset
                    dataset = dataset_repo.load(exp["dataset"])
                    print(
                        f"[DEBUG] Dataset carregado: {len(dataset.sequences)} sequências"
                    )

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
                    # Adiciona resultado de erro
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
        else:
            # Suporte para estrutura antiga 'algorithms' + 'datasets'
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
                        # Adiciona resultado de erro
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

    def _create_dataset_from_config(
        self, dataset_config: Dict[str, Any], dataset_repo
    ) -> Dataset:
        """Cria dataset a partir da configuração."""
        dataset_type = dataset_config["tipo"]
        params = dataset_config.get("parametros", {})

        if dataset_type == "synthetic":
            from src.domain import SyntheticDatasetGenerator

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
                raise ValueError(
                    "Campo 'filename' obrigatório para dataset do tipo 'file'"
                )
            # Remover "saved_datasets/" se estiver presente, pois o dataset_repo já aponta para essa pasta
            if filename.startswith("saved_datasets/"):
                filename = filename[len("saved_datasets/") :]
            return dataset_repo.load(filename)
        elif dataset_type == "entrez":
            # TODO: Implementar suporte a Entrez
            raise NotImplementedError("Tipo de dataset 'entrez' ainda não implementado")
        else:
            raise ValueError(f"Tipo de dataset '{dataset_type}' não suportado")

    def execute_optimization(
        self, algorithm_name: str, dataset: Dataset, optimization_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Executa otimização de hiperparâmetros.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para otimização
            optimization_config: Configuração da otimização

        Returns:
            Dict[str, Any]: Resultados da otimização
        """
        try:
            from src.infrastructure import (
                DomainAlgorithmRegistry,
                FileDatasetRepository,
            )
            from src.infrastructure.orchestrators.optimization_orchestrator import (
                OptimizationOrchestrator,
            )

            # Configurar registry e repository
            algorithm_registry = DomainAlgorithmRegistry()
            dataset_repository = FileDatasetRepository("./datasets")

            # Salvar o dataset temporariamente no repositório para o orquestrador
            dataset_id = f"temp_optimization_{int(time.time())}"
            dataset_repository.save(dataset, dataset_id)

            # Criar configuração completa para o orquestrador
            full_config = {
                "algorithm": algorithm_name,
                "dataset": dataset_id,
                "optimization": optimization_config,
                "export": optimization_config.get("export", {"enabled": True}),
                "monitoring": optimization_config.get("monitoring", {}),
                "resources": optimization_config.get("resources", {}),
                "plots": optimization_config.get("plots", {}),
            }

            # Criar orquestrador
            orchestrator = OptimizationOrchestrator(
                algorithm_registry=algorithm_registry,
                dataset_repository=dataset_repository,
                config=full_config,
            )

            # Executar otimização
            results = orchestrator.run_optimization()

            # Limpar dataset temporário
            try:
                dataset_repository.delete(dataset_id)
            except:
                pass

            return results

        except Exception as e:
            raise AlgorithmExecutionError(f"Erro na otimização: {e}") from e

    def get_execution_status(self, execution_id: str) -> str:
        """
        Obtém status de execução.

        Args:
            execution_id: ID da execução

        Returns:
            str: Status da execução
        """
        execution = self._executions.get(execution_id)
        return execution["status"] if execution else "not_found"

    def cancel_execution(self, execution_id: str) -> bool:
        """
        Cancela execução (não implementado para executor sequencial).

        Args:
            execution_id: ID da execução

        Returns:
            bool: Sempre False (não suportado)
        """
        return False

    def _generate_param_combinations(
        self, param_space: Dict[str, List]
    ) -> List[Dict[str, Any]]:
        """
        Gera combinações de parâmetros para otimização.

        Args:
            param_space: Espaço de parâmetros

        Returns:
            List: Combinações de parâmetros
        """
        if not param_space:
            return [{}]

        # Implementação simples - produto cartesiano
        combinations = [{}]

        for param_name, values in param_space.items():
            new_combinations = []
            for combination in combinations:
                for value in values:
                    new_combo = combination.copy()
                    new_combo[param_name] = value
                    new_combinations.append(new_combo)
            combinations = new_combinations

        return combinations

    def execute_sensitivity_analysis(
        self, algorithm_name: str, dataset: Dataset, sensitivity_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Executa análise de sensibilidade.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para análise
            sensitivity_config: Configuração da análise

        Returns:
            Dict[str, Any]: Resultados da análise de sensibilidade
        """
        try:
            from src.infrastructure import DomainAlgorithmRegistry
            from src.infrastructure.orchestrators.sensitivity_orchestrator import (
                SensitivityOrchestrator,
            )

            # Configurar registry
            algorithm_registry = DomainAlgorithmRegistry()

            # Criar orquestrador de sensibilidade
            orchestrator = SensitivityOrchestrator(algorithm_registry, self)

            # Executar análise
            results = orchestrator.execute_sensitivity_analysis(
                algorithm_name, dataset, sensitivity_config
            )

            return results

        except Exception as e:
            raise AlgorithmExecutionError(
                f"Erro na análise de sensibilidade: {e}"
            ) from e
