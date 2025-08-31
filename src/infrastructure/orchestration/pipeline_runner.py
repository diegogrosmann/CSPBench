"""PipelineRunner com suporte a controle (pause/cancel/resume) via WorkManager."""

from __future__ import annotations
import time
import os
from typing import Any, Dict
from src.domain.status import BaseStatus
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)
from src.application.services.work_service import get_work_service
from src.infrastructure.persistence.work_state import (
    WorkStatePersistence,
    WorkScopedPersistence,
)
from src.domain.config import (
    CSPBenchConfig,
    ExperimentTask,
    OptimizationTask,
    SensitivityTask,
    TasksGroup,
)
from src.domain.dataset import Dataset

from src.infrastructure.execution_control import ExecutionController, ExecutionLimitError
from src.infrastructure.logging_config import get_logger
from src.application.services.dataset_service import load_dataset
from src.application.ports.repositories import ExecutionEngine

from .experiment_executor import ExperimentExecutor
from .optimization_executor import OptimizationExecutor
from .sensitivity_executor import SensitivityExecutor

# Create module logger
logger = get_logger("CSPBench.PipelineRunner")


class PipelineRunner:
    def __init__(self, work_store: WorkScopedPersistence):
        self.work_store: WorkScopedPersistence = work_store
        self.work_id = work_store.work_id
        # ExecutionController will be initialized in run() with proper resources config
        self.execution_controller: ExecutionController = None
        self.datasets_cache: dict[str, Dataset] = {}
        
        logger.info(f"PipelineRunner inicializado para work_id: {self.work_id}")

    def run(self, config: CSPBenchConfig) -> None:
        """Executa o pipeline seguindo o novo fluxo reorganizado."""
        logger.info("Iniciando execução do pipeline reorganizado")
        logger.info(
            f"Configuração: {len(config.tasks.items)} tarefas do tipo {config.tasks.type}"
        )
        
        # Initialize ExecutionController with resources configuration
        self.execution_controller = ExecutionController(
            work_id=self.work_id, 
            resources=config.resources
        )
        
        # Apply resource limits to main process
        if config.resources:
            if config.resources.memory:
                self.execution_controller.apply_memory_limits()
            if config.resources.cpu:
                self.execution_controller.apply_cpu_limits()

        try:
            self.work_store.update_work_status(BaseStatus.RUNNING)

            # Fase 1: Gerar todos os datasets primeiro
            logger.info("=== FASE 1: Gerando datasets ===")
            self.execution_controller.check_batch_timeout()  # Check timeout before each phase
            self._generate_all_datasets(config)

            # Fase 2: Gerar todas as combinações
            logger.info("=== FASE 2: Gerando combinações ===")
            self.execution_controller.check_batch_timeout()  # Check timeout before each phase
            self._generate_pipeline_combinations(config)

            time.sleep(2)  # Pequena pausa para garantir que as combinações sejam processadas

            # Fase 3: Executar combinações
            logger.info("=== FASE 3: Executando combinações ===")
            self.execution_controller.check_batch_timeout()  # Check timeout before each phase
            status = self._execute_combinations(config)

            if status == BaseStatus.RUNNING:
                logger.error("Pipeline ainda em execução após finalização")
                self.work_store.update_work_status(BaseStatus.FAILED)
            else:
                logger.info("Pipeline executado com sucesso")
                self.work_store.update_work_status(status)

        except ExecutionLimitError as timeout_exc:
            logger.error(f"Pipeline interrompido por timeout: {timeout_exc}")
            if self.work_store:
                self.work_store.work_error(timeout_exc)
                self.work_store.update_work_status(BaseStatus.CANCELED, error=str(timeout_exc))
        except Exception as e:
            logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
            if self.work_store:
                self.work_store.work_error(e)
                self.work_store.update_work_status(BaseStatus.FAILED, error=str(e))
            raise e
        finally:
            # Cleanup execution controller
            if self.execution_controller:
                logger.debug("Limpando ExecutionController")
                self.execution_controller.cleanup()
                logger.info("Pipeline finalizado e recursos liberados")
            # Finalization export (raw db + manifest + full_results + summary)
            try:
                from src.infrastructure.persistence.work_state.core import WorkStatePersistence  # local import to avoid cycles
                from src.infrastructure.export.finalization_service import (
                    FinalizationConfig,
                    FinalizationService,
                )
                # Infer DB path and output path from work_store
                store: WorkStatePersistence = self.work_store._store  # type: ignore[attr-defined]
                db_path = store.db_path  # Path object
                # Attempt to get output directory from work row (output_path field) if present
                output_dir = None
                try:
                    output_dir = store.get_work_output_path(self.work_id)  # may return Path or None
                except Exception:
                    output_dir = None
                if output_dir is None:
                    # fallback: derive from environment OUTPUT_BASE_DIRECTORY
                    from pathlib import Path
                    base = os.environ.get("OUTPUT_BASE_DIRECTORY", "./data/outputs")
                    output_dir = Path(base) / self.work_id
                output_dir.mkdir(parents=True, exist_ok=True)
                cfg = FinalizationConfig(
                    work_id=self.work_id,
                    db_path=db_path,
                    output_dir=output_dir,
                    tool_version=FinalizationService.detect_tool_version(),
                )
                FinalizationService(cfg).run()
            except Exception as fe:  # pragma: no cover - defensive
                logger.error(f"Erro durante finalização de exportação: {fe}")

    def _generate_all_datasets(self, config: CSPBenchConfig) -> None:
        """Fase 1: Gera todos os datasets necessários."""
        # Coletar todos os dataset IDs únicos utilizados nas tasks
        used_dataset_ids = set()

        for task in config.tasks.items:
            # task.datasets agora é List[str] contendo IDs de datasets
            used_dataset_ids.update(task.datasets)

        logger.info(
            f"Encontrados {len(used_dataset_ids)} dataset IDs únicos utilizados: {sorted(used_dataset_ids)}"
        )

        # Verificar se todos os IDs referenciados existem na configuração
        missing_datasets = used_dataset_ids - set(config.datasets.keys())
        if missing_datasets:
            raise ValueError(
                f"Dataset IDs referenciados mas não definidos: {sorted(missing_datasets)}"
            )

        # Processar cada dataset utilizado
        for dataset_id in sorted(used_dataset_ids):

            logger.info(f"Processando dataset: {dataset_id}")

            try:
                # Obter a configuração do dataset do dicionário config.datasets
                dataset_config = config.datasets[dataset_id]
                self._process_dataset(dataset_config)

                logger.info(
                    f"Dataset {dataset_id} processado."
                )

            except Exception as e:
                logger.error(f"Erro ao processar dataset {dataset_id}: {e}")
                raise

    def _generate_pipeline_combinations(self, config: CSPBenchConfig) -> None:
        """Fase 2: Gera todas as combinações possíveis do pipeline."""
        try:
            requeued = self.work_store.init_combination()
            if requeued:
                logger.info(
                    "Combinações em andamento/pausadas/canceladas reiniciadas para 'queued'. Verificando novas combinações..."
                )
                return  # Não gerar novas combinações se houver em andamento

            # Sempre gerar combinações: tanto para casos novos quanto para adicionar novas após reset
            combinations = []

            for task in config.tasks.items:
                total_sequences = None

                if isinstance(task, ExperimentTask):
                    total_sequences = task.repetitions
                elif isinstance(task, OptimizationTask):
                    total_sequences = task.config.get("trials", 50) if task.config else 50
                elif isinstance(task, SensitivityTask):
                    total_sequences = task.config.get("samples", 100) if task.config else 100

                # task.datasets agora é List[str] contendo IDs de datasets
                for dataset_id in task.datasets:

                    # task.algorithms agora é List[str] contendo IDs de algorithm presets
                    for preset_id in task.algorithms:
                        # Obter o preset do dicionário config.algorithms
                        preset = config.algorithms.get(preset_id)
                        if not preset:
                            self.work_store.preset_error(
                                preset_id, "Algorithm não encontrado na configuração"
                            )
                            logger.warning(
                                f"Algorithm preset '{preset_id}' não encontrado na configuração"
                            )
                            continue

                        for alg in preset.items:
                            combination = {
                                "task_id": task.id,
                                "dataset_id": dataset_id,
                                "preset_id": preset_id,
                                "algorithm_id": alg.name,
                                "mode": task.type,
                                "total_sequences": total_sequences,
                            }
                            combinations.append(combination)

            logger.info(f"Geradas {len(combinations)} combinações para processamento")

            # Submeter todas as combinações (INSERT OR IGNORE garante que duplicatas sejam ignoradas)
            if combinations:
                inserted_count = self.work_store.submit_combinations(combinations)
                logger.info(f"{inserted_count} novas combinações inseridas no banco")
            else:
                logger.warning("Nenhuma combinação válida foi gerada")
        except Exception as e:
            self.work_store.combination_error("N/A", e)
            raise e

    def _execute_combinations(self, config: CSPBenchConfig) -> BaseStatus:
        """Fase 3: Executa todas as combinações pendentes."""
        if not self.work_store:
            logger.error("Work store não disponível para execução")
            return BaseStatus.FAILED

        execution_statuses = []  # Lista para armazenar os status das execuções

        while True:
            # Se controle indicar que não está mais em RUNNING, parar
            if self.execution_controller.check_status() != BaseStatus.RUNNING:
                return self.execution_controller.check_status()

            # Obter próxima combinação pendente
            combination = self.work_store.get_next_pending_combination()
            logger.debug(
                "[PipelineRunner] Próxima combinação retornada=%s",
                combination,
            )
            if not combination:
                logger.info("Todas as combinações foram processadas")
                break

            work_combination = CombinationScopedPersistence(
                self.work_store.store,  # store base
                combination["id"],
            )

            # Executar combinação e armazenar o status
            status = self._execute_single_combination(
                combination, config, work_combination
            )
            execution_statuses.append(status)

        # Verificar a lista de status e determinar o resultado final
        if not execution_statuses:
            # Nenhuma combinação foi executada
            return BaseStatus.COMPLETED

        # Se existe falha, retorna falha
        if BaseStatus.FAILED in execution_statuses:
            return BaseStatus.FAILED

        # Se existe erro, retorna erro
        if BaseStatus.ERROR in execution_statuses:
            return BaseStatus.ERROR

        # Se todos completos, retorna completo
        if all(status == BaseStatus.COMPLETED for status in execution_statuses):
            return BaseStatus.COMPLETED

        # Se outra coisa, retorna erro
        return BaseStatus.ERROR

    def _execute_single_combination(
        self,
        combination: dict[str, Any],
        config: CSPBenchConfig,
        work_combination: CombinationScopedPersistence,
    ) -> BaseStatus:
        """Executa uma única combinação."""
        task_id = combination["task_id"]
        dataset_id = combination["dataset_id"]
        preset_id = combination["preset_id"]
        algorithm_id = combination["algorithm_id"]

        logger.info(
            f"Executando combinação: {task_id}/{dataset_id}/{preset_id}/{algorithm_id}"
        )

        try:
            # Marcar como running
            work_combination.update_combination_status(BaseStatus.RUNNING)

            # Encontrar objetos de configuração
            task = self._find_task(config, task_id)
            if not task:
                raise ValueError(f"Task não encontrada: {task_id}")

            alg = self._find_algorithm(config, preset_id, algorithm_id)
            if not alg:
                raise ValueError(f"Algoritmo não encontrado: {algorithm_id}")

            dataset_obj = self.datasets_cache.get(dataset_id)
            if not dataset_obj:
                raise ValueError(f"Dataset não encontrado no cache: {dataset_id}")

            # Executar algoritmo
            status = self._execute_algorithm(
                work_combination, task, alg, dataset_obj, config
            )

            work_combination.update_combination_status(status)

            logger.info(
                f"Combinação {task_id}/{dataset_id}/{preset_id}/{algorithm_id} concluída: {status}"
            )

            return status

        except Exception as e:
            logger.error(f"Erro na execução da combinação: {e}")
            work_combination.update_combination_status(BaseStatus.FAILED)
            work_combination.record_error(e)
            if self.work_store:
                self.work_store.algorithm_error(algorithm_id, e)
            return BaseStatus.FAILED

    def _execute_algorithm(
        self, work_combination, task, alg, dataset_obj, config
    ) -> BaseStatus:
        """Execute a single algorithm with proper result capture."""
        try:

            engine: ExecutionEngine = self._get_executor(work_combination, config, task)

            result_status = engine.run(
                task=task,
                dataset_obj=dataset_obj,
                alg=alg,
            )

            return result_status

        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise e

    def _process_dataset(self, dataset_config) -> None:
        """Processa um dataset (cache ou gera novo)."""
        dataset_id = getattr(dataset_config, "id", None)

        # Verificar cache primeiro
        if (
            self.work_store
            and dataset_id
            and self.work_store.store.has_dataset(dataset_id)
        ):
            try:
                dataset_obj = self.work_store.get_dataset(dataset_id)
                if dataset_obj is None or not dataset_obj.sequences:
                    logger.warning(
                        f"Dataset {dataset_id} encontrado no cache, mas está vazio"
                    )
                else:
                    self.datasets_cache[dataset_id] = dataset_obj
                    return
            except Exception as e:
                logger.warning(f"Erro ao carregar dataset do cache {dataset_id}: {e}")

        # Sem cache ou dataset vazio, gerar novo
        try:
            # Gerar/carregar dataset normalmente
            logger.info(f"Gerando dataset: {dataset_id}")
            resolver, params = load_dataset(dataset_config)

            meta = {
                "batch_params": dataset_config if dataset_config else {},
                "dataset_params": params if params else {},
                "dataset_statistics": resolver.get_statistics() if resolver else {}
            }

            self.work_store.submit_dataset(id= dataset_id, dataset_obj=resolver, meta=meta)
            self.datasets_cache[dataset_id] = resolver

        except Exception as e:
            logger.error(f"Erro ao gerar/salvar dataset {dataset_id}: {e}")
            if self.work_store:
                self.work_store.dataset_error(dataset_id, e)
            raise ValueError(f"Erro ao gerar/salvar dataset {dataset_id}: {e}")

    def _find_task(self, config: CSPBenchConfig, task_id: str):
        """Encontra uma task pelo ID."""
        for task in config.tasks.items:
            if task.id == task_id:
                return task
        return None

    def _find_algorithm(
        self, config: CSPBenchConfig, preset_id: str, algorithm_id: str
    ):
        """Encontra um algoritmo dentro de um preset usando a nova estrutura."""
        preset = config.algorithms.get(preset_id)
        if not preset:
            return None

        for alg in preset.items:
            if alg.name == algorithm_id:
                return alg
        return None

    def _get_executor(self, work_combination, batch_config, task) -> ExecutionEngine:
        """Factory method para obter engine de execução apropriada baseada no tipo da task."""
        if isinstance(task, ExperimentTask):
            return ExperimentExecutor(
                combination_store=work_combination,
                execution_controller=self.execution_controller,
                batch_config=batch_config,
            )
        elif isinstance(task, OptimizationTask):
            return OptimizationExecutor(
                combination_store=work_combination,
                execution_controller=self.execution_controller,
                batch_config=batch_config,
            )
        elif isinstance(task, SensitivityTask):
            return SensitivityExecutor(
                combination_store=work_combination,
                execution_controller=self.execution_controller,
                batch_config=batch_config,
            )
        else:
            raise ValueError(f"Tipo de task não suportado: {type(task)}")
