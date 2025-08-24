"""PipelineRunner com suporte a controle (pause/cancel/resume) via WorkManager."""

from __future__ import annotations
from typing import Any
from application.services.work_service import get_work_service
from infrastructure.persistence.work_state_persistence import WorkStatePersistence
from src.domain.config import (
    CSPBenchConfig,
    ExperimentTask,
    OptimizationTask,
    SensitivityTask,
    TasksGroup,
)
from src.domain.dataset import Dataset

from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.application.services.dataset_service import load_dataset
from src.application.ports.repositories import ExecutionEngine

from .experiment_executor import ExperimentExecutor
from .optimization_executor import OptimizationExecutor
from .sensitivity_executor import SensitivityExecutor

import time

# Create module logger
logger = get_logger("CSPBench.PipelineRunner")


class PipelineRunner:
    def __init__(self, work_id: str, store: WorkStatePersistence):
        self.work_id = work_id
        self.store = store
        self.current_state = self.load_pipeline_state() if work_id else None
        self.wm = get_work_service()
        self.execution_controller: ExecutionController | None = None

        logger.info(f"PipelineRunner inicializado para work_id: {work_id}")
        if self.current_state:
            logger.info(
                f"Estado do pipeline carregado: {self.current_state.get('pipeline_status', 'unknown')}"
            )
        else:
            logger.info("Nenhum estado anterior encontrado - pipeline novo")

    def _check_control(self) -> str:
        """Verifica status de controle do pipeline (running/paused/canceled)."""
        if not self.work_id:
            return "running"
        status = self.wm.get_status(self.work_id)
        current_status = status or "running"
        logger.debug(f"Status de controle verificado: {current_status}")
        return current_status

    def load_pipeline_state(self) -> dict[str, Any] | None:
        """Carrega estado salvo do pipeline."""
        if not self.store or not self.work_id:
            logger.debug("Nenhuma store ou work_id disponível para carregar estado")
            return None

        logger.debug(
            f"Tentando carregar estado do pipeline para work_id: {self.work_id}"
        )
        state = self.store.load_pipeline_state(self.work_id)

        if state:
            logger.info(f"Estado do pipeline carregado com sucesso: {state.keys()}")
        else:
            logger.info("Nenhum estado de pipeline encontrado")

        return state

    def save_current_state(
        self,
        task_idx: int,
        dataset_idx: int,
        preset_idx: int,
        alg_idx: int,
        config: CSPBenchConfig,
        status: str = "running",
    ) -> None:
        """Salva estado atual no banco."""
        if not self.store or not self.work_id:
            return

        self.store.save_pipeline_state(
            work_id=self.work_id,
            task_index=task_idx,
            dataset_index=dataset_idx,
            preset_index=preset_idx,
            algorithm_index=alg_idx,
            status=status,
            config=config,
        )

    def pause_and_exit(
        self,
        task_idx: int,
        dataset_idx: int,
        preset_idx: int,
        alg_idx: int,
        config: CSPBenchConfig,
    ) -> None:
        """Pausa pipeline e encerra processo."""
        self.save_current_state(
            task_idx, dataset_idx, preset_idx, alg_idx, config, "paused"
        )
        return

    def run(self, config: CSPBenchConfig) -> None:
        logger.info("Iniciando execução do pipeline")
        logger.info(
            f"Configuração: {len(config.tasks.items)} tarefas do tipo {config.tasks.type}"
        )

        # Log pipeline start
        if self.store and self.work_id:
            self.store.pipeline_started(self.work_id, config)

        # Create ExecutionController with resources and control function
        self.execution_controller = ExecutionController(
            resources=config.resources, check_control=self._check_control
        )
        logger.info("ExecutionController criado com sucesso")

        try:
            # Verificar se deve retomar de estado pausado
            if (
                self.current_state
                and self.current_state.get("pipeline_status") == "paused"
            ):
                logger.info("Estado pausado detectado - retomando execução")

                # Determinar índices de retomada
                start_task = self.current_state["current_task_index"]
                start_dataset = self.current_state["current_dataset_index"]
                start_preset = self.current_state["current_preset_index"]
                start_algorithm = self.current_state["current_algorithm_index"]

                logger.info(
                    f"Retomando do ponto: task={start_task}, dataset={start_dataset}, preset={start_preset}, algorithm={start_algorithm}"
                )

                # Atualizar status para running
                if self.store and self.work_id:
                    self.store.update_pipeline_status(self.work_id, "running")
                    logger.debug("Status do pipeline atualizado para 'running'")
            else:
                # Iniciar do começo
                start_task = start_dataset = start_preset = start_algorithm = 0
                logger.info("Iniciando pipeline do começo")

            # Executar pipeline unificado
            logger.info("Iniciando execução unificada do pipeline")
            self._execute_pipeline(
                config, start_task, start_dataset, start_preset, start_algorithm
            )
            logger.info("Pipeline executado com sucesso")

        except Exception as e:
            logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
            if self.store and self.work_id:
                self.store.log(self.work_id, "error", f"Pipeline error: {e}", {"error_type": e.__class__.__name__})
            raise
        finally:
            # Cleanup execution controller
            if self.execution_controller:
                logger.debug("Limpando ExecutionController")
                self.execution_controller.cleanup()
                logger.info("Pipeline finalizado e recursos liberados")

    def _execute_pipeline(
        self,
        config: CSPBenchConfig,
        start_task: int = 0,
        start_dataset: int = 0,
        start_preset: int = 0,
        start_algorithm: int = 0,
    ) -> None:
        """Executa pipeline a partir dos índices especificados."""
        logger.info(
            f"Iniciando _execute_pipeline com índices: task={start_task}, dataset={start_dataset}, preset={start_preset}, algorithm={start_algorithm}"
        )

        tasks_group: TasksGroup = config.tasks
        tasks = tasks_group.items

        logger.info(f"Total de tarefas a processar: {len(tasks)}")

        # Loop principal unificado
        for task_idx, task in enumerate(tasks):
            if task_idx < start_task:
                logger.debug(f"Pulando tarefa {task_idx} (já processada)")
                continue  # Pular tasks já processadas

            logger.info(f"Processando tarefa {task_idx}: {task.id} ({task.type})")

            if self.work_id and self._check_control() == "canceled":
                logger.warning(
                    f"Pipeline cancelado antes do início da tarefa {task.id}"
                )
                break

            # Log task start
            if self.store and self.work_id:
                self.store.task_started(self.work_id, task.id, {"type": task.type, "index": task_idx})

            # Verificar pause no início de cada task
            if self.work_id and self._check_control() == "paused":
                logger.info(f"Pipeline pausado no início da tarefa {task.id}")
                self.pause_and_exit(
                    task_idx,
                    0 if task_idx > start_task else start_dataset,
                    0,
                    0,
                    config,
                )
                return

            logger.debug(f"Iniciando execução da tarefa: {task.name}")

            dataset_start = start_dataset if task_idx == start_task else 0
            for dataset_idx, dataset in enumerate(task.datasets):
                if task_idx == start_task and dataset_idx < dataset_start:
                    continue  # Pular datasets já processadas

                # Verificar pause antes de cada dataset
                if self.work_id and self._check_control() == "paused":
                    self.pause_and_exit(task_idx, dataset_idx, 0, 0, config)
                    return

                dataset_obj = self._process_dataset(dataset)

                preset_start = (
                    start_preset
                    if (task_idx == start_task and dataset_idx == start_dataset)
                    else 0
                )
                for preset_idx, preset in enumerate(task.algorithms):
                    if (
                        task_idx == start_task
                        and dataset_idx == start_dataset
                        and preset_idx < preset_start
                    ):
                        continue  # Pular presets já processados

                    # Verificar pause antes de cada preset
                    if self.work_id and self._check_control() == "paused":
                        self.pause_and_exit(
                            task_idx, dataset_idx, preset_idx, 0, config
                        )
                        return

                    alg_start = (
                        start_algorithm
                        if (
                            task_idx == start_task
                            and dataset_idx == start_dataset
                            and preset_idx == start_preset
                        )
                        else 0
                    )
                    for alg_idx, alg in enumerate(preset.items):
                        if (
                            task_idx == start_task
                            and dataset_idx == start_dataset
                            and preset_idx == start_preset
                            and alg_idx < alg_start
                        ):
                            continue  # Pular algoritmos já processados

                        # Verificar pause antes de cada algorithm
                        if self.work_id and self._check_control() == "paused":
                            self.pause_and_exit(
                                task_idx, dataset_idx, preset_idx, alg_idx, config
                            )
                            return

                        if self._check_control() == "canceled":
                            if self.store and self.work_id:
                                self.store.log(self.work_id, "warning", "Pipeline canceled mid-execution", {"task": task.id})
                            break

                        # Salvar estado atual antes da execução
                        self.save_current_state(
                            task_idx, dataset_idx, preset_idx, alg_idx, config
                        )

                        # Executar algoritmo
                        self._execute_algorithm(
                            task, alg, dataset_obj, config, self._check_control
                        )
                    else:
                        continue
                    break
            
            # Log task finish
            if self.store and self.work_id:
                self.store.task_finished(self.work_id, task.id, {"status": "ok"})

        # Pipeline completado - limpar estado
        if self.store and self.work_id:
            self.store.clear_pipeline_state(self.work_id)
            self.store.pipeline_finished(self.work_id, True)

    def _process_dataset(self, dataset) -> Dataset:
        """Processa um dataset (cache ou gera novo)."""
        dataset_id = getattr(dataset, "id", None)

        # Verificar cache primeiro
        if self.store and dataset_id and self.store.has_dataset(dataset_id):
            strings: list[str] = self.store.get_dataset_strings(dataset_id)
            if self.store and self.work_id:
                self.store.log(
                    self.work_id,
                    "info",
                    f"Using cached dataset: {dataset_id} ({len(strings)} sequences)",
                    {"dataset_id": dataset_id, "cached": True},
                )
            # Create Dataset object from cached strings
            dataset_name = getattr(dataset, "name", "cached_dataset")
            dataset_obj = Dataset(name=dataset_name, sequences=strings)
            if hasattr(dataset, "id"):
                dataset_obj.id = dataset.id
            return dataset_obj
        else:
            # Gerar/carregar dataset normalmente
            resolver, parameters = load_dataset(dataset)
            try:
                # Return the Dataset object directly instead of just strings
                dataset_obj = resolver
                if hasattr(dataset, "id") and not hasattr(dataset_obj, "id"):
                    dataset_obj.id = dataset.id
                if hasattr(dataset, "name") and not hasattr(dataset_obj, "name"):
                    dataset_obj.name = dataset.name
            except Exception as e:
                if self.store and self.work_id:
                    self.store.error(self.work_id, f"dataset:{getattr(dataset, 'id', 'unknown')}", e)
                dataset_name = getattr(dataset, "name", "error_dataset")
                dataset_obj = Dataset(name=dataset_name, sequences=[])
                if hasattr(dataset, "id"):
                    dataset_obj.id = dataset.id

            # Persistir dataset completo para cache futuro
            try:
                if self.store and dataset_id and dataset_obj.sequences:
                    self.store.persist_complete_dataset(
                        dataset, dataset_obj.sequences, parameters
                    )
                    if self.store and self.work_id:
                        self.store.log(
                            self.work_id,
                            "info",
                            f"Cached new dataset: {dataset_id} ({len(dataset_obj.sequences)} sequences)",
                            {"dataset_id": dataset_id, "cached": False},
                        )
            except Exception as e:
                if self.store and self.work_id:
                    self.store.log(
                        self.work_id,
                        "warning",
                        f"Failed to cache dataset {dataset_id}: {e}",
                        {"dataset_id": dataset_id, "error": str(e)},
                    )

        return dataset_obj

    def _get_executor(self, task) -> ExecutionEngine:
        """Factory method para obter engine de execução apropriada baseada no tipo da task."""
        if isinstance(task, ExperimentTask):
            return ExperimentExecutor()
        elif isinstance(task, OptimizationTask):
            return OptimizationExecutor()
        elif isinstance(task, SensitivityTask):
            return SensitivityExecutor()
        else:
            raise ValueError(f"Tipo de task não suportado: {type(task)}")

    def _execute_algorithm(self, task, alg, dataset_obj, config, check_control):
        """Execute a single algorithm with proper result capture."""
        try:
            engine = self._get_executor(task)

            # Pass ExecutionController to the execution engine if it has the attribute
            if hasattr(engine, "execution_controller"):
                engine.execution_controller = self.execution_controller

            result = engine.run(
                task=task,
                dataset_obj=dataset_obj,
                alg=alg,
                resources=config.resources,
                work_id=self.work_id,
                system_config=config.system,
                check_control=check_control,
                store=self.store,
            )

            # Store result in database if available
            if result and self.store and self.work_id:
                try:
                    # Results are already stored by the execution engines
                    self.store.log(
                        self.work_id,
                        "info",
                        f"Algorithm execution completed for {alg.name}",
                        {
                            "algorithm": alg.name,
                            "dataset": (
                                dataset_obj.id
                                if hasattr(dataset_obj, "id")
                                else "unknown"
                            ),
                            "status": result.get("status", "unknown"),
                        },
                    )
                except Exception as e:
                    self.store.log(
                        self.work_id,
                        "warning",
                        f"Failed to log result for {alg.name}: {e}",
                        {"algorithm": alg.name, "error": str(e)},
                    )

            return result

        except KeyboardInterrupt:
            raise
        except Exception as e:
            if self.store and self.work_id:
                self.store.log(
                    self.work_id,
                    "error",
                    f"Failed to execute {alg.name}: {str(e)}",
                    {
                        "algorithm": alg.name,
                        "dataset": (
                            dataset_obj.id if hasattr(dataset_obj, "id") else "unknown"
                        ),
                        "error": str(e),
                    },
                )
            return None
