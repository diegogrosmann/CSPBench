"""ExperimentExecutor: executa ExperimentTask.

Paraleliza repetições via ProcessPool se ``max_workers > 1`` usando ``ProcessPoolExecutor``.
Implementa controles de recursos, resume functionality e salvamento completo de dados.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from src.application.ports.repositories import AbstractExecutionEngine
from src.domain.config import AlgParams, CSPBenchConfig, ExperimentTask
from src.domain.dataset import Dataset
from src.domain.distance import create_distance_calculator
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.infrastructure.monitoring.persistence_monitor import PersistenceMonitor
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
    ExecutionScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
    WorkScopedPersistence,
)

from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.ExperimentExecutor")


import random
import time

def _worker_exec(
    experiment_unit_id: str,
    seed: int | None,
    strings: list[str],
    alphabet: str,
    distance_method: str,
    use_cache: bool,
    params: dict[str, Any],
    work_id: str,  # Alterado de db_path para work_id
    internal_jobs: int,
    algorithm_name: str,
    cpu_config: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], str]:
    """Função executada em subprocesso (precisa ser picklable)."""
    import psutil

    from src.domain.status import BaseStatus as _BaseStatus
    from src.infrastructure.execution_control import ExecutionController
    from src.infrastructure.logging_config import get_logger
    from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
        WorkScopedPersistence,
    )

    logger = get_logger("CSPBench.WorkerProcess")

    logger.debug(f"[WORKER][{experiment_unit_id}] Iniciando worker_exec.")

    # Timer aleatório para evitar que vários processos iniciem juntos
    delay = random.uniform(0, 3)
    logger.debug(f"[WORKER][{experiment_unit_id}] Aguardando {delay:.2f}s antes de iniciar.")
    time.sleep(delay)

    # Apply CPU configuration to this worker process
    if cpu_config:
        try:
            current_process = psutil.Process()

            # Apply CPU affinity if specified
            if cpu_config.get("affinity"):
                current_process.cpu_affinity(cpu_config["affinity"])
                logger.info(f"[WORKER][{experiment_unit_id}] CPU affinity set to: {cpu_config['affinity']}")

            # Apply CPU priority if max_workers is limited
            if cpu_config.get("exclusive_cores", False):
                current_process.nice(5)  # Slightly lower priority than main process
                logger.info(f"[WORKER][{experiment_unit_id}] CPU priority adjusted for exclusive cores")

        except Exception as e:
            logger.warning(f"[WORKER][{experiment_unit_id}] Cannot apply CPU configuration: {e}")

    logger.warning(f"[WORKER][{experiment_unit_id}] Aqui 1")

    try:
        # Usar factory method para criar work_scoped diretamente
        work_scoped = WorkScopedPersistence.submit(work_id)
        execution_store = work_scoped.for_execution(experiment_unit_id)
        logger.warning(f"[WORKER][{experiment_unit_id}] Aqui 1.2")
    except Exception as e:
        logger.error(f"[WORKER][{experiment_unit_id}] Erro ao recriar store: {e}")
        raise

    logger.warning(f"[WORKER][{experiment_unit_id}] Aqui 2")

    # Controller with work_id for status checks
    dummy_controller = ExecutionController(work_id=work_id)
    dummy_controller._internal_jobs = internal_jobs  # type: ignore[attr-defined]

    logger.warning(f"[WORKER][{experiment_unit_id}] Aqui 3")

    # Apply CPU configuration manually if provided (since we can't pass ResourcesConfig to subprocess)
    if cpu_config:

        logger.warning(f"[WORKER][{experiment_unit_id}] Aqui 4")
        dummy_controller._exclusive_cores = cpu_config.get("exclusive_cores", False)
        if "max_workers" in cpu_config:
            dummy_controller._max_workers = cpu_config["max_workers"]
            dummy_controller._current_workers = cpu_config["max_workers"]

        logger.warning(f"[WORKER][{experiment_unit_id}] Aqui 5")
        # Try to apply CPU affinity in worker process (best effort)
        try:
            if cpu_config.get("exclusive_cores", False) and hasattr(
                dummy_controller, "apply_memory_limits"
            ):
                dummy_controller.apply_memory_limits()
        except Exception as cpu_exc:
            logger.warning(
                f"[WORKER][{experiment_unit_id}] CPU configuration parcialmente aplicada: {cpu_exc}"
            )

    logger.warning(f"[WORKER][{experiment_unit_id}] Aqui 6")
    # Create monitor with controller for cancellation checks
    monitor = PersistenceMonitor(execution_store, execution_controller=dummy_controller)

    distance_calc = create_distance_calculator(
        distance_method=distance_method, strings=strings, use_cache=use_cache
    )
    try:
        execution_store.update_execution_status(_BaseStatus.RUNNING)
        monitor.on_progress(0.0, f"Iniciando instância do algoritmo")

        result = run_algorithm(
            algorithm_name=algorithm_name,
            strings=strings,
            alphabet=alphabet,
            distance_calculator=distance_calc,
            execution_controller=dummy_controller,
            monitor=monitor,
            seed=seed,
            params=params,
        )

        monitor.on_progress(1.0, f"[{experiment_unit_id}] Instância do algoritmo finalizada")
        status_value = result.get("status")
        if hasattr(status_value, "value"):
            status_value = status_value.value
        objective = None
        alg_result = result.get("algorithm_result", None)
        if alg_result is not None:
            max_distance = alg_result.get("max_distance", None)
            if max_distance is not None:
                objective = max_distance

        execution_store.update_execution_status(
            status=status_value,
            result=alg_result,
            params=result.get("actual_params", None),
            objective=objective,
        )
        logger.debug(f"[WORKER][{experiment_unit_id}] worker_exec finalizado com status: {status_value}")
        return result, experiment_unit_id
    except Exception as exc:  # pragma: no cover - defensive
        execution_store.update_execution_status("failed", result={"error": str(exc)})
        logger.error(f"[WORKER][{experiment_unit_id}] worker_exec falhou: {exc}")
        return {
            "status": _BaseStatus.FAILED,
            "error": str(exc),
        }, experiment_unit_id


class ExperimentExecutor(AbstractExecutionEngine):
    _combination_store: CombinationScopedPersistence = None
    _execution_controller: ExecutionController = None
    _batch_config: CSPBenchConfig = None

    def run(
        self,
        task: ExperimentTask,
        dataset_obj: Dataset,
        alg: AlgParams,
    ) -> BaseStatus:
        """Executa o experimento (simplificado)."""

        repetitions = max(1, int(getattr(task, "repetitions", 1)))
        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        dataset_id = getattr(dataset_obj, "id", None)

        global_seed = getattr(self._batch_config.system, "global_seed", None)
        if global_seed is None:
            global_seed = alg.params.get("seed", None)

        results: list[dict[str, Any]] = []

        try:
            max_workers = self._execution_controller.get_worker_config()["cpu"][
                "max_workers"
            ]
        except Exception:
            max_workers = 1

        if max_workers > 1:
            # Execução paralela: submeter todas as tarefas de uma vez

            # Get CPU configuration for worker processes
            cpu_config = self._execution_controller.create_worker_config()

            # Add memory configuration to be passed to workers
            worker_config = cpu_config.copy()
            if self._batch_config.resources and self._batch_config.resources.memory:
                worker_config["max_memory_gb"] = (
                    self._batch_config.resources.memory.max_memory_gb
                )

            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Preparar todas as tarefas
                futures = []
                for r in range(1, repetitions + 1):
                    if self._execution_controller.check_status() != BaseStatus.RUNNING:
                        break

                    # Construir unit_id padronizado: type:task:dataset:config:name:r
                    config_id = (
                        getattr(self._combination_store, "_preset_id", None)
                        or "default"
                    )
                    config_id = str(config_id).replace(":", "_")
                    experiment_unit_id = (
                        f"experiment:{task.id}:{dataset_id}:{config_id}:{alg.name}:{r}"
                    )

                    ex = self._combination_store.get_executions(
                        unit_id=experiment_unit_id
                    )
                    if ex and ex[0]["status"] in (
                        BaseStatus.COMPLETED.value,
                        "completed",
                        BaseStatus.FAILED.value,
                        "failed",
                    ):
                        continue

                    self._combination_store.submit_execution(
                        unit_id=experiment_unit_id, sequencia=r
                    )

                    # Calcular seed para esta repetição
                    seed = None
                    if global_seed is not None:
                        seed = global_seed
                        if self._batch_config.system.seed_increment:
                            seed += r

                    future = executor.submit(
                        _worker_exec,
                        experiment_unit_id,
                        seed,
                        strings,
                        alphabet,
                        self._batch_config.system.distance_method,
                        self._batch_config.system.enable_distance_cache,
                        alg.params,
                        self._combination_store.work_id,  # Usar work_id em vez de db_path
                        self._execution_controller.internal_jobs,
                        alg.name,  # Adicionar nome do algoritmo
                        worker_config.get(
                            "cpu"
                        ),  # Pass full worker configuration including memory
                    )
                    futures.append((r, future))

                # Coletar resultados
                for r, future in futures:
                    result, unit_id = future.result()

                    results.append({unit_id: result})

        else:
            # Execução sequencial

            # Get CPU configuration for sequential execution
            cpu_config = self._execution_controller.create_worker_config()

            # Add memory configuration to be passed to workers
            worker_config = cpu_config.copy()
            if self._batch_config.resources and self._batch_config.resources.memory:
                worker_config["max_memory_gb"] = (
                    self._batch_config.resources.memory.max_memory_gb
                )

            for r in range(1, repetitions + 1):
                status = self._execution_controller.check_status()
                if status != BaseStatus.RUNNING:
                    return status

                # Construir unit_id padronizado: type:task:dataset:config:name:r
                config_id = (
                    getattr(self._combination_store, "_preset_id", None) or "default"
                )
                config_id = str(config_id).replace(":", "_")
                experiment_unit_id = (
                    f"experiment:{task.id}:{dataset_id}:{config_id}:{alg.name}:{r}"
                )
                ex = self._combination_store.get_executions(unit_id=experiment_unit_id)
                if ex and ex[0]["status"] in (
                    BaseStatus.COMPLETED.value,
                    "completed",
                    BaseStatus.FAILED.value,
                    "failed",
                ):
                    continue

                self._combination_store.submit_execution(
                    unit_id=experiment_unit_id, sequencia=r
                )

                # Calcular seed para esta repetição
                seed = None
                if global_seed is not None:
                    seed = global_seed
                    if self._batch_config.system.seed_increment:
                        seed += r

                result, _ = _worker_exec(
                    experiment_unit_id,
                    seed,
                    strings,
                    alphabet,
                    self._batch_config.system.distance_method,
                    self._batch_config.system.enable_distance_cache,
                    alg.params,
                    self._combination_store.work_id,  # Usar work_id em vez de db_path
                    self._execution_controller.internal_jobs,
                    alg.name,  # Adicionar nome do algoritmo
                    worker_config.get(
                        "cpu"
                    ),  # Pass full worker configuration including memory
                )
                results.append({experiment_unit_id: result})

        status = self._execution_controller.check_status()
        if status != BaseStatus.RUNNING:
            return status

        # Avaliar resultados (considera enum BaseStatus)
        flat_results = [list(item.values())[0] for item in results]
        failed_repetitions = [
            r
            for r in flat_results
            if r.get("status")
            not in (BaseStatus.COMPLETED, BaseStatus.COMPLETED.value, "completed")
        ]
        final_result = BaseStatus.ERROR if failed_repetitions else BaseStatus.COMPLETED
        return final_result
