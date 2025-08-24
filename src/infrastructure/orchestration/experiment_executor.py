"""ExperimentExecutor: executa ExperimentTask.

Paraleliza repetições via ProcessPool se ``max_workers > 1`` usando ``ProcessPoolExecutor``.
Implementa controles de recursos, resume functionality e salvamento completo de dados.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Optional

from src.domain.distance import create_distance_calculator
from src.domain.config import AlgParams, CSPBenchConfig, ExperimentTask
from src.domain.dataset import Dataset
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.application.ports.repositories import AbstractExecutionEngine
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
    ExecutionScopedPersistence,
)
from src.infrastructure.monitoring.persistence_monitor import PersistenceMonitor
from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
    WorkScopedPersistence,
)

from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.ExperimentExecutor")


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

        # Caminho do banco para recriar persistência em subprocessos
        db_path = Path(self._combination_store.store.db_path)  # type: ignore[attr-defined]

        def _worker_exec(
            experiment_unit_id: str,
            seed: int | None,
            strings: list[str],
            alphabet: str,
            distance_method: str,
            use_cache: bool,
            params: dict[str, Any],
            db_path: str,
            internal_jobs: int,
        ) -> tuple[dict[str, Any], str]:
            """Função executada em subprocesso (precisa ser picklable)."""
            from src.infrastructure.persistence.work_state.core import (
                WorkStatePersistence,
            )  # local import
            from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
                ExecutionScopedPersistence,
            )
            from src.infrastructure.monitoring.persistence_monitor import (
                PersistenceMonitor,
            )
            from src.infrastructure.execution_control import ExecutionController
            from src.domain.distance import create_distance_calculator
            from src.domain.status import BaseStatus as _BaseStatus

            # Recriar store e wrappers
            store = WorkStatePersistence(Path(db_path))
            execution_store = ExecutionScopedPersistence(store, experiment_unit_id)
            monitor = PersistenceMonitor(execution_store)

            # Controller mínimo para fornecer internal_jobs (sem controles cruzados)
            dummy_controller = ExecutionController(resources=None, check_control=None)
            dummy_controller._internal_jobs = internal_jobs  # type: ignore[attr-defined]

            distance_calc = create_distance_calculator(
                distance_method=distance_method, strings=strings, use_cache=use_cache
            )
            try:
                execution_store.update_execution_status(BaseStatus.RUNNING)
                monitor.on_progress(0.0, "Iniciando instância do algoritmo")
                result = run_algorithm(
                    algorithm_name=alg.name,
                    strings=strings,
                    alphabet=alphabet,
                    distance_calculator=distance_calc,
                    execution_controller=dummy_controller,
                    monitor=monitor,
                    seed=seed,
                    params=params,
                )

                monitor.on_progress(1.0, "Instância do algoritmo finalizada")
                status_value = result.get("status")
                if hasattr(status_value, "value"):
                    status_value = status_value.value
                execution_store.update_execution_status(
                    status_value, 
                    result=result.get('algorithm_result', None), 
                    params=result.get("actual_params"),
                    objective=result.get('algorithm_result', None).get('max_distance', None)
                )
                return result, experiment_unit_id
            except Exception as exc:  # pragma: no cover - defensive
                execution_store.update_execution_status(
                    "failed", result={"error": str(exc)}
                )
                return {
                    "status": _BaseStatus.FAILED,
                    "error": str(exc),
                }, experiment_unit_id

        results: list[dict[str, Any]] = []

        try:
            max_workers = self._execution_controller.get_worker_config()["cpu"]["max_workers"]
        except Exception:
            max_workers = 1

        if max_workers > 1:
            # Execução paralela: submeter todas as tarefas de uma vez
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Preparar todas as tarefas
                futures = []
                for r in range(1, repetitions + 1):
                    if self._execution_controller.check_status() != BaseStatus.RUNNING:
                        break

                    experiment_unit_id = (
                        f"experiment:{task.id}:{dataset_id}:{alg.name}:{r}"
                    )
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
                        str(db_path),
                        self._execution_controller.internal_jobs,
                    )
                    futures.append((r, future))

                # Coletar resultados
                for r, future in futures:
                    result, unit_id = future.result()

                    results.append({unit_id: result})

        else:
            # Execução sequencial
            for r in range(1, repetitions + 1):
                if self._execution_controller.check_status() != BaseStatus.RUNNING:
                    break

                experiment_unit_id = f"experiment:{task.id}:{dataset_id}:{alg.name}:{r}"
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
                        str(db_path),
                        self._execution_controller.internal_jobs,
                )
                results.append({experiment_unit_id: result})

        # Avaliar resultados (considera enum BaseStatus)
        flat_results = [list(item.values())[0] for item in results]
        failed_repetitions = [
            r
            for r in flat_results
            if r.get("status") not in (BaseStatus.COMPLETED, BaseStatus.COMPLETED.value, "completed")
        ]
        final_result = BaseStatus.ERROR if failed_repetitions else BaseStatus.COMPLETED
        return final_result
