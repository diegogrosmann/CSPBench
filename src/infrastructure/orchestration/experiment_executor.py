"""ExperimentExecutor: executa ExperimentTask.

Paraleliza repetições via ProcessPool se ``max_workers > 1`` usando ``ProcessPoolExecutor``.
Implementa controles de recursos, resume functionality e salvamento completo de dados.
"""

from __future__ import annotations

import time
from concurrent.futures import ProcessPoolExecutor
import hashlib
from typing import Any, Callable

from src.domain.config import AlgParams, ExperimentTask, ResourcesConfig, SystemConfig
from src.domain.dataset import Dataset
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.application.ports.repositories import AbstractExecutionEngine, AbstractStore
from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.ExperimentExecutor")


def _worker_exec(algorithm_name: str, strings: list[str], alphabet: str, params: dict[str, Any], unit_id: str, repetition: int, work_id: str | None = None):
    """Execução simples em worker (picklable)."""
    try:
        result = run_algorithm(algorithm_name, strings, alphabet, params)
        return result
    except Exception as exc:  # pragma: no cover - defensive
        return {"status": "error", "error_message": str(exc), "objective": float("inf")}


class ExperimentExecutor(AbstractExecutionEngine):
    def __init__(self, store: AbstractStore | None = None):
        self.store = store
        self.execution_controller: ExecutionController | None = None
        logger.info("ExperimentExecutor inicializado")

    def run(self, task: ExperimentTask, dataset_obj: Dataset, alg: AlgParams, resources: ResourcesConfig | None, work_id: str | None = None, system_config: SystemConfig | None = None, check_control: Callable[[], str] | None = None, store: AbstractStore | None = None) -> dict[str, Any]:
        """Executa o experimento (simplificado)."""
        store = store or self.store

        repetitions = max(1, int(getattr(task, "repetitions", 1)))
        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        dataset_id = getattr(dataset_obj, "id", None)

        results: list[dict[str, Any]] = []

        # Decide paralelismo simples
        max_workers = 1
        if resources and getattr(resources, "cpu", None):
            max_workers = getattr(resources.cpu, "max_workers", 1) or 1

        if store and work_id:
            store.log(work_id, "info", f"Starting experiment: {alg.name}", {
                "task": task.id, "dataset": dataset_id, "repetitions": repetitions
            })

        if max_workers <= 1 or repetitions == 1:
            for r in range(repetitions):
                if check_control and check_control() == "canceled":
                    break
                unit_id = f"{task.id}:{dataset_id}:{alg.name}:rep{r}"
                try:
                    if store and work_id:
                        store.unit_started(work_id, unit_id, {"sequencia": r})
                    if store:
                        try:
                            store.record_execution_start(
                                unit_id=unit_id,
                                task_id=task.id,
                                dataset_id=dataset_id,
                                algorithm=alg.name,
                                mode="experiment",
                                sequencia=r,
                                params=dict(alg.params)
                            )
                        except Exception:
                            pass
                    res = run_algorithm(alg.name, strings, alphabet, dict(alg.params))
                    if store and work_id:
                        store.unit_finished(work_id, unit_id, res)
                    results.append(res)
                    if store:
                        try:
                            store.record_execution_end(unit_id, res.get("status", "ok"), res, res.get("objective"))
                        except Exception:
                            pass
                except Exception as exc:
                    if store and work_id:
                        store.error(work_id, unit_id, exc)
                    results.append({"status": "error", "error_message": str(exc), "objective": float("inf")})
        else:
            with ProcessPoolExecutor(max_workers=max_workers) as pool:
                futures = []
                for r in range(repetitions):
                    unit_id = f"{task.id}:{dataset_id}:{alg.name}:rep{r}"
                    params = dict(alg.params)
                    futures.append(pool.submit(_worker_exec, alg.name, strings, alphabet, params, unit_id, r, work_id))
                for fut in futures:
                    try:
                        res = fut.result()
                        results.append(res)
                    except Exception as exc:
                        results.append({"status": "error", "error_message": str(exc), "objective": float("inf")})

        success = [r for r in results if r.get("status") == "ok"]
        best = min((r.get("objective", float("inf")) for r in success), default=float("inf"))

        overall = "completed" if len(success) == len(results) and results else ("partial" if len(success) else "failed")

        if store and work_id:
            store.log(work_id, "info", "ExperimentExecutor finished", {"task": task.id, "total": len(results)})

        return {"status": overall, "total_units": len(results), "success": len(success), "best_objective": best}