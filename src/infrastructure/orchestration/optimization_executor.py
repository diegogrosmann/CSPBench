"""OptimizationExecutor: simplified executor for optimization tasks."""

from __future__ import annotations

from typing import Any, Callable

from src.application.ports.repositories import AbstractExecutionEngine, AbstractStore
from src.domain.config import OptimizationTask, ResourcesConfig, SystemConfig, AlgParams
from src.domain.dataset import Dataset
from src.infrastructure.logging_config import get_logger
from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.OptimizationExecutor")


class OptimizationExecutor(AbstractExecutionEngine):
    def __init__(self, store: AbstractStore | None = None):
        self.store = store

    def run(self, task: OptimizationTask, dataset_obj: Dataset, alg: AlgParams, resources: ResourcesConfig | None, work_id: str | None = None, system_config: SystemConfig | None = None, check_control: Callable[[], str] | None = None, store: AbstractStore | None = None) -> dict[str, Any]:
        """Executa uma tentativa simplificada de otimização."""
        store = store or self.store

        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        params = dict(getattr(task, "params", {}) or {})

        unit_id = f"{task.id}:{getattr(dataset_obj, 'id', None)}:{getattr(task, 'name', 'opt')}"
        
        if store and work_id:
            store.unit_started(work_id, unit_id, {"sequencia": 0})

        if store:
            try:
                store.record_execution_start(
                    unit_id=unit_id,
                    task_id=task.id,
                    dataset_id=getattr(dataset_obj, 'id', None),
                    algorithm=getattr(task, "algorithm", getattr(alg, "name", "algorithm")),
                    mode="optimization",
                    sequencia=0,
                    params=params
                )
            except Exception:
                pass

        try:
            result = run_algorithm(getattr(task, "algorithm", getattr(alg, "name", "algorithm")), strings, alphabet, params)
            if store and work_id:
                store.unit_finished(work_id, unit_id, result)
            if store:
                try:
                    store.record_execution_end(unit_id, result.get("status", "ok"), result, result.get("objective"))
                except Exception:
                    pass
            return {"status": "completed", "result": result}
        except Exception as exc:
            if store and work_id:
                store.error(work_id, unit_id, exc)
            return {"status": "error", "error": str(exc)}