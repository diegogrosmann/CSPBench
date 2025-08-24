"""OptimizationExecutor: simplified executor for optimization tasks."""

from __future__ import annotations

from typing import Any, Callable

from src.application.ports.repositories import AbstractExecutionEngine
from src.infrastructure.persistence.work_state import (
    WorkStatePersistence,
    WorkScopedPersistence,
    UnitScopedPersistence,
)
from src.domain.config import OptimizationTask, ResourcesConfig, SystemConfig, AlgParams
from src.domain.dataset import Dataset
from src.infrastructure.logging_config import get_logger
from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.OptimizationExecutor")


class OptimizationExecutor(AbstractExecutionEngine):
    def __init__(self, work_store: WorkScopedPersistence | None = None):
        self.work_store = work_store

    def run(
        self,
        task: OptimizationTask,
        dataset_obj: Dataset,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        work_id: str | None = None,
        system_config: SystemConfig | None = None,
        check_control: Callable[[], str] | None = None,
        work_store: WorkScopedPersistence | None = None,
    ) -> dict[str, Any]:
        """Executa uma tentativa simplificada de otimização."""
        work_store = work_store or self.work_store

        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        params = dict(getattr(task, "params", {}) or {})

        unit_id = f"{task.id}:{getattr(dataset_obj, 'id', None)}:{getattr(task, 'name', 'opt')}"

        unit_store = None
        if work_store:
            unit_store = work_store.create_unit_scoped(unit_id)
            unit_store.started({"sequencia": 0})

        combination_id = None
        if work_store:
            try:
                # Obter combination_id se disponível
                combination_id = work_store.get_combination_id(
                    task.id,
                    getattr(dataset_obj, "id", "unknown"),
                    getattr(task, "preset_id", "default"),
                    getattr(task, "algorithm", getattr(alg, "name", "algorithm")),
                    "optimization",
                )
                if unit_store:
                    unit_store.record_execution_start(
                        combination_id=combination_id or 1,
                        sequencia=0,
                        params=params,
                    )
            except Exception:
                pass

        try:
            result = run_algorithm(
                getattr(task, "algorithm", getattr(alg, "name", "algorithm")),
                strings,
                alphabet,
                params,
                store=work_store.store if work_store else None,
                work_id=work_store.work_id if work_store else work_id,
                unit_id=unit_id,
            )
            if unit_store:
                unit_store.finished(result)
                try:
                    unit_store.record_execution_end(
                        result.get("status", "ok"),
                        result,
                        result.get("objective"),
                    )
                except Exception:
                    pass
            return {"status": "completed", "result": result}
        except Exception as exc:
            if unit_store:
                unit_store.error(exc)
            return {"status": "error", "error": str(exc)}
