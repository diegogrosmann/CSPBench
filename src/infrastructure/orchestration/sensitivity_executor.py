"""SensitivityExecutor: simple sampling + variance ratio placeholder.
Assume task.parameters define base params; we vary each numeric param +-delta frac.
"""

from __future__ import annotations

import random
import time
from statistics import variance
from typing import Any, Callable

from src.domain.config import AlgParams, SensitivityTask, ResourcesConfig, SystemConfig
from src.domain.dataset import Dataset
from src.application.ports.repositories import AbstractExecutionEngine
from src.infrastructure.persistence.work_state import (
    WorkScopedPersistence,
)
from .algorithm_runner import run_algorithm


class SensitivityExecutor(AbstractExecutionEngine):
    def run(
        self,
        task: SensitivityTask,
        dataset_obj: Dataset,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        work_id: str | None = None,
        system_config: SystemConfig | None = None,
        check_control: Callable[[], str] | None = None,
        work_store: WorkScopedPersistence | None = None,
    ) -> dict[str, Any]:
        """Run simple parameter sampling to estimate rough sensitivity via variance."""
        base_params = dict(alg.params)

        # Extract data from Dataset object
        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        dataset_id = getattr(
            dataset_obj,
            "id",
            f"ds_{__import__('hashlib').md5(''.join(strings).encode()).hexdigest()[:8]}",
        )

        # Extract global_seed from system_config
        global_seed = None
        if system_config and hasattr(system_config, "global_seed"):
            global_seed = system_config.global_seed

        rng = random.Random(global_seed)
        param_space = task.parameters or {}
        samples = int(task.config.get("samples", 16) if task.config else 16)
        results: list[dict[str, Any]] = []

        sensitivity_unit_id = f"sensitivity:{task.id}:{dataset_id}:{alg.name}"
        sensitivity_store = None
        if work_store:
            sensitivity_store = work_store.create_unit_scoped(sensitivity_unit_id)
            sensitivity_store.started(
                {
                    "task_id": task.id,
                    "dataset_id": dataset_id,
                    "algorithm_name": alg.name,
                    "samples": samples,
                    "mode": "sensitivity",
                }
            )

        for i in range(samples):
            if check_control:
                status = check_control()
                if status == "canceled":
                    if sensitivity_store:
                        sensitivity_store.warning(
                            "Sensitivity analysis canceled",
                            {"task": task.id},
                        )
                    break

                while status == "paused":
                    time.sleep(0.3)
                    status = check_control()
                    if status == "canceled":
                        if sensitivity_store:
                            sensitivity_store.warning(
                                "Sensitivity analysis canceled during pause",
                                {"task": task.id},
                            )
                        break
                if status == "canceled":
                    break

            unit_id = f"{task.id}:{dataset_id}:{alg.name}:sample{i}"
            unit_store = None
            if work_store:
                unit_store = work_store.create_unit_scoped(unit_id)
                unit_store.started({"sequencia": i})

            sample_params = dict(base_params)
            for pname, spec in param_space.items():
                if not spec:
                    continue
                if "low" in spec and "high" in spec:
                    sample_params[pname] = rng.uniform(spec["low"], spec["high"])

            combination_id = None
            if work_store:
                try:
                    # Obter combination_id se disponível
                    combination_id = work_store.get_combination_id(
                        task.id,
                        dataset_id,
                        getattr(alg, "preset_id", "default"),
                        alg.name,
                        "sensitivity",
                    )
                    if unit_store:
                        unit_store.record_execution_start(
                            combination_id=combination_id or 1,
                            sequencia=i,
                            params=sample_params,
                        )
                except Exception:  # noqa: BLE001
                    pass

            res = run_algorithm(
                alg.name,
                strings,
                alphabet,
                sample_params,
                store=work_store.store if work_store else None,
                work_id=work_store.work_id if work_store else work_id,
                unit_id=unit_id,
            )

            if unit_store:
                unit_store.finished(res)

            results.append(
                {
                    "params": sample_params,
                    "objective": res.get("objective", float("inf")),
                }
            )

            if unit_store:
                try:
                    unit_store.record_execution_end(
                        res.get("status", "ok"), res, res.get("objective")
                    )
                except Exception:  # noqa: BLE001
                    pass

        param_scores: dict[str, float] = {}
        values = [r["objective"] for r in results]
        var_all = variance(values) if len(values) > 1 else 0.0
        for pname in param_space.keys():
            # Placeholder: each param gets same share currently if variance >0
            param_scores[pname] = var_all

        if sensitivity_store:
            sensitivity_store.finished(
                {
                    "samples": samples,
                    "param_scores": param_scores,
                    "status": "completed" if samples > 0 else "failed",
                }
            )

        return {
            "status": "completed" if samples > 0 else "failed",
            "samples": samples,
            "param_scores": param_scores,
        }
