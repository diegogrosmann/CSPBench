"""SensitivityExecutor: simple sampling + variance ratio placeholder.
Assume task.parameters define base params; we vary each numeric param +-delta frac.
"""

from __future__ import annotations

import random
import time
from statistics import variance
from typing import Any, Callable

from src.domain.config import AlgParams, SensitivityTask, ResourcesConfig
from src.infrastructure.monitoring.monitor_interface import Monitor, NoOpMonitor
from .algorithm_runner import run_algorithm


class SensitivityExecutor:
    def run(
        self,
        task: SensitivityTask,
        dataset_obj,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        monitor: Monitor | None,
        writer=None,
        global_seed: int | None = None,
        check_control: Callable[[], str] | None = None,
        store=None,
    ) -> dict[str, Any]:
        """Run simple parameter sampling to estimate rough sensitivity via variance."""
        monitor = monitor or NoOpMonitor()
        base_params = dict(alg.params)
        if isinstance(dataset_obj, list):
            strings = dataset_obj
            alphabet = (
                "".join(sorted({c for s in strings for c in s})) if strings else ""
            )
            dataset_id = (
                f"ds_{__import__('hashlib').md5(''.join(strings).encode()).hexdigest()[:8]}"
                if strings
                else "ds_unknown"
            )
        else:
            strings = getattr(dataset_obj, "sequences", [])
            alphabet = (
                dataset_obj.metadata.get("alphabet")
                if hasattr(dataset_obj, "metadata")
                else ""
            )
            if not alphabet and strings:
                alphabet = "".join(sorted({c for s in strings for c in s}))
            dataset_id = getattr(dataset_obj, "id", "dataset")
        rng = random.Random(global_seed)
        param_space = task.parameters or {}
        samples = int(task.config.get("samples", 16) if task.config else 16)
        results: list[dict[str, Any]] = []

        for i in range(samples):
            if check_control:
                status = check_control()
                if status == "canceled":
                    monitor.log("warning", "Sensitivity canceled", {"task": task.id})
                    break
                while status == "paused":
                    time.sleep(0.3)
                    status = check_control()
                    if status == "canceled":
                        monitor.log(
                            "warning", "Sensitivity canceled", {"task": task.id}
                        )
                        break
                if status == "canceled":
                    break
            unit_id = f"{task.id}:{dataset_id}:{alg.name}:sample{i}"
            monitor.unit_started(unit_id, {"sample": i})
            sample_params = dict(base_params)
            for pname, spec in param_space.items():
                if not spec:
                    continue
                if "low" in spec and "high" in spec:
                    sample_params[pname] = rng.uniform(spec["low"], spec["high"])
            if store:
                try:
                    store.record_execution_start(
                        unit_id=unit_id,
                        task_id=task.id,
                        dataset_id=dataset_id,
                        algorithm=alg.name,
                        mode="sensitivity",
                        sample=i,
                        params=sample_params,
                    )
                except Exception:  # noqa: BLE001
                    pass
            res = run_algorithm(alg.name, strings, alphabet, sample_params)
            monitor.unit_finished(unit_id, res)
            results.append(
                {
                    "params": sample_params,
                    "objective": res.get("objective", float("inf")),
                }
            )
            if writer:
                writer.append(
                    {
                        "task_id": task.id,
                        "dataset_id": dataset_id,
                        "algorithm": alg.name,
                        "unit_id": unit_id,
                        "sample": i,
                        "mode": "sensitivity",
                        **{k: v for k, v in res.items() if k != "metadata"},
                    }
                )
            if store:
                try:
                    store.record_execution_end(
                        unit_id, res.get("status", "ok"), res, res.get("objective")
                    )
                except Exception:  # noqa: BLE001
                    pass

        param_scores: dict[str, float] = {}
        values = [r["objective"] for r in results]
        var_all = variance(values) if len(values) > 1 else 0.0
        for pname in param_space.keys():
            # Placeholder: each param gets same share currently if variance >0
            param_scores[pname] = var_all

        monitor.log(
            "info",
            "SensitivityExecutor finished",
            {"task": task.id, "samples": samples},
        )
        return {"samples": samples, "param_scores": param_scores}
