"""OptimizationExecutor: random sampler inicial.
Assume task.parameters descreve espaços {param: {type: 'int'|'float'|'choice', low, high, choices}} e config.trials.
"""

from __future__ import annotations

import random
import time
from typing import Any, Callable

from src.domain.config import AlgParams, OptimizationTask, ResourcesConfig
from src.infrastructure.monitoring.monitor_interface import Monitor, NoOpMonitor
from .algorithm_runner import run_algorithm


class OptimizationExecutor:
    def run(
        self,
        task: OptimizationTask,
        dataset_obj,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        monitor: Monitor | None,
        writer=None,
        global_seed: int | None = None,
        check_control: Callable[[], str] | None = None,
        store=None,
    ) -> dict[str, Any]:
        """Run random search optimization trials for a given algorithm/dataset pair."""
        monitor = monitor or NoOpMonitor()
        trials = int(task.config.get("trials", 10) if task.config else 10)
        rng = random.Random(global_seed)
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
        param_space = task.parameters or {}
        best = float("inf")
        best_params: dict[str, Any] = {}

        for t in range(trials):
            if check_control:
                status = check_control()
                if status == "canceled":
                    monitor.log("warning", "Optimization canceled", {"task": task.id})
                    break
                while status == "paused":
                    time.sleep(0.3)
                    status = check_control()
                    if status == "canceled":
                        monitor.log(
                            "warning", "Optimization canceled", {"task": task.id}
                        )
                        break
                if status == "canceled":
                    break
            unit_id = f"{task.id}:{dataset_id}:{alg.name}:trial{t}"
            monitor.unit_started(unit_id, {"trial": t})
            sample_params = dict(alg.params)
            for pname, spec in param_space.items():
                if not spec:
                    continue
                ptype = spec.get("type")
                if ptype == "int":
                    sample_params[pname] = rng.randint(
                        spec.get("low", 0), spec.get("high", 10)
                    )
                elif ptype == "float":
                    sample_params[pname] = rng.uniform(
                        spec.get("low", 0.0), spec.get("high", 1.0)
                    )
                elif ptype == "choice":
                    sample_params[pname] = rng.choice(spec.get("choices", [0]))
            if store:
                try:
                    store.record_execution_start(
                        unit_id=unit_id,
                        task_id=task.id,
                        dataset_id=dataset_id,
                        algorithm=alg.name,
                        mode="optimization",
                        trial=t,
                        params=sample_params,
                    )
                except Exception:  # noqa: BLE001
                    pass
            res = run_algorithm(alg.name, strings, alphabet, sample_params)
            if res.get("status") == "ok" and res.get("objective", float("inf")) < best:
                best = res["objective"]
                best_params = sample_params
            monitor.unit_finished(unit_id, res)
            if writer:
                writer.append(
                    {
                        "task_id": task.id,
                        "dataset_id": dataset_id,
                        "algorithm": alg.name,
                        "unit_id": unit_id,
                        "trial": t,
                        "mode": "optimization",
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

        summary = {"trials": trials, "best_objective": best, "best_params": best_params}
        monitor.log(
            "info", "OptimizationExecutor finished", {**summary, "task": task.id}
        )
        return summary
