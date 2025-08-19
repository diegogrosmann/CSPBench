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
from src.infrastructure.monitoring.monitor_interface import Monitor, NoOpMonitor
from src.application.ports.repositories import AbstractExecutionEngine
from .algorithm_runner import run_algorithm


class SensitivityExecutor(AbstractExecutionEngine):
    def run(
        self,
        task: SensitivityTask,
        dataset_obj: Dataset,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        monitor: Monitor | None = None,
        system_config: SystemConfig | None = None,
        check_control: Callable[[], str] | None = None,
        store=None,
    ) -> dict[str, Any]:
        """Run simple parameter sampling to estimate rough sensitivity via variance."""
        monitor = monitor or NoOpMonitor()
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
        return {
            "status": "completed" if samples > 0 else "failed",
            "samples": samples,
            "param_scores": param_scores,
        }
