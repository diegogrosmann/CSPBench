"""ExperimentExecutor: executa ExperimentTask.

Paraleliza repetições via ProcessPool se ``internal_jobs > 1`` usando ``ProcessPoolExecutor``.
"""

from __future__ import annotations

import time
from concurrent.futures import ProcessPoolExecutor
import hashlib
from typing import Any, Callable

from src.domain.config import AlgParams, ExperimentTask, ResourcesConfig
from src.infrastructure.monitoring.monitor_interface import Monitor, NoOpMonitor
from src.infrastructure.resource_control import ResourceController
from .algorithm_runner import run_algorithm

_WORKER_CONTROLLER_CFG = None


def _worker_init(controller_cfg):  # executa em cada worker
    global _WORKER_CONTROLLER_CFG
    _WORKER_CONTROLLER_CFG = controller_cfg
    try:
        # Recria controller local e aplica limites CPU/memória
        rc = ResourceController(controller_cfg)
        rc.apply_cpu_limits()
        rc.apply_memory_limits()
    except Exception:  # noqa: BLE001
        pass


# Função nível módulo (picklable) para ProcessPool


def _worker_exec(
    algorithm_name: str, strings: list[str], alphabet: str, params: dict[str, Any]
):
    from .algorithm_runner import run_algorithm as _run

    return _run(algorithm_name, strings, alphabet, params)


class ExperimentExecutor:
    def run(
        self,
        task: ExperimentTask,
        dataset_obj,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        monitor: Monitor | None,
        writer=None,
        global_seed: int | None = None,
        check_control: Callable[[], str] | None = None,
        store=None,
    ) -> dict[str, Any]:
        monitor = monitor or NoOpMonitor()
        repetitions = max(1, task.repetitions)
        # Accept list[str] or dataset-like object
        if isinstance(dataset_obj, list):
            strings = dataset_obj
            alphabet = ""
            dataset_id = (
                f"ds_{hashlib.md5(''.join(strings).encode()).hexdigest()[:8]}"
                if strings
                else "ds_unknown"
            )
        else:
            strings = getattr(dataset_obj, "sequences", []) or (
                dataset_obj.metadata.get("sequences")
                if hasattr(dataset_obj, "metadata")
                else []
            )
            if not strings and hasattr(dataset_obj, "sequences"):
                strings = dataset_obj.sequences
            alphabet = getattr(dataset_obj, "alphabet", None) or (
                dataset_obj.metadata.get("alphabet")
                if hasattr(dataset_obj, "metadata")
                else ""
            )
            if hasattr(dataset_obj, "metadata") and "alphabet" in dataset_obj.metadata:
                alphabet = dataset_obj.metadata["alphabet"]
            dataset_id = getattr(dataset_obj, "id", "dataset")
        # Fallback se dataset_obj for object simplificado
        if not alphabet and strings:
            alphabet = "".join(sorted({c for s in strings for c in s}))

        per_item_timeout = None
        total_timeout = None
        max_workers = 1
        if resources is not None:
            per_item_timeout = resources.timeouts.timeout_per_item
            total_timeout = resources.timeouts.timeout_total_batch
            cpu_cfg = resources.cpu
            max_workers = max(
                1,
                min(cpu_cfg.internal_jobs, cpu_cfg.max_cores or cpu_cfg.internal_jobs),
            )

        start_task = time.time()
        results: list[dict[str, Any]] = []

        controller_cfg = {
            "cpu": {"max_cores": max_workers, "affinity": None},
            "memory": {
                "max_memory_gb": (
                    resources.memory.max_memory_gb
                    if resources and resources.memory
                    else None
                )
            },
            "parallel": {},
            "timeouts": {
                "timeout_per_algorithm": per_item_timeout or 3600,
                "timeout_total_batch": total_timeout or 86400,
            },
        }
        ResourceController(controller_cfg)  # aplica limites master (best-effort)

        monitor.log(
            "info",
            f"ExperimentExecutor start alg={alg.name} reps={repetitions}",
            {"task": task.id},
        )

        def submit_seq(rep_index: int):
            params = dict(alg.params)
            if global_seed is not None:
                params.setdefault("seed", (global_seed + rep_index) & 0x7FFFFFFF)
            unit_id = f"{task.id}:{dataset_id}:{alg.name}:rep{rep_index}"
            if store:
                try:
                    store.record_execution_start(
                        unit_id=unit_id,
                        task_id=task.id,
                        dataset_id=dataset_id,
                        algorithm=alg.name,
                        mode="experiment",
                        repetition=rep_index,
                        params=params,
                    )
                except Exception:  # noqa: BLE001
                    pass
            monitor.unit_started(unit_id, {"repetition": rep_index})
            result = run_algorithm(alg.name, strings, alphabet, params)
            monitor.unit_finished(unit_id, result)
            if writer:
                writer.append(
                    {
                        "task_id": task.id,
                        "dataset_id": dataset_id,
                        "algorithm": alg.name,
                        "unit_id": unit_id,
                        "mode": "experiment",
                        **{k: v for k, v in result.items() if k not in {"metadata"}},
                    }
                )
            if store:
                try:
                    store.record_execution_end(
                        unit_id,
                        result.get("status", "ok"),
                        result,
                        result.get("objective"),
                    )
                except Exception:  # noqa: BLE001
                    pass
            results.append(result)

        def _wait_ok() -> bool:
            if not check_control:
                return True
            status = check_control()
            if status == "canceled":
                return False
            while status == "paused":
                time.sleep(0.3)
                status = check_control()
                if status == "canceled":
                    return False
            return True

        if max_workers == 1 or repetitions == 1:
            for r in range(repetitions):
                if total_timeout and (time.time() - start_task) > total_timeout:
                    monitor.log("warning", "Task timeout reached", {"task": task.id})
                    break
                if not _wait_ok():
                    monitor.log(
                        "warning", "Task canceled (sequential)", {"task": task.id}
                    )
                    break
                submit_seq(r)
        else:
            with ProcessPoolExecutor(
                max_workers=max_workers,
                initializer=_worker_init,
                initargs=(controller_cfg,),
            ) as pool:
                futures = []
                for r in range(repetitions):
                    if total_timeout and (time.time() - start_task) > total_timeout:
                        monitor.log(
                            "warning", "Task timeout reached", {"task": task.id}
                        )
                        break
                    if not _wait_ok():
                        monitor.log(
                            "warning",
                            "Task canceled before dispatch (parallel)",
                            {"task": task.id},
                        )
                        break
                    params = dict(alg.params)
                    if global_seed is not None:
                        params.setdefault("seed", (global_seed + r) & 0x7FFFFFFF)
                    unit_id = f"{task.id}:{dataset_id}:{alg.name}:rep{r}"
                    if store:
                        try:
                            store.record_execution_start(
                                unit_id=unit_id,
                                task_id=task.id,
                                dataset_id=dataset_id,
                                algorithm=alg.name,
                                mode="experiment",
                                repetition=r,
                                params=params,
                            )
                        except Exception:  # noqa: BLE001
                            pass
                    monitor.unit_started(unit_id, {"repetition": r})
                    f = pool.submit(_worker_exec, alg.name, strings, alphabet, params)
                    futures.append((unit_id, f))
                for unit_id, fut in futures:
                    if not _wait_ok():
                        monitor.log(
                            "warning",
                            "Task canceled during result collection",
                            {"task": task.id},
                        )
                        fut.cancel()
                        continue
                    try:
                        res = fut.result(timeout=per_item_timeout)
                    except Exception as e:  # noqa: BLE001
                        # Tentativa de cancelamento; se já iniciou execução em outro processo não garante kill.
                        fut.cancel()
                        # Hard kill best-effort: ProcessPoolExecutor não expõe PID diretamente; opção seria migrar para multiprocessing manual.
                        res = {
                            "status": "error",
                            "error_type": e.__class__.__name__,
                            "error_message": f"timeout/cancel: {e}",
                            "objective": float("inf"),
                        }
                        monitor.error(unit_id, e)
                    monitor.unit_finished(unit_id, res)
                    if writer:
                        writer.append(
                            {
                                "task_id": task.id,
                                "dataset_id": dataset_id,
                                "algorithm": alg.name,
                                "unit_id": unit_id,
                                "mode": "experiment",
                                **{
                                    k: v
                                    for k, v in res.items()
                                    if k not in {"metadata"}
                                },
                            }
                        )
                    if store:
                        try:
                            store.record_execution_end(
                                unit_id,
                                res.get("status", "ok"),
                                res,
                                res.get("objective"),
                            )
                        except Exception:  # noqa: BLE001
                            pass
                    results.append(res)
        # Agregar estatísticas simples
        success = [r for r in results if r["status"] == "ok"]
        best = min((r["objective"] for r in success), default=float("inf"))
        agg = {
            "total_units": len(results),
            "success": len(success),
            "failures": len(results) - len(success),
            "best_objective": best,
        }
        monitor.log("info", "ExperimentExecutor finished", {**agg, "task": task.id})
        return agg
