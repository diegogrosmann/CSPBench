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
from src.infrastructure.monitoring.monitor_interface import Monitor, NoOpMonitor
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.application.ports.repositories import AbstractExecutionEngine, AbstractStore
from .algorithm_runner import run_algorithm

# Create module logger
logger = get_logger("CSPBench.ExperimentExecutor")

_WORKER_CONTROLLER_CFG = None
_WORKER_MONITOR = None


def _worker_init(controller_cfg, monitor_instance=None):  # executa em cada worker
    global _WORKER_CONTROLLER_CFG, _WORKER_MONITOR
    _WORKER_CONTROLLER_CFG = controller_cfg
    _WORKER_MONITOR = monitor_instance or NoOpMonitor()
    
    # Create worker-specific logger
    worker_logger = get_logger("CSPBench.ExperimentExecutor.Worker")
    worker_logger.debug("Inicializando worker do ExperimentExecutor")
    
    try:
        # Recria controller local e aplica limites CPU/memória
        from src.infrastructure.execution_control import ExecutionController
        ec = ExecutionController.from_config(controller_cfg)
        ec.apply_cpu_limits()
        ec.apply_memory_limits()
        worker_logger.info("Controles de recursos aplicados no worker")
    except Exception as e:  # noqa: BLE001
        worker_logger.warning(f"Falha ao aplicar controles de recursos no worker: {e}")


# Função nível módulo (picklable) para ProcessPool


def _worker_exec(
    algorithm_name: str, 
    strings: list[str], 
    alphabet: str, 
    params: dict[str, Any],
    unit_id: str,
    repetition: int
):
    global _WORKER_MONITOR
    from .algorithm_runner import run_algorithm as _run

    # Create worker-specific logger
    worker_logger = get_logger("CSPBench.ExperimentExecutor.Worker")
    worker_logger.info(f"Worker iniciando execução: algoritmo={algorithm_name}, repetição={repetition}")

    monitor = _WORKER_MONITOR or NoOpMonitor()
    
    # Agora o unit_started é chamado quando realmente inicia no worker
    monitor.unit_started(unit_id, {"repetition": repetition})
    
    try:
        worker_logger.debug(f"Executando algoritmo {algorithm_name} com {len(strings)} strings")
        result = _run(algorithm_name, strings, alphabet, params)
        
        worker_logger.info(f"Algoritmo {algorithm_name} (rep {repetition}) executado: status={result.get('status')}")
        monitor.unit_finished(unit_id, result)
        return result
    except Exception as e:
        worker_logger.error(f"Erro na execução do worker: {e}", exc_info=True)
        error_result = {
            "status": "error",
            "error_type": e.__class__.__name__,
            "error_message": str(e),
            "objective": float("inf"),
        }
        monitor.unit_finished(unit_id, error_result)
        monitor.error(unit_id, e)
        return error_result


class ExperimentExecutor(AbstractExecutionEngine):
    def __init__(self, store: AbstractStore | None = None):
        """Initialize ExperimentExecutor with typed store."""
        self.store = store
        self.execution_controller: ExecutionController | None = None
        logger.info("ExperimentExecutor inicializado")

    def run(
        self,
        task: ExperimentTask,
        dataset_obj: Dataset,
        alg: AlgParams,
        resources: ResourcesConfig | None,
        monitor: Monitor | None = None,
        system_config: SystemConfig | None = None,
        check_control: Callable[[], str] | None = None,
        store: AbstractStore | None = None,
    ) -> dict[str, Any]:
        logger.info(f"Iniciando experimento: {task.name}")
        logger.info(f"Dataset: {dataset_obj.name}, Algoritmo: {alg.name}, Repetições: {task.repetitions}")
        
        # Use store from parameter or instance
        store = store or self.store
        monitor = monitor or NoOpMonitor()
        repetitions = max(1, task.repetitions)
        
        logger.debug(f"Configuração do experimento: repetições={repetitions}")
        
        # Extract global_seed from system_config
        global_seed = None
        if system_config and hasattr(system_config, 'global_seed'):
            global_seed = system_config.global_seed
            logger.debug(f"Seed global configurada: {global_seed}")
        
        # Extract data from Dataset object
        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        dataset_id = getattr(dataset_obj, "id", None)

        # Check for completed repetitions (resume functionality)
        start_repetition = 0
        completed_results = []
        if store and dataset_id:
            completed_reps = store.get_completed_repetitions(task.id, dataset_id, alg.name)
            start_repetition = len(completed_reps)
            completed_results = [rep["result"] for rep in completed_reps]
            if start_repetition > 0:
                monitor.log("info", f"Resuming from repetition {start_repetition}", {"task": task.id})

        per_item_timeout = None
        total_timeout = None
        max_workers = None
        if resources is not None:
            per_item_timeout = resources.timeouts.timeout_per_item
            total_timeout = resources.timeouts.timeout_total_batch
            cpu_cfg = resources.cpu
            # Let ResourceController handle None values and set defaults
            max_workers = cpu_cfg.max_workers or cpu_cfg.internal_jobs

        # Apply resource controls if ExecutionController is available
        if self.execution_controller:
            self.execution_controller.apply_cpu_limits()
            self.execution_controller.apply_memory_limits()
            # Get max_workers from ExecutionController
            max_workers = self.execution_controller.max_workers
        else:
            # Fallback to original behavior
            if resources:
                self._apply_resource_controls(resources, max_workers)
            # ResourceController will set defaults, but we need a fallback for executor parallelization
            if max_workers is None:
                import os
                max_workers = os.cpu_count() or 1
            else:
                max_workers = max(1, max_workers)

        start_task = time.time()
        results: list[dict[str, Any]] = completed_results.copy()

        # Use ExecutionController configuration if available
        if self.execution_controller:
            controller_cfg = self.execution_controller.get_worker_config()
        else:
            controller_cfg = {
                "cpu": {
                    "max_workers": resources.cpu.max_workers if resources else None,
                    "exclusive_cores": resources.cpu.exclusive_cores if resources else False,
                    "internal_jobs": resources.cpu.internal_jobs if resources else 1,
                    "affinity": None  # Let ExecutionController calculate this
                },
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

        monitor.log(
            "info",
            f"ExperimentExecutor start alg={alg.name} reps={repetitions}",
            {"task": task.id},
        )

        def submit_seq(rep_index: int):
            # Prepare parameters including internal_jobs for algorithm
            params = dict(alg.params)
            if global_seed is not None:
                params.setdefault("seed", (global_seed + rep_index) & 0x7FFFFFFF)
            
            # Pass internal_jobs to algorithm (not for executor parallelization)
            if self.execution_controller:
                internal_jobs = self.execution_controller.internal_jobs
                if internal_jobs:
                    params.setdefault("internal_jobs", internal_jobs)
            elif resources and resources.cpu.internal_jobs:
                params.setdefault("internal_jobs", resources.cpu.internal_jobs)
            
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
            if store:
                try:
                    store.record_execution_end(
                        unit_id,
                        result.get("status", "ok"),
                        result,  # Now includes extended data (history, metadata, actual_params)
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
            for r in range(start_repetition, repetitions):
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
                initargs=(controller_cfg, monitor),  # Passa monitor para workers
            ) as pool:
                futures = []
                for r in range(start_repetition, repetitions):
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
                    
                    # Pass internal_jobs to algorithm (not for executor parallelization)
                    if resources and resources.cpu.internal_jobs:
                        params.setdefault("internal_jobs", resources.cpu.internal_jobs)
                    
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
                    # Submit job to pool (will be queued if no workers available)
                    f = pool.submit(_worker_exec, alg.name, strings, alphabet, params, unit_id, r)
                    futures.append((unit_id, f))
                
                # Collect results (may complete in different order than submission)
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
                    # unit_finished já foi chamado dentro do worker
                    if store:
                        try:
                            store.record_execution_end(
                                unit_id,
                                res.get("status", "ok"),
                                res,  # Now includes extended data (history, metadata, actual_params)
                                res.get("objective"),
                            )
                        except Exception:  # noqa: BLE001
                            pass
                    results.append(res)
        # Agregar estatísticas simples
        success = [r for r in results if r["status"] == "ok"]
        best = min((r["objective"] for r in success), default=float("inf"))
        
        # Determine overall status based on results
        if len(success) == len(results):
            overall_status = "completed"
        elif len(success) > 0:
            overall_status = "partial"
        else:
            overall_status = "failed"
        
        agg = {
            "status": overall_status,
            "total_units": len(results),
            "success": len(success),
            "failures": len(results) - len(success),
            "best_objective": best,
        }
        monitor.log("info", "ExperimentExecutor finished", {**agg, "task": task.id})
        return agg
    
    def _apply_resource_controls(self, resources: ResourcesConfig, max_workers: int | None) -> None:
        """Apply resource controls by creating temporary ExecutionController."""
        try:
            if resources:
                # Create temporary ExecutionController for resource control
                temp_controller = ExecutionController(resources=resources)
                temp_controller.apply_cpu_limits()
                temp_controller.apply_memory_limits()
        except Exception:  # noqa: BLE001
            pass  # Best effort resource control
