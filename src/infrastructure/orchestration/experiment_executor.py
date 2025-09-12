"""ExperimentExecutor: executa ExperimentTask.

Paraleliza repetições via ProcessPool se ``max_workers > 1`` usando ``ProcessPoolExecutor``.
Implementa controles de recursos, resume functionality e salvamento completo de dados.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from src.application.ports.repositories import AbstractExecutionEngine
from src.domain.config import AlgParams, CSPBenchConfig, ExperimentTask
from src.domain.dataset import Dataset
from src.domain.distance import create_distance_calculator
from src.domain.status import BaseStatus
from src.infrastructure.execution_control import ExecutionController
from src.infrastructure.logging_config import get_logger
from src.infrastructure.monitoring.persistence_monitor import PersistenceMonitor
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import (
    CombinationScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
    ExecutionScopedPersistence,
)
from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
    WorkScopedPersistence,
)

from .algorithm_runner import run_algorithm

logger = get_logger("CSPBench.ExperimentExecutor")


import random
import time

def _worker_exec(
    experiment_unit_id: str,
    seed: int | None,
    strings: list[str],
    alphabet: str,
    distance_method: str,
    use_cache: bool,
    params: dict[str, Any],
    work_id: str,  # Alterado de db_path para work_id
    internal_jobs: int,
    algorithm_name: str,
    cpu_config: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], str]:
    """Função executada em subprocesso (precisa ser picklable)."""
    import psutil

    from src.domain.status import BaseStatus as _BaseStatus
    from src.infrastructure.execution_control import ExecutionController
    from src.infrastructure.logging_config import get_logger
    from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
        WorkScopedPersistence,
    )

    logger = get_logger("CSPBench.WorkerProcess")

    logger.debug(f"[WORKER][{experiment_unit_id}] Iniciando worker_exec.")

    # Apply CPU configuration to this worker process
    if cpu_config:
        try:
            current_process = psutil.Process()

            # Apply CPU affinity if specified
            if cpu_config.get("affinity"):
                current_process.cpu_affinity(cpu_config["affinity"])
                logger.info(f"[WORKER][{experiment_unit_id}] CPU affinity set to: {cpu_config['affinity']}")

            # Apply CPU priority if max_workers is limited
            if cpu_config.get("exclusive_cores", False):
                current_process.nice(5)  # Slightly lower priority than main process
                logger.info(f"[WORKER][{experiment_unit_id}] CPU priority adjusted for exclusive cores")

        except Exception as e:
            logger.warning(f"[WORKER][{experiment_unit_id}] Cannot apply CPU configuration: {e}")

    try:
        # Usar factory method para criar work_scoped diretamente
        work_scoped = WorkScopedPersistence.submit(work_id)
        execution_store = work_scoped.for_execution(experiment_unit_id)
    except Exception as e:
        logger.error(f"[WORKER][{experiment_unit_id}] Erro ao recriar store: {e}")
        raise

    # Controller with work_id for status checks
    dummy_controller = ExecutionController(work_id=work_id)
    dummy_controller._internal_jobs = internal_jobs  # type: ignore[attr-defined]

    # Apply CPU configuration manually if provided (since we can't pass ResourcesConfig to subprocess)
    if cpu_config:
        dummy_controller._exclusive_cores = cpu_config.get("exclusive_cores", False)
        if "max_workers" in cpu_config:
            dummy_controller._max_workers = cpu_config["max_workers"]
            dummy_controller._current_workers = cpu_config["max_workers"]

        # Try to apply CPU affinity in worker process (best effort)
        try:
            if cpu_config.get("exclusive_cores", False) and hasattr(
                dummy_controller, "apply_memory_limits"
            ):
                dummy_controller.apply_memory_limits()
        except Exception as cpu_exc:
            logger.warning(
                f"[WORKER][{experiment_unit_id}] CPU configuration parcialmente aplicada: {cpu_exc}"
            )

    # Create monitor with controller for cancellation checks
    monitor = PersistenceMonitor(execution_store, execution_controller=dummy_controller)

    distance_calc = create_distance_calculator(
        distance_method=distance_method, strings=strings, use_cache=use_cache
    )
    try:
        execution_store.update_execution_status(_BaseStatus.RUNNING)
        monitor.on_progress(0.0, f"Iniciando instância do algoritmo")

        result = run_algorithm(
            algorithm_name=algorithm_name,
            strings=strings,
            alphabet=alphabet,
            distance_calculator=distance_calc,
            execution_controller=dummy_controller,
            monitor=monitor,
            seed=seed,
            params=params,
        )

        monitor.on_progress(1.0, f"[{experiment_unit_id}] Instância do algoritmo finalizada")
        status_value = result.get("status")
        if hasattr(status_value, "value"):
            status_value = status_value.value
        objective = None
        alg_result = result.get("algorithm_result", None)
        if alg_result is not None:
            max_distance = alg_result.get("max_distance", None)
            if max_distance is not None:
                objective = max_distance

        execution_store.update_execution_status(
            status=status_value,
            result=alg_result,
            params=result.get("actual_params", None),
            objective=objective,
        )
        logger.debug(f"[WORKER][{experiment_unit_id}] worker_exec finalizado com status: {status_value}")
        return result, experiment_unit_id
    except Exception as exc:  # pragma: no cover - defensive
        execution_store.update_execution_status("failed", result={"error": str(exc)})
        logger.error(f"[WORKER][{experiment_unit_id}] worker_exec falhou: {exc}")
        return {
            "status": _BaseStatus.FAILED,
            "error": str(exc),
        }, experiment_unit_id


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
        """Executa o experimento seguindo o novo processo passo a passo."""

        repetitions = max(1, int(getattr(task, "repetitions", 1)))
        strings = dataset_obj.sequences
        alphabet = dataset_obj.alphabet
        dataset_id = getattr(dataset_obj, "id", None)

        global_seed = getattr(self._batch_config.system, "global_seed", None)
        if global_seed is None:
            global_seed = alg.params.get("seed", None)

        logger.info(f"Iniciando experimento para {alg.name} no dataset {dataset_id} com {repetitions} repetições")

        # PASSO 1: Identificar todas as execuções necessárias
        logger.info("PASSO 1: Identificando execuções necessárias")
        all_needed_executions = self._identify_needed_executions(
            task, dataset_id, alg, repetitions
        )
        
        # PASSO 2: Listar execuções já finalizadas 
        logger.info("PASSO 2: Verificando execuções já finalizadas")
        completed_executions = self._get_completed_executions(all_needed_executions)
        
        # PASSO 3: Determinar execuções pendentes
        logger.info("PASSO 3: Determinando execuções pendentes")
        pending_executions = self._get_pending_executions(
            all_needed_executions, completed_executions
        )
        
        if not pending_executions:
            logger.info("Nenhuma execução pendente encontrada - todas já foram finalizadas")
            return BaseStatus.COMPLETED
            
        logger.info(f"Encontradas {len(pending_executions)} execuções pendentes de {len(all_needed_executions)} totais")
        
        # PASSO 4: Persistir execuções pendentes em lote
        logger.info("PASSO 4: Persistindo execuções pendentes em lote")
        self._persist_pending_executions(pending_executions)
        
        # PASSO 5: Executar as execuções pendentes
        logger.info("PASSO 5: Executando execuções pendentes")
        results = self._execute_pending_executions(
            pending_executions, strings, alphabet, global_seed
        )

        status = self._execution_controller.check_status()
        if status != BaseStatus.RUNNING:
            return status

        # Avaliar resultados (considera enum BaseStatus)
        if not results:
            return BaseStatus.COMPLETED
            
        flat_results = [list(item.values())[0] for item in results]
        failed_repetitions = [
            r
            for r in flat_results
            if r.get("status")
            not in (BaseStatus.COMPLETED, BaseStatus.COMPLETED.value, "completed")
        ]
        final_result = BaseStatus.ERROR if failed_repetitions else BaseStatus.COMPLETED
        return final_result

    def _identify_needed_executions(
        self, task: ExperimentTask, dataset_id: str, alg: AlgParams, repetitions: int
    ) -> list[dict[str, Any]]:
        """PASSO 1: Identifica todas as execuções que precisam ser realizadas."""
        needed_executions = []
        
        config_id = (
            getattr(self._combination_store, "_preset_id", None) or "default"
        )
        config_id = str(config_id).replace(":", "_")
        
        for r in range(1, repetitions + 1):
            experiment_unit_id = (
                f"experiment:{task.id}:{dataset_id}:{config_id}:{alg.name}:{r}"
            )
            
            execution_data = {
                "unit_id": experiment_unit_id,
                "sequencia": r,
                "task_id": task.id,
                "dataset_id": dataset_id,
                "config_id": config_id,
                "algorithm_name": alg.name,
                "algorithm_params": alg.params
            }
            needed_executions.append(execution_data)
            
        logger.debug(f"Identificadas {len(needed_executions)} execuções necessárias")
        return needed_executions

    def _get_completed_executions(self, all_executions: list[dict[str, Any]]) -> set[str]:
        """PASSO 2: Lista todas as execuções que já estão com status finalizado."""
        # Status finais que consideramos como "completos"
        final_statuses = [
            BaseStatus.COMPLETED.value,
            "completed",
            BaseStatus.FAILED.value,
            "failed", 
            BaseStatus.ERROR.value,
            "error",
        ]
        
        # Buscar todas as execuções da combinação atual de uma só vez
        all_existing_executions = self._combination_store.get_executions()
        
        # Filtrar apenas execuções com status finais
        completed_executions = [
            ex for ex in all_existing_executions 
            if ex.get("status") in final_statuses
        ]
        
        # Criar um set com os unit_ids das execuções já finalizadas
        completed_unit_ids = {ex["unit_id"] for ex in completed_executions}
        
        # Filtrar apenas as que estão na nossa lista de execuções necessárias
        needed_unit_ids = {ex["unit_id"] for ex in all_executions}
        relevant_completed = completed_unit_ids.intersection(needed_unit_ids)
        
        logger.debug(f"Encontradas {len(relevant_completed)} execuções já finalizadas de {len(completed_executions)} com status final, de {len(all_existing_executions)} totais no banco")
        return relevant_completed

    def _get_pending_executions(
        self, all_executions: list[dict[str, Any]], completed: set[str]
    ) -> list[dict[str, Any]]:
        """PASSO 3: Compara e retorna apenas as execuções que não foram realizadas."""
        pending = []
        
        for execution in all_executions:
            if execution["unit_id"] not in completed:
                pending.append(execution)
                
        logger.debug(f"Determinadas {len(pending)} execuções pendentes")
        return pending

    def _persist_pending_executions(self, pending_executions: list[dict[str, Any]]) -> None:
        """PASSO 4: Persiste a lista de execuções pendentes em lote."""
        if not pending_executions:
            return
            
        # Preparar dados para inserção em lote
        executions_to_insert = []
        for execution in pending_executions:
            execution_data = {
                "unit_id": execution["unit_id"],
                "sequencia": execution["sequencia"],
                # combination_id será adicionado automaticamente pelo submit_executions
            }
            executions_to_insert.append(execution_data)
            
        # Usar inserção em lote para melhor performance
        inserted_count = self._combination_store.submit_executions(executions_to_insert)
        logger.info(f"Persistidas {inserted_count} execuções pendentes em lote")

    def _execute_pending_executions(
        self, 
        pending_executions: list[dict[str, Any]], 
        strings: list[str], 
        alphabet: str,
        global_seed: int | None
    ) -> list[dict[str, Any]]:
        """PASSO 5: Executa a lista de execuções pendentes (sequencial ou paralela)."""
        if not pending_executions:
            return []
            
        results: list[dict[str, Any]] = []

        try:
            max_workers = self._execution_controller.get_worker_config()["cpu"][
                "max_workers"
            ]
        except Exception:
            max_workers = 1

        if max_workers > 1:
            # Execução paralela
            logger.info(f"Executando {len(pending_executions)} execuções em paralelo com {max_workers} workers")
            results = self._execute_parallel(pending_executions, strings, alphabet, global_seed)
        else:
            # Execução sequencial
            logger.info(f"Executando {len(pending_executions)} execuções sequencialmente")
            results = self._execute_sequential(pending_executions, strings, alphabet, global_seed)
            
        return results

    def _execute_parallel(
        self, 
        pending_executions: list[dict[str, Any]], 
        strings: list[str], 
        alphabet: str,
        global_seed: int | None
    ) -> list[dict[str, Any]]:
        """Executa as execuções pendentes em paralelo."""
        results = []
        
        # Get CPU configuration for worker processes
        cpu_config = self._execution_controller.create_worker_config()
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        try:
            max_workers = self._execution_controller.get_worker_config()["cpu"][
                "max_workers"
            ]
        except Exception:
            max_workers = 1

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            
            for execution in pending_executions:
                if self._execution_controller.check_status() != BaseStatus.RUNNING:
                    break

                # Calcular seed para esta repetição
                seed = None
                if global_seed is not None:
                    seed = global_seed
                    if self._batch_config.system.seed_increment:
                        seed += execution["sequencia"]

                future = executor.submit(
                    _worker_exec,
                    execution["unit_id"],
                    seed,
                    strings,
                    alphabet,
                    self._batch_config.system.distance_method,
                    self._batch_config.system.enable_distance_cache,
                    execution["algorithm_params"],
                    self._combination_store.work_id,
                    self._execution_controller.internal_jobs,
                    execution["algorithm_name"],
                    worker_config.get("cpu"),
                )
                futures.append((execution["sequencia"], future))

            # Coletar resultados
            for sequencia, future in futures:
                result, unit_id = future.result()
                results.append({unit_id: result})
                
        return results

    def _execute_sequential(
        self, 
        pending_executions: list[dict[str, Any]], 
        strings: list[str], 
        alphabet: str,
        global_seed: int | None
    ) -> list[dict[str, Any]]:
        """Executa as execuções pendentes sequencialmente."""
        results = []
        
        # Get CPU configuration for sequential execution
        cpu_config = self._execution_controller.create_worker_config()
        worker_config = cpu_config.copy()
        if self._batch_config.resources and self._batch_config.resources.memory:
            worker_config["max_memory_gb"] = (
                self._batch_config.resources.memory.max_memory_gb
            )

        for execution in pending_executions:
            status = self._execution_controller.check_status()
            if status != BaseStatus.RUNNING:
                break

            # Calcular seed para esta repetição
            seed = None
            if global_seed is not None:
                seed = global_seed
                if self._batch_config.system.seed_increment:
                    seed += execution["sequencia"]

            result, _ = _worker_exec(
                execution["unit_id"],
                seed,
                strings,
                alphabet,
                self._batch_config.system.distance_method,
                self._batch_config.system.enable_distance_cache,
                execution["algorithm_params"],
                self._combination_store.work_id,
                self._execution_controller.internal_jobs,
                execution["algorithm_name"],
                worker_config.get("cpu"),
            )
            results.append({execution["unit_id"]: result})
            
        return results
