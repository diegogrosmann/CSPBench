"""
Unified Executor

Executes processing, optimization (with trials and repetitions) and analysis
using a single event model for display. Designed to be integrated incrementally
without breaking existing orchestrators.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple
import statistics
import time

from src.application.monitoring.progress_events import DisplayEvent, UnifiedPhase
from src.infrastructure.logging_config import get_logger
from pathlib import Path

# Optional Optuna import (only used when available and configured)
try:
    import optuna
    from optuna.pruners import MedianPruner, SuccessiveHalvingPruner
    from optuna.samplers import CmaEsSampler, RandomSampler, TPESampler
except Exception:  # pragma: no cover - optuna may be optional in some environments
    optuna = None  # type: ignore
    MedianPruner = SuccessiveHalvingPruner = TPESampler = RandomSampler = CmaEsSampler = None  # type: ignore


def _seed_from(*parts: Any, modulo: int = 2**31 - 1) -> int:
    """Deterministic seed from parts (supports None) within 32-bit range."""
    s = ":".join(str(p) for p in parts)
    return (hash(s) % modulo) + 1


@dataclass
class TrialStats:
    mean: float
    std: float
    values: List[float]
    duration_s: float


class UnifiedExecutor:
    def run_optimization(self, cfg: Dict[str, Any]) -> Dict[str, Any]:
        """Executa a fase de otimização, suportando Optuna e fallback simples."""
        opt = cfg.get("optimization", {})
        tasks = opt.get("tasks", [])
        direction_default = opt.get("direction", "minimize")
        use_optuna = optuna is not None and opt.get("use_optuna", True)
        repetitions_default = int(opt.get("repetitions", 1))
        trials_default = int(opt.get("trials", 10))
        results: List[Dict[str, Any]] = []
        # Determinar apenas as tasks que realmente têm trabalho
        effective_tasks: List[Dict[str, Any]] = []
        for t in tasks:
            datasets = t.get("datasets", []) or []
            alg_confs = t.get("algorithm_config", []) or []
            raw_params = t.get("parameters", {}) or {}
            has_work = False
            for conf_id in alg_confs:
                algo_names, _ = self._resolve_algorithms(cfg, conf_id)
                for algo_name in algo_names:
                    if self._extract_algo_param_space(raw_params, algo_name):
                        has_work = True
                        break
                if has_work:
                    break
            if has_work and datasets:
                effective_tasks.append(t)

        total_tasks = len(effective_tasks)
        for task_index, task in enumerate(effective_tasks, start=1):
            task_id = task.get("id", f"opt_task_{task_index}")
            task_name = task.get("name", task_id)
            direction = task.get("direction", direction_default)
            trials = int(task.get("trials", task.get("n_trials", trials_default)))
            repetitions = int(task.get("repetitions", repetitions_default))
            datasets = task.get("datasets", [])
            alg_confs = task.get("algorithm_config", [])
            per_task_optuna_cfg = task.get("optuna_config", {})

            for ds_idx, ds_id in enumerate(datasets, start=1):
                dataset = self._load_dataset_by_id(cfg, ds_id)
                for cfg_idx, conf_id in enumerate(alg_confs, start=1):
                    algo_names, base_params = self._resolve_algorithms(cfg, conf_id)
                    for algo_name in algo_names:
                        raw_params = task.get("parameters", {})
                        algo_param_space = self._extract_algo_param_space(raw_params, algo_name)
                        if not algo_param_space:
                            self._emit(DisplayEvent(
                                phase=UnifiedPhase.OPTIMIZATION,
                                message="algorithm skipped",
                                progress=0.0,
                                dataset_id=ds_id,
                                algorithm_name=algo_name,
                                task_id=task_id,
                                payload={
                                    "reason": "no parameters defined",
                                    "task_id": task_id,
                                    "task_name": task_name,
                                    "task_index": task_index,
                                    "total_tasks": total_tasks,
                                    "dataset_index": ds_idx,
                                    "total_datasets": len(datasets),
                                    "config_index": cfg_idx,
                                    "total_configs": len(alg_confs),
                                },
                            ))
                            continue
                        if use_optuna:
                            study = self._create_study(
                                cfg,
                                opt.get("optuna_defaults", {}),
                                per_task_optuna_cfg,
                                direction,
                                f"{task.get('study_name', task_id)}_{algo_name}_{ds_id}",
                            )
                            def objective(trial):
                                try:
                                    params = self._suggest_params_optuna(
                                        trial, algo_param_space, base_params.get(algo_name, {})
                                    )
                                    values: List[float] = []
                                    t_start = time.time()
                                    trial_no = trial.number + 1
                                    self._emit(DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="trial start",
                                        progress=0.0,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        trial_no=trial_no,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "params": params,
                                            "direction": direction,
                                            "total_trials": trials,
                                            "repetitions": repetitions,
                                            "dataset_index": ds_idx,
                                            "total_datasets": len(datasets),
                                            "config_index": cfg_idx,
                                            "total_configs": len(alg_confs),
                                        },
                                    ))
                                    for rep in range(1, repetitions + 1):
                                        seed = self._derive_seed(cfg, ds_id, algo_name, trial_no, rep)
                                        value, meta = self._run_algorithm(algo_name, dataset, {**params, "seed": seed})
                                        values.append(float(value))
                                        try:
                                            trial.report(float(statistics.fmean(values)), step=rep)
                                        except Exception:
                                            pass
                                        self._emit(DisplayEvent(
                                            phase=UnifiedPhase.OPTIMIZATION,
                                            message="trial rep",
                                            progress=rep / repetitions,
                                            dataset_id=ds_id,
                                            algorithm_name=algo_name,
                                            task_id=task_id,
                                            trial_no=trial_no,
                                            rep_idx=rep,
                                            payload={
                                                "task_id": task_id,
                                                "task_name": task_name,
                                                "task_index": task_index,
                                                "total_tasks": total_tasks,
                                                "partial_value": value,
                                                "seed": seed,
                                                "repetitions": repetitions,
                                                "dataset_index": ds_idx,
                                                "total_datasets": len(datasets),
                                                "config_index": cfg_idx,
                                                "total_configs": len(alg_confs),
                                            },
                                        ))
                                    mean_v = float(statistics.fmean(values)) if values else float("inf")
                                    std_v = float(statistics.pstdev(values)) if len(values) > 1 else 0.0
                                    t_dur = time.time() - t_start
                                    try:
                                        trial.set_user_attr("values", values)
                                        trial.set_user_attr("std", std_v)
                                    except Exception:
                                        pass
                                    self._emit(DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="trial end",
                                        progress=1.0,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        trial_no=trial_no,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "value_mean": mean_v,
                                            "value_std": std_v,
                                            "duration_s": t_dur,
                                            "direction": direction,
                                            "total_trials": trials,
                                            "repetitions": repetitions,
                                            "dataset_index": ds_idx,
                                            "total_datasets": len(datasets),
                                            "config_index": cfg_idx,
                                            "total_configs": len(alg_confs),
                                        },
                                    ))
                                    return mean_v
                                except KeyboardInterrupt:
                                    self._emit(DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="optimization cancelled",
                                        progress=0.0,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "cancelled": True,
                                        },
                                    ))
                                    raise
                            parallel_cfg = cfg.get("resources", {}).get("parallel", {}) or {}
                            n_jobs = parallel_cfg.get("n_jobs") or parallel_cfg.get("max_workers") or 1
                            timeout = int(task.get("timeout_per_trial", 300)) * trials
                            try:
                                study.optimize(
                                    objective,
                                    n_trials=trials,
                                    timeout=timeout,
                                    n_jobs=n_jobs,
                                    show_progress_bar=False,
                                )
                            except KeyboardInterrupt:
                                self._emit(DisplayEvent(
                                    phase=UnifiedPhase.OPTIMIZATION,
                                    message="optimization cancelled",
                                    progress=0.0,
                                    dataset_id=ds_id,
                                    algorithm_name=algo_name,
                                    task_id=task_id,
                                    payload={
                                        "task_id": task_id,
                                        "task_name": task_name,
                                        "task_index": task_index,
                                        "total_tasks": total_tasks,
                                        "cancelled": True,
                                    },
                                ))
                                break
                            study_res = {
                                "task_id": task_id,
                                "dataset": ds_id,
                                "algorithm": algo_name,
                                "direction": study.direction.name if hasattr(study, "direction") else direction,
                                "best_value": getattr(study, "best_value", None),
                                "best_params": getattr(study, "best_params", None),
                                "n_trials": len(study.trials),
                                "trials": [
                                    {
                                        "number": t.number,
                                        "value": t.value,
                                        "params": t.params,
                                        "state": t.state.name if hasattr(t.state, "name") else str(t.state),
                                        "duration": self._safe_trial_duration(t),
                                        "user_attrs": getattr(t, "user_attrs", {}),
                                    }
                                    for t in study.trials
                                ],
                            }
                            study_res["status"] = "success"
                            results.append(study_res)
                            self._optuna_contexts.append({
                                "study": study,
                                "results": study_res,
                                "algo": algo_name,
                                "dataset": ds_id,
                                "task_id": task_id,
                            })
                        else:
                            best_value: Optional[float] = None
                            best_params: Optional[Dict[str, Any]] = None
                            trial_history: List[Dict[str, Any]] = []
                            for trial_no in range(1, trials + 1):
                                params = self._suggest_params(algo_param_space, base_params.get(algo_name, {}), trial_no)
                                values: List[float] = []
                                t_start = time.time()
                                self._emit(DisplayEvent(
                                    phase=UnifiedPhase.OPTIMIZATION,
                                    message="trial start",
                                    progress=0.0,
                                    dataset_id=ds_id,
                                    algorithm_name=algo_name,
                                    task_id=task_id,
                                    trial_no=trial_no,
                                    payload={
                                        "task_id": task_id,
                                        "task_name": task_name,
                                        "task_index": task_index,
                                        "total_tasks": total_tasks,
                                        "params": params,
                                        "direction": direction,
                                        "total_trials": trials,
                                        "repetitions": repetitions,
                                        "dataset_index": ds_idx,
                                        "total_datasets": len(datasets),
                                        "config_index": cfg_idx,
                                        "total_configs": len(alg_confs),
                                    },
                                ))
                                for rep in range(1, repetitions + 1):
                                    seed = self._derive_seed(cfg, ds_id, algo_name, trial_no, rep)
                                    value, meta = self._run_algorithm(algo_name, dataset, {**params, "seed": seed})
                                    values.append(float(value))
                                    self._emit(DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="trial rep",
                                        progress=rep / repetitions,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        trial_no=trial_no,
                                        rep_idx=rep,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "partial_value": value,
                                            "seed": seed,
                                            "repetitions": repetitions,
                                            "dataset_index": ds_idx,
                                            "total_datasets": len(datasets),
                                            "config_index": cfg_idx,
                                            "total_configs": len(alg_confs),
                                        },
                                    ))
                                mean_v = float(statistics.fmean(values)) if values else float("inf")
                                std_v = float(statistics.pstdev(values)) if len(values) > 1 else 0.0
                                t_dur = time.time() - t_start
                                improved = (
                                    best_value is None
                                    or (direction == "minimize" and mean_v < best_value)
                                    or (direction == "maximize" and mean_v > best_value)
                                )
                                if improved:
                                    best_value, best_params = mean_v, params
                                trial_info = {
                                    "trial": trial_no,
                                    "params": params,
                                    "mean": mean_v,
                                    "std": std_v,
                                    "values": values,
                                    "duration_s": t_dur,
                                }
                                trial_history.append(trial_info)
                                self._emit(DisplayEvent(
                                    phase=UnifiedPhase.OPTIMIZATION,
                                    message="trial end",
                                    progress=1.0,
                                    dataset_id=ds_id,
                                    algorithm_name=algo_name,
                                    task_id=task_id,
                                    trial_no=trial_no,
                                    payload={
                                        "task_id": task_id,
                                        "task_name": task_name,
                                        "task_index": task_index,
                                        "total_tasks": total_tasks,
                                        "value_mean": mean_v,
                                        "value_std": std_v,
                                        "best_so_far": best_value,
                                        "duration_s": t_dur,
                                        "direction": direction,
                                        "total_trials": trials,
                                        "repetitions": repetitions,
                                        "dataset_index": ds_idx,
                                        "total_datasets": len(datasets),
                                        "config_index": cfg_idx,
                                        "total_configs": len(alg_confs),
                                    },
                                ))
                            results.append({
                                "task_id": task_id,
                                "dataset": ds_id,
                                "algorithm": algo_name,
                                "direction": direction,
                                "best_value": best_value,
                                "best_params": best_params,
                                "trials": trial_history,
                                "status": "success",
                            })
        return {"optimizations": results}
    """Single orchestrator for processing, optimization and analysis phases."""

    def __init__(self, algorithm_registry, dataset_repository, monitoring_service=None):
        self._algos = algorithm_registry
        self._datasets = dataset_repository
        self._monitor = monitoring_service
        self._logger = get_logger(__name__)
        # Store optuna studies for analysis/report generation
        self._optuna_contexts: List[Dict[str, Any]] = []
        # Reduce Optuna verbosity in console and suppress experimental warnings
        if optuna is not None:
            try:
                import warnings
                from optuna.exceptions import ExperimentalWarning as OptunaExperimentalWarning
                from optuna import logging as optuna_logging

                # Only warnings and above from Optuna
                optuna_logging.set_verbosity(optuna_logging.WARNING)
                # Hide experimental feature warnings from Optuna
                warnings.filterwarnings("ignore", category=OptunaExperimentalWarning)
            except Exception:
                # Best-effort; ignore if not available
                pass

    # Public entrypoint
    def run_pipeline(self, batch_config: Dict[str, Any]) -> Dict[str, Any]:
        results: Dict[str, Any] = {}

        task_type = batch_config.get("task", {}).get("type", "execution")
        if task_type == "execution":
            results["processing"] = self.run_processing(batch_config)
        elif task_type == "optimization":
            results["optimization"] = self.run_optimization(batch_config)
            # Optional post-analysis hook (reuse existing report generator later)
            results["analysis"] = self.run_analysis(batch_config, results)
        elif task_type == "sensitivity":
            # Future: unify sensitivity here as well
            pass
        return results

    # Processing phase
    def run_processing(self, cfg: Dict[str, Any]) -> Dict[str, Any]:
        execs = cfg.get("execution", {}).get("executions", [])
        # ...existing code...

        # Determine only the tasks that actually have work (at least one algo with params)
        effective_tasks: List[Dict[str, Any]] = []
        for t in tasks:
            datasets = t.get("datasets", []) or []
            alg_confs = t.get("algorithm_config", []) or []
            raw_params = t.get("parameters", {}) or {}
            has_work = False
            for conf_id in alg_confs:
                algo_names, _ = self._resolve_algorithms(cfg, conf_id)
                for algo_name in algo_names:
                    if self._extract_algo_param_space(raw_params, algo_name):
                        has_work = True
                        break
                if has_work:
                    break
            if has_work and datasets:
                effective_tasks.append(t)

        total_tasks = len(effective_tasks)
        try:
            import os as _os
                # Debug removido
        except Exception:
            pass
        for task_index, task in enumerate(effective_tasks, start=1):
            task_id = task.get("id", "opt_task")
            task_name = task.get("name", task_id)
            direction = task.get("direction", direction_default)
            trials = int(task.get("trials", task.get("n_trials", 10)))
            repetitions = int(task.get("repetitions", 1))
            datasets = task.get("datasets", [])
            alg_confs = task.get("algorithm_config", [])
            per_task_optuna_cfg = task.get("optuna_config", {})

            for ds_idx, ds_id in enumerate(datasets, start=1):
                dataset = self._load_dataset_by_id(cfg, ds_id)
                for cfg_idx, conf_id in enumerate(alg_confs, start=1):
                    algo_names, base_params = self._resolve_algorithms(cfg, conf_id)
                    for algo_name in algo_names:
                        # Extract parameter space for this algorithm only
                        raw_params = task.get("parameters", {})
                        algo_param_space = self._extract_algo_param_space(raw_params, algo_name)

                        # Skip algorithms with no parameter space defined in this task
                        if not algo_param_space:
                            # Emit a low-verbosity event for monitoring callbacks
                            self._emit(
                                DisplayEvent(
                                    phase=UnifiedPhase.OPTIMIZATION,
                                    message="algorithm skipped",
                                    progress=0.0,
                                    dataset_id=ds_id,
                                    algorithm_name=algo_name,
                                    task_id=task_id,
                                    payload={
                                        "reason": "no parameters defined",
                                        "task_id": task_id,
                                        "task_name": task_name,
                                        "task_index": task_index,
                                        "total_tasks": total_tasks,
                                        "dataset_index": ds_idx,
                                        "total_datasets": len(datasets),
                                        "config_index": cfg_idx,
                                        "total_configs": len(alg_confs),
                                    },
                                )
                            )
                            continue

                        # Optuna-backed optimization
                        if use_optuna:
                            study = self._create_study(
                                cfg,
                                opt.get("optuna_defaults", {}),
                                per_task_optuna_cfg,
                                direction,
                                f"{task.get('study_name', task_id)}_{algo_name}_{ds_id}",
                            )

                            def objective(trial):
                                # Params via optuna
                                params = self._suggest_params_optuna(
                                    trial, algo_param_space, base_params.get(algo_name, {})
                                )
                                values: List[float] = []
                                t_start = time.time()
                                trial_no = trial.number + 1

                            # Debug removido

                                self._emit(
                                    DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="trial start",
                                        progress=0.0,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        trial_no=trial_no,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "params": params,
                                            "direction": direction,
                                            "total_trials": trials,
                                            "repetitions": repetitions,
                                            "dataset_index": ds_idx,
                                            "total_datasets": len(datasets),
                                            "config_index": cfg_idx,
                                            "total_configs": len(alg_confs),
                                        },
                                    )
                                )

                                for rep in range(1, repetitions + 1):
                                    seed = self._derive_seed(
                                        cfg, ds_id, algo_name, trial_no, rep
                                    )
                                    value, meta = self._run_algorithm(
                                        algo_name, dataset, {**params, "seed": seed}
                                    )
                                    values.append(float(value))
                                    # Optional: report intermediate
                                    try:
                                        trial.report(
                                            float(statistics.fmean(values)), step=rep
                                        )
                                    except Exception:
                                        pass
                                    self._emit(
                                        DisplayEvent(
                                            phase=UnifiedPhase.OPTIMIZATION,
                                            message="trial rep",
                                            progress=rep / repetitions,
                                            dataset_id=ds_id,
                                            algorithm_name=algo_name,
                                            task_id=task_id,
                                            trial_no=trial_no,
                                            rep_idx=rep,
                                            payload={
                                                "task_id": task_id,
                                                "task_name": task_name,
                                                "task_index": task_index,
                                                "total_tasks": total_tasks,
                                                "partial_value": value,
                                                "seed": seed,
                                                "repetitions": repetitions,
                                                "dataset_index": ds_idx,
                                                "total_datasets": len(datasets),
                                                "config_index": cfg_idx,
                                                "total_configs": len(alg_confs),
                                            },
                                        )
                                    )

                                mean_v = (
                                    float(statistics.fmean(values)) if values else float("inf")
                                )
                                std_v = (
                                    float(statistics.pstdev(values)) if len(values) > 1 else 0.0
                                )
                                t_dur = time.time() - t_start

                                # Set attrs for richer reporting
                                try:
                                    trial.set_user_attr("values", values)
                                    trial.set_user_attr("std", std_v)
                                except Exception:
                                    pass

                                self._emit(
                                    DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="trial end",
                                        progress=1.0,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        trial_no=trial_no,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "value_mean": mean_v,
                                            "value_std": std_v,
                                            "duration_s": t_dur,
                                            "direction": direction,
                                            "total_trials": trials,
                                            "repetitions": repetitions,
                                            "dataset_index": ds_idx,
                                            "total_datasets": len(datasets),
                                            "config_index": cfg_idx,
                                            "total_configs": len(alg_confs),
                                        },
                                    )
                                )
                                return mean_v

                            # Optimize
                            parallel_cfg = (
                                cfg.get("resources", {}).get("parallel", {}) or {}
                            )
                            n_jobs = (
                                parallel_cfg.get("n_jobs")
                                or parallel_cfg.get("max_workers")
                                or 1
                            )
                            timeout = int(task.get("timeout_per_trial", 300)) * trials
                            study.optimize(
                                objective,
                                n_trials=trials,
                                timeout=timeout,
                                n_jobs=n_jobs,
                                show_progress_bar=False,
                            )

                            # Consolidate results from study
                            study_res = {
                                "task_id": task_id,
                                "dataset": ds_id,
                                "algorithm": algo_name,
                                "direction": study.direction.name
                                if hasattr(study, "direction")
                                else direction,
                                "best_value": getattr(study, "best_value", None),
                                "best_params": getattr(study, "best_params", None),
                                "n_trials": len(study.trials),
                                "trials": [
                                    {
                                        "number": t.number,
                                        "value": t.value,
                                        "params": t.params,
                                        "state": t.state.name
                                        if hasattr(t.state, "name")
                                        else str(t.state),
                                        "duration": self._safe_trial_duration(t),
                                        "user_attrs": getattr(t, "user_attrs", {}),
                                    }
                                    for t in study.trials
                                ],
                            }
                            study_res["status"] = "success"
                            results.append(study_res)
                            # Store for analysis phase (reports)
                            self._optuna_contexts.append(
                                {
                                    "study": study,
                                    "results": study_res,
                                    "algo": algo_name,
                                    "dataset": ds_id,
                                    "task_id": task_id,
                                }
                            )

                        else:
                            # Fallback: simple mid-point sampler (no Optuna available)
                            best_value: Optional[float] = None
                            best_params: Optional[Dict[str, Any]] = None
                            trial_history: List[Dict[str, Any]] = []

                            for trial_no in range(1, trials + 1):
                                params = self._suggest_params(
                                    algo_param_space, base_params.get(algo_name, {}), trial_no
                                )
                                values: List[float] = []
                                t_start = time.time()

                                self._emit(
                                    DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="trial start",
                                        progress=0.0,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        trial_no=trial_no,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "params": params,
                                            "direction": direction,
                                            "total_trials": trials,
                                            "repetitions": repetitions,
                                            "dataset_index": ds_idx,
                                            "total_datasets": len(datasets),
                                            "config_index": cfg_idx,
                                            "total_configs": len(alg_confs),
                                        },
                                    )
                                )

                                for rep in range(1, repetitions + 1):
                                    seed = self._derive_seed(
                                        cfg, ds_id, algo_name, trial_no, rep
                                    )
                                    value, meta = self._run_algorithm(
                                        algo_name, dataset, {**params, "seed": seed}
                                    )
                                    values.append(float(value))

                                    self._emit(
                                        DisplayEvent(
                                            phase=UnifiedPhase.OPTIMIZATION,
                                            message="trial rep",
                                            progress=rep / repetitions,
                                            dataset_id=ds_id,
                                            algorithm_name=algo_name,
                                            task_id=task_id,
                                            trial_no=trial_no,
                                            rep_idx=rep,
                                            payload={
                                                "task_id": task_id,
                                                "task_name": task_name,
                                                "task_index": task_index,
                                                "total_tasks": total_tasks,
                                                "partial_value": value,
                                                "seed": seed,
                                                "repetitions": repetitions,
                                                "dataset_index": ds_idx,
                                                "total_datasets": len(datasets),
                                                "config_index": cfg_idx,
                                                "total_configs": len(alg_confs),
                                            },
                                        )
                                    )

                                mean_v = (
                                    float(statistics.fmean(values)) if values else float("inf")
                                )
                                std_v = (
                                    float(statistics.pstdev(values)) if len(values) > 1 else 0.0
                                )
                                t_dur = time.time() - t_start
                                improved = (
                                    best_value is None
                                    or (direction == "minimize" and mean_v < best_value)
                                    or (direction == "maximize" and mean_v > best_value)
                                )
                                if improved:
                                    best_value, best_params = mean_v, params

                                trial_info = {
                                    "trial": trial_no,
                                    "params": params,
                                    "mean": mean_v,
                                    "std": std_v,
                                    "values": values,
                                    "duration_s": t_dur,
                                }
                                trial_history.append(trial_info)

                                self._emit(
                                    DisplayEvent(
                                        phase=UnifiedPhase.OPTIMIZATION,
                                        message="trial end",
                                        progress=1.0,
                                        dataset_id=ds_id,
                                        algorithm_name=algo_name,
                                        task_id=task_id,
                                        trial_no=trial_no,
                                        payload={
                                            "task_id": task_id,
                                            "task_name": task_name,
                                            "task_index": task_index,
                                            "total_tasks": total_tasks,
                                            "value_mean": mean_v,
                                            "value_std": std_v,
                                            "best_so_far": best_value,
                                            "duration_s": t_dur,
                                            "direction": direction,
                                            "total_trials": trials,
                                            "repetitions": repetitions,
                                            "dataset_index": ds_idx,
                                            "total_datasets": len(datasets),
                                            "config_index": cfg_idx,
                                            "total_configs": len(alg_confs),
                                        },
                                    )
                                )

                            results.append(
                                {
                                    "task_id": task_id,
                                    "dataset": ds_id,
                                    "algorithm": algo_name,
                                    "direction": direction,
                                    "best_value": best_value,
                                    "best_params": best_params,
                                    "trials": trial_history,
                                    "status": "success",
                                }
                            )

        return {"optimizations": results}

    # Analysis phase (placeholder: can call existing report generator)
    def run_analysis(self, cfg: Dict[str, Any], results: Dict[str, Any]) -> Dict[str, Any]:
        # Generate optimization reports if we have optuna studies
        start = time.time()
        artifacts: List[Dict[str, Any]] = []
        self._emit(DisplayEvent(
            phase=UnifiedPhase.ANALYSIS,
            message="analysis start",
            progress=0.0,
            payload={"artifacts": []},
        ))

        if self._optuna_contexts:
            try:
                from src.infrastructure.orchestrators.optimization_report_generator import (
                    OptimizationReportGenerator,
                )
                # Compute destination using SessionManager and settings.yaml if available
                destination = self._resolve_output_destination()
                org = OptimizationReportGenerator(destination=destination, config=cfg)
                for ctx in self._optuna_contexts:
                    study = ctx["study"]
                    res = ctx["results"]
                    org.generate_all_reports(study, res)
                    artifacts.append({
                        "dataset": ctx["dataset"],
                        "algorithm": ctx["algo"],
                        "task_id": ctx["task_id"],
                        "reports_path": destination,
                    })
            except Exception as e:
                self._logger.warning(f"Falha ao gerar relatórios de otimização: {e}")

        dur = time.time() - start
        self._emit(DisplayEvent(
            phase=UnifiedPhase.ANALYSIS,
            message="analysis end",
            progress=1.0,
            payload={"duration_s": dur, "artifacts_count": len(artifacts)},
        ))
        return {"artifacts": artifacts}

    # Helpers
    def _emit(self, event: DisplayEvent) -> None:
        if self._monitor:
            try:
                self._monitor.emit_event(event)
            except Exception:
                pass

    def _safe_trial_duration(self, t) -> Optional[float]:
        """Return trial duration in seconds if available and valid."""
        try:
            dur = getattr(t, "duration", None)
            if dur is None:
                return None
            return float(dur.total_seconds())
        except Exception:
            return None

    def _load_dataset_by_id(self, cfg: Dict[str, Any], dataset_id: str):
        # Find dataset config by id in batch cfg Section 2
        dcfgs = cfg.get("datasets", []) or []
        dcfg = next((d for d in dcfgs if d.get("id") == dataset_id), None)
        if not dcfg:
            # Fallback: allow repository to try loading by raw id
            return self._datasets.load(dataset_id)

        dtype = dcfg.get("type")
        params = dcfg.get("parameters", {})
        # File dataset: use filename relative to repository base
        if dtype == "file":
            fname = params.get("filename") or dataset_id
            return self._datasets.load(fname)
        elif dtype == "synthetic":
            # Build synthetic dataset using domain generator
            from src.domain import SyntheticDatasetGenerator, Dataset
            n = int(params.get("n", 10))
            L = int(params.get("L", 20))
            alphabet = params.get("alphabet", "ACGT")
            seed = params.get("seed")
            noise = float(params.get("noise", 0.0))
            gen = SyntheticDatasetGenerator()
            if noise and noise > 0:
                import random
                rng = random.Random(seed)
                center = "".join(rng.choice(alphabet) for _ in range(L))
                return gen.generate_from_center(center=center, n=n, noise_rate=noise, alphabet=alphabet, seed=seed)
            else:
                return gen.generate_random(n=n, length=L, alphabet=alphabet, seed=seed)
        elif dtype == "entrez":
            # Fetch sequences from NCBI Entrez when available
            try:
                from src.infrastructure.persistence.entrez_dataset_repository import (
                    NCBIEntrezDatasetRepository,
                )
                from src.domain import Dataset

                entrez_repo = NCBIEntrezDatasetRepository()
                if not entrez_repo.is_available():
                    raise RuntimeError(
                        "Entrez repository not available. Check NCBI_EMAIL and network."
                    )

                query = params.get("query")
                if not query:
                    raise ValueError("Field 'query' is required for entrez dataset")
                db = params.get("db", "nucleotide")
                retmax = params.get("retmax", 20)

                # Avoid passing duplicate keys present in params
                extra = {k: v for k, v in (params or {}).items() if k not in {"query", "db", "retmax"}}
                sequences, meta = entrez_repo.fetch_dataset(
                    query=query, db=db, retmax=retmax, **extra
                )
                if not sequences:
                    raise RuntimeError("No sequences returned from Entrez query")

                dataset = Dataset(
                    sequences=sequences,
                    metadata={
                        "type": "entrez",
                        "query": query,
                        "db": db,
                        "retmax": retmax,
                        "n_obtained": len(sequences),
                        "L": len(sequences[0]) if sequences else 0,
                        **(meta or {}),
                    },
                )
                return dataset
            except Exception as e:
                # Fallback with clear error
                raise RuntimeError(f"Failed to load Entrez dataset '{dataset_id}': {e}")
        else:
            # Unknown type, fallback to repository
            return self._datasets.load(dataset_id)

    def _resolve_algorithms(self, cfg: Dict[str, Any], config_id: str) -> Tuple[List[str], Dict[str, Dict[str, Any]]]:
        alg_cfgs = cfg.get("algorithms", [])
        for entry in alg_cfgs:
            if entry.get("id") == config_id:
                names = entry.get("algorithms", [])
                params = entry.get("algorithm_params", {})
                return names, params
        # Fallback: assume config_id is directly the algorithm name
        return [config_id], {}

    def _create_study(
        self,
        cfg: Dict[str, Any],
        global_optuna: Dict[str, Any],
        per_task: Dict[str, Any],
        direction: str,
        study_name: str,
    ):
        if optuna is None:  # safety
            raise RuntimeError("Optuna não disponível")

        # Merge configs: per_task overrides global
        ocfg = {**(global_optuna or {}), **(per_task or {})}
        sampler_name = ocfg.get("sampler", "TPESampler")
        pruner_name = ocfg.get("pruner", "MedianPruner")

        # Sampler
        sampler = None
        try:
            if sampler_name == "TPESampler" and TPESampler is not None:
                sampler = TPESampler(
                    n_startup_trials=ocfg.get("n_startup_trials", ocfg.get("startup_trials", 10)),
                    n_ei_candidates=ocfg.get("n_ei_candidates", ocfg.get("ei_candidates", 24)),
                    multivariate=ocfg.get("multivariate", False),
                )
            elif sampler_name == "RandomSampler" and RandomSampler is not None:
                sampler = RandomSampler(seed=ocfg.get("seed"))
            elif sampler_name == "CmaEsSampler" and CmaEsSampler is not None:
                sampler = CmaEsSampler(
                    n_startup_trials=ocfg.get("n_startup_trials", 1),
                    restart_strategy=ocfg.get("restart_strategy", "ipop"),
                )
        except Exception:
            sampler = None

        # Pruner
        pruner = None
        try:
            if pruner_name == "MedianPruner" and MedianPruner is not None:
                pruner = MedianPruner(
                    n_startup_trials=ocfg.get("n_startup_trials", ocfg.get("startup_trials", 5)),
                    n_warmup_steps=ocfg.get("n_warmup_steps", ocfg.get("warmup_steps", 10)),
                    interval_steps=ocfg.get("interval_steps", 5),
                )
            elif pruner_name == "SuccessiveHalvingPruner" and SuccessiveHalvingPruner is not None:
                pruner = SuccessiveHalvingPruner(
                    min_resource=ocfg.get("min_resource", 1),
                    reduction_factor=ocfg.get("reduction_factor", 4),
                    min_early_stopping_rate=ocfg.get("min_early_stopping_rate", 0),
                )
        except Exception:
            pruner = None

        storage_url = ocfg.get("storage")
        try:
            study = optuna.create_study(
                study_name=study_name,
                direction=direction,
                sampler=sampler,
                pruner=pruner,
                storage=storage_url,
                load_if_exists=True,
            )
        except Exception:
            # Fallback to in-memory
            study = optuna.create_study(
                study_name=study_name,
                direction=direction,
                sampler=sampler,
                pruner=pruner,
            )
        return study

    def _suggest_params_optuna(self, trial, space: Dict[str, Any], base: Dict[str, Any]) -> Dict[str, Any]:
        params = dict(base)
        for name, spec in (space or {}).items():
            t = spec.get("type", "uniform")
            if t == "categorical":
                choices = spec.get("choices") or spec.get("values") or []
                params[name] = trial.suggest_categorical(name, choices)
            elif t in ("int", "integer"):
                low = int(spec.get("low", 0))
                high = int(spec.get("high", 1))
                step = spec.get("step")
                if step is None or int(step) <= 0:
                    params[name] = trial.suggest_int(name, low, high)
                else:
                    params[name] = trial.suggest_int(name, low, high, step=int(step))
            elif t in ("loguniform", "log_uniform"):
                low, high = float(spec.get("low", 1e-3)), float(spec.get("high", 1.0))
                params[name] = trial.suggest_float(name, low, high, log=True)
            else:  # float/uniform
                low, high = float(spec.get("low", 0.0)), float(spec.get("high", 1.0))
                step = spec.get("step")
                params[name] = trial.suggest_float(name, low, high, step=step)
        return params

    def _extract_algo_param_space(self, task_params: Dict[str, Any], algo_name: str) -> Dict[str, Any]:
        """Return parameter space dict for given algorithm from task parameters.

        Accepts both shapes:
        - {algo_name: {param: spec, ...}, other_algo: {...}}
        - {param: spec, ...} (already scoped for this algorithm)
        """
        if not isinstance(task_params, dict):
            return {}
        if algo_name in task_params and isinstance(task_params[algo_name], dict):
            return task_params[algo_name]
        # If looks like already a parameter space (values are spec dicts with 'type' or bounds)
        if any(
            isinstance(v, dict) and ("type" in v or "low" in v or "high" in v)
            for v in task_params.values()
        ):
            return task_params
        return {}

    def _resolve_output_destination(self) -> str:
        """Resolve destination folder consistent with session manager and output config."""
        try:
            from src.infrastructure.session_manager import SessionManager
            import yaml
            settings_path = Path("config/settings.yaml")
            if settings_path.exists():
                with open(settings_path, encoding="utf-8") as f:
                    settings = yaml.safe_load(f)
            else:
                settings = {}
            session_manager = SessionManager(settings)
            # If no session exists, create it (keeps consistent behavior)
            if not session_manager.get_session_folder():
                session_manager.create_session()
            # Determine base dir from unified output config or fallback
            output_cfg = settings.get("output", {})
            if output_cfg:
                base_directory = output_cfg.get("base_directory", "outputs/{session}")
                base_result_dir = base_directory.replace("{session}", "") if "{session}" in base_directory else base_directory
            else:
                base_result_dir = (
                    settings.get("infrastructure", {}).get("result", {}).get("base_result_dir", "./outputs/results")
                )
            session_folder = session_manager.get_session_folder() or "default_session"
            destination = str(Path(base_result_dir) / session_folder)
            Path(destination).mkdir(parents=True, exist_ok=True)
            return destination
        except Exception:
            # Fallback to outputs directory
            destination = "outputs"
            Path(destination).mkdir(parents=True, exist_ok=True)
            return destination

    def _suggest_params(self, space: Dict[str, Any], base: Dict[str, Any], trial_no: int) -> Dict[str, Any]:
        # Minimal sampler: if space empty, return base; otherwise choose midpoints/deterministic choices
        params = dict(base)
        for name, spec in (space or {}).items():
            t = spec.get("type", "uniform")
            if t == "categorical":
                choices = spec.get("choices") or spec.get("values") or []
                if choices:
                    idx = (trial_no - 1) % len(choices)
                    params[name] = choices[idx]
            elif t in ("int", "integer"):
                low, high = spec.get("low", 0), spec.get("high", 1)
                step = max(1, int(spec.get("step", 1)))
                mid = low + ((high - low) // (2 * step)) * step
                params[name] = int(mid)
            else:
                low, high = spec.get("low", 0.0), spec.get("high", 1.0)
                params[name] = float((low + high) / 2.0)
        return params

    def _derive_seed(self, cfg: Dict[str, Any], dataset_id: str, algorithm_name: str, trial_no: int, rep_idx: int) -> int:
        global_seed = (cfg.get("system", {}) or {}).get("global_seed")
        return _seed_from(global_seed or 0, dataset_id, algorithm_name, trial_no, rep_idx)

    def _run_algorithm(self, algorithm_name: str, dataset, params: Dict[str, Any]) -> Tuple[float, Dict[str, Any]]:
        # Instantiate algorithm from registry and run
        algo_cls = self._algos.get_algorithm(algorithm_name)
        instance = algo_cls(strings=dataset.sequences, alphabet=dataset.alphabet, **params)
        start = time.time()
        result = instance.run()
        dur = time.time() - start
        if isinstance(result, tuple) and len(result) >= 3:
            best_string, max_distance, metadata = result[:3]
            value = float(max_distance)
            meta = {"best_string": str(best_string), "metadata": metadata, "time_s": dur}
        else:
            value = 0.0
            meta = {"raw_result": str(result), "time_s": dur}
        return value, meta
