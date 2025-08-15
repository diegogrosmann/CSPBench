"""PipelineRunner com suporte a controle (pause/cancel) via WorkManager."""

from __future__ import annotations

import time
from src.domain.config import (
    CSPBenchConfig,
    ExperimentTask,
    OptimizationTask,
    SensitivityTask,
    TasksGroup,
    SyntheticDataset,
    FileDataset,
    EntrezDataset,
)
from pathlib import Path
from src.infrastructure.monitoring.monitor_interface import Monitor, NoOpMonitor
from src.infrastructure.persistence.work_state_store import WorkStateStore
from src.application.work.manager import get_work_manager
from .experiment_executor import ExperimentExecutor
from .optimization_executor import OptimizationExecutor
from .sensitivity_executor import SensitivityExecutor


class PipelineRunner:
    def __init__(self, monitor: Monitor, work_id: str, store: WorkStateStore):
        self.monitor = monitor
        self.work_id = work_id
        self.store = store

    def run(self, config: CSPBenchConfig) -> None:
        tasks_group: TasksGroup = config.tasks
        tasks = tasks_group.items
        self.monitor.pipeline_started(config.metadata.name, len(tasks))
        wm = get_work_manager()

        def check_control() -> str:
            if not self.work_id:
                return "running"
            status = wm.get_status(self.work_id)
            return status or "running"

        for task in tasks:
            if self.work_id and check_control() == "canceled":
                self.monitor.log(
                    "warning", "Pipeline canceled before task start", {"task": task.id}
                )
                break
            while self.work_id and check_control() == "paused":
                time.sleep(0.3)
            self.monitor.task_started(task.id, {"type": task.type, "name": task.name})
            for dataset in task.datasets:
                # Resolve dataset into strings using DatasetResolver
                resolver = DatasetResolver(datasets_base_path="datasets")
                try:
                    strings: list[str] = resolver.to_strings(dataset)
                except Exception as e:  # noqa: BLE001
                    self.monitor.error(
                        f"dataset:{getattr(dataset, 'id', 'unknown')}", e
                    )
                    strings = []

                # Persist in work storage when available
                try:
                    if self.store and getattr(dataset, "id", None):
                        self.store.ensure_dataset(dataset)
                        if strings:
                            self.store.add_sequences(dataset.id, strings)
                except Exception:
                    pass
                for preset in task.algorithms:
                    for alg in preset.items:
                        status = check_control()
                        if status == "canceled":
                            self.monitor.log(
                                "warning",
                                "Pipeline canceled mid-execution",
                                {"task": task.id},
                            )
                            break
                        while status == "paused":
                            time.sleep(0.3)
                            status = check_control()
                        if isinstance(task, ExperimentTask):
                            executor = ExperimentExecutor()
                        elif isinstance(task, OptimizationTask):
                            executor = OptimizationExecutor()
                        elif isinstance(task, SensitivityTask):
                            executor = SensitivityExecutor()
                        else:
                            continue
                        # Pass only the strings to the executor as requested
                        executor.run(
                            task=task,  # type: ignore[arg-type]
                            dataset_obj=strings,  # only strings
                            alg=alg,
                            resources=config.resources,
                            monitor=self.monitor,
                            global_seed=(
                                config.system.global_seed if config.system else None
                            ),
                            check_control=check_control,
                            store=self.store,
                        )
                    else:
                        continue
                    break
            self.monitor.task_finished(task.id, {"status": "ok"})
        self.monitor.pipeline_finished(True)
