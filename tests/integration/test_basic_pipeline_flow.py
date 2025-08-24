import time
from pathlib import Path
from src.application.services.execution_manager import ExecutionManager
from src.domain.config import (
    CSPBenchConfig,
    MetadataConfig,
    SyntheticDatasetConfig,
    AlgorithmsPresetConfig,
    AlgParamsConfig,
    ExperimentTasksConfig,
    ExperimentTaskConfig,
    TasksGroupConfig,
    OutputConfig,
    ResultsConfig,
    ResultsFormats,  # usar classe existente
    ResourcesConfig,
    CPUConfig,
    MemoryConfig,
    TimeoutsConfig,
    SystemConfig,
)
from src.domain.status import BaseStatus


def build_minimal_config():
    metadata = MetadataConfig(
        name="basic-test",
        version="1.0",
        author="test",
        description="integration basic",
        creation_date="2025-01-01",
        tags=["test"],
    )
    datasets = {
        "ds1": SyntheticDatasetConfig(
            id="ds1",
            type="synthetic",
            name="ds1",
            alphabet="ACGT",
            n=2,
            L=4,
        )
    }
    algorithms = {
        "preset1": AlgorithmsPresetConfig(
            id="preset1",
            name="preset1",
            description="preset",
            items=[AlgParamsConfig(name="Baseline", params={})],
        )
    }
    tasks = ExperimentTasksConfig(
        type="experiment",
        items=[
            ExperimentTaskConfig(
                id="task1",
                name="task1",
                type="experiment",
                datasets=["ds1"],
                algorithms=["preset1"],
                repetitions=2,
            )
        ],
    )
    tasks_group = TasksGroupConfig(
        type="experiment", items=tasks.items
    )  # wrapper compat

    # ResultsFormats exige flags; definir defaults mínimos
    output = OutputConfig(
        logging=False,
        results=ResultsConfig(
            formats=ResultsFormats(csv=False, json=False, parquet=False, pickle=False),
            partial_results=False,
        ),
    )
    resources = ResourcesConfig(
        cpu=CPUConfig(exclusive_cores=False, max_workers=1, internal_jobs=1),
        memory=MemoryConfig(max_memory_gb=1.0),
        timeouts=TimeoutsConfig(timeout_per_item=10, timeout_total_batch=60),
    )
    system = SystemConfig(
        global_seed=123, distance_method="hamming", enable_distance_cache=False
    )

    return CSPBenchConfig(
        metadata=metadata,
        datasets=datasets,
        algorithms=algorithms,
        tasks=tasks_group,
        output=output,
        resources=resources,
        system=system,
    )


def test_basic_pipeline_flow(tmp_path, monkeypatch):
    # Definir diretório de saída
    monkeypatch.setenv("OUTPUT_BASE_DIRECTORY", str(tmp_path))

    config = build_minimal_config()
    manager = ExecutionManager()
    work_id = manager.execute(config=config)

    # Esperar término (timeout defensivo)
    start = time.time()
    final_status = None
    while time.time() - start < 30:  # aumentar um pouco o timeout
        # Carrega estado do work service interno
        status = manager._work_service.get_status(work_id)  # type: ignore[attr-defined]
        if status in (
            BaseStatus.COMPLETED.value,
            BaseStatus.FAILED.value,
            BaseStatus.ERROR.value,
        ):
            final_status = status
            break
        time.sleep(0.5)

    assert final_status == BaseStatus.COMPLETED.value

    # Verifica que o banco foi criado e contém combinações
    state_db = Path(tmp_path) / work_id / "state.db"
    assert state_db.exists()

    # Conecta leitura rápida
    import sqlite3

    conn = sqlite3.connect(state_db)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM combinations")
    combos = cur.fetchone()[0]
    assert combos == 1
    cur.execute("SELECT COUNT(*) FROM executions")
    execs = cur.fetchone()[0]
    assert execs == 2  # repetições
    conn.close()
