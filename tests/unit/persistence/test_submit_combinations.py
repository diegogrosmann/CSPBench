import tempfile
from pathlib import Path

from src.infrastructure.persistence.work_state import WorkPersistence
from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
    WorkScopedPersistence,
)


def test_submit_combinations_inserts_all(tmp_path: Path):
    db_dir = tmp_path / "work"
    db_dir.mkdir()
    store = WorkPersistence(f"sqlite:///{db_dir / 'state.db'}")
    work_id = "W1"

    # inserir work mínimo
    class DummyWork:
        def to_dict(self):
            return {
                "id": work_id,
                "config_json": "{}",
                "status": "queued",
                "created_at": 0.0,
                "updated_at": 0.0,
                "output_path": str(db_dir),
                "error": None,
                "extra_json": "{}",
            }

    store.submit_work(DummyWork())
    scoped = WorkScopedPersistence(work_id, store)

    combinations = [
        {
            "task_id": "t1",
            "dataset_id": "d1",
            "preset_id": "p1",
            "algorithm_id": "Baseline",
            "mode": "experiment",
            "total_sequences": 2,
        },
        {
            "task_id": "t2",
            "dataset_id": "d1",
            "preset_id": "p1",
            "algorithm_id": "Baseline",
            "mode": "experiment",
            "total_sequences": 3,
        },
    ]

    inserted = scoped.submit_combinations(combinations)
    assert inserted == 2

    # Repetir (INSERT OR IGNORE) deve inserir zero novas
    again = scoped.submit_combinations(combinations)
    assert again == 0
