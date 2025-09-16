import json
import os
import sqlite3
from pathlib import Path
from unittest.mock import Mock

from src.infrastructure.export.finalization_service import (
    FinalizationConfig,
    FinalizationService,
)


def create_mock_work_store(db_path: Path):
    """Create a mock work_store that provides access to the SQLite database."""
    # Keep a persistent connection to prevent "Cannot operate on a closed database" error
    persistent_conn = sqlite3.connect(db_path)

    def get_connection():
        return persistent_conn

    mock_engine = Mock()
    mock_engine.raw_connection = get_connection

    def get_work_export_data(work_id: str):
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT * FROM work WHERE id = ?", (work_id,))
        row = cur.fetchone()
        if row:
            return {"id": row[0], "status": row[1]}
        return {}

    def get_combinations_for_export(work_id: str):
        # No combinations table in minimal test db
        return []

    def get_executions_for_export(work_id: str):
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT * FROM executions")
        rows = cur.fetchall()
        results = []
        for row in rows:
            results.append(
                {
                    "id": row[0],
                    "combination_id": row[1],
                    "status": row[2],
                    "params_json": row[3],
                    "result_json": row[4],
                    "objective": row[5],
                    "unit_id": row[6],
                    "sequencia": row[7],
                }
            )
        return results

    def get_execution_progress_for_export(work_id: str):
        return []

    def get_events_for_export(work_id: str):
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT * FROM events")
        rows = cur.fetchall()
        results = []
        for row in rows:
            results.append(
                {
                    "id": row[0],
                    "event_type": row[1],
                    "entity_data_json": row[2],
                    "timestamp": row[3],
                }
            )
        return results

    def get_datasets_for_export(work_id: str):
        return []

    def get_dataset_sequences_for_export(work_id: str):
        return []

    def get_optimization_executions_for_export(work_id: str):
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute("SELECT * FROM executions WHERE unit_id LIKE 'optimization:%'")
        rows = cur.fetchall()
        results = []
        for row in rows:
            results.append(
                {
                    "id": row[0],
                    "combination_id": row[1],
                    "status": row[2],
                    "params_json": row[3],
                    "result_json": row[4],
                    "objective": row[5],
                    "unit_id": row[6],
                    "sequencia": row[7],
                }
            )
        return results

    def get_sensitivity_events_for_export(work_id: str):
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute(
            "SELECT * FROM events WHERE event_type = 'progress' AND entity_data_json LIKE '%\"unit_id\": \"sensitivity_analysis\"%'"
        )
        rows = cur.fetchall()
        results = []
        for row in rows:
            results.append(
                {
                    "id": row[0],
                    "event_type": row[1],
                    "entity_data_json": row[2],
                    "timestamp": row[3],
                }
            )
        return results

    mock_store = Mock()
    mock_store._engine = mock_engine
    mock_store.get_work_export_data = get_work_export_data
    mock_store.get_combinations_for_export = get_combinations_for_export
    mock_store.get_executions_for_export = get_executions_for_export
    mock_store.get_execution_progress_for_export = get_execution_progress_for_export
    mock_store.get_events_for_export = get_events_for_export
    mock_store.get_datasets_for_export = get_datasets_for_export
    mock_store.get_dataset_sequences_for_export = get_dataset_sequences_for_export
    mock_store.get_optimization_executions_for_export = (
        get_optimization_executions_for_export
    )
    mock_store.get_sensitivity_events_for_export = get_sensitivity_events_for_export

    mock_work_store = Mock()
    mock_work_store.store = mock_store

    return mock_work_store


def create_minimal_db(tmp_path: Path):
    import sqlite3

    db_path = tmp_path / "work.db"
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    # minimal tables used
    cur.execute(
        "CREATE TABLE executions (id INTEGER PRIMARY KEY, combination_id INTEGER, status TEXT, params_json TEXT, result_json TEXT, objective REAL, unit_id TEXT, sequencia INTEGER)"
    )
    cur.execute(
        "CREATE TABLE events (id INTEGER PRIMARY KEY, event_type TEXT, entity_data_json TEXT, timestamp TEXT)"
    )
    cur.execute("CREATE TABLE work (id TEXT, status TEXT)")
    # insert work row
    cur.execute("INSERT INTO work (id, status) VALUES (?, ?)", ("w1", "completed"))
    # optimization executions
    for i in range(3):
        cur.execute(
            "INSERT INTO executions (combination_id, status, params_json, result_json, objective, unit_id, sequencia) VALUES (?,?,?,?,?,?,?)",
            (
                1,
                "completed",
                json.dumps({"lr": 0.1 * (i + 1)}),
                json.dumps({"metric": 1.0 / (i + 1)}),
                1.0 / (i + 1),
                f"optimization:taskA:dsX:preset1:AlgoY:{i}",
                i,
            ),
        )
    # sensitivity event
    sens_payload = {
        "unit_id": "sensitivity_analysis",
        "message": "Análise de sensibilidade concluída usando método morris",
        "context": {
            "method": "variance",
            "num_samples": 10,
            "param_scores": {"lr": 0.5, "batch_size": 0.2},
        },
    }
    cur.execute(
        "INSERT INTO events (event_type, entity_data_json, timestamp) VALUES (?,?,?)",
        ("progress", json.dumps(sens_payload), "2025-01-01T00:00:00Z"),
    )
    conn.commit()
    conn.close()
    return db_path


def test_finalization_generates_expected_artifacts(tmp_path):
    db_path = create_minimal_db(tmp_path)
    work_store = create_mock_work_store(db_path)
    output_dir = tmp_path / "out"
    cfg = FinalizationConfig(work_id="w1", output_dir=output_dir)
    svc = FinalizationService(cfg, work_store=work_store)
    svc.run()

    # core files
    assert (output_dir / "manifest.json").exists()
    assert (output_dir / "full_results.json").exists()
    assert (output_dir / "summary.md").exists()
    # raw db dumps include executions table
    assert (output_dir / "raw_db" / "executions.csv").exists()

    manifest = json.loads((output_dir / "manifest.json").read_text())
    assert manifest["work_id"] == "w1"
    # optimization section
    optuna_studies = manifest["artefacts"]["optimization"]
    assert len(optuna_studies) == 1
    assert optuna_studies[0]["trial_count"] == 3

    # sensitivity section
    sensitivity = manifest["artefacts"]["sensitivity"]
    assert len(sensitivity) >= 1

    # plots were generated
    study_dir = next(output_dir.glob("optuna/taskA__dsX__preset1__AlgoY"))
    assert (study_dir / "objective_vs_trial.png").exists()
    assert (study_dir / "best_objective_vs_trial.png").exists()

    sens_pngs = list((output_dir / "sensitivity" / "native").glob("*_scores.png"))
    assert sens_pngs, "Expected at least one sensitivity scores plot"

    # full results structure
    full_results = json.loads((output_dir / "full_results.json").read_text())
    assert full_results["optimization"]["studies"][0]["trial_count"] == 3
    assert full_results["sensitivity"]["analyses"], "Analyses list should not be empty"
