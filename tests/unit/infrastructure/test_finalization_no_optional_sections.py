import json
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

    mock_store = Mock()
    mock_store._engine = mock_engine

    # Mock the export query methods to return empty results
    mock_store.get_work_export_data.return_value = []
    mock_store.get_combinations_for_export.return_value = []
    mock_store.get_executions_for_export.return_value = []
    mock_store.get_execution_progress_for_export.return_value = []
    mock_store.get_events_for_export.return_value = []
    mock_store.get_datasets_for_export.return_value = []
    mock_store.get_dataset_sequences_for_export.return_value = []
    mock_store.get_optimization_executions_for_export.return_value = []
    mock_store.get_sensitivity_events_for_export.return_value = []

    mock_work_store = Mock()
    mock_work_store.store = mock_store

    return mock_work_store


def create_db_no_optional(tmp_path: Path):
    import sqlite3

    db_path = tmp_path / "work.db"
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute(
        "CREATE TABLE executions (id INTEGER PRIMARY KEY, combination_id INTEGER, status TEXT, params_json TEXT, result_json TEXT, objective REAL, unit_id TEXT, sequencia INTEGER)"
    )
    cur.execute(
        "CREATE TABLE events (id INTEGER PRIMARY KEY, event_type TEXT, entity_data_json TEXT, timestamp TEXT)"
    )
    cur.execute("CREATE TABLE work (id TEXT, status TEXT)")
    cur.execute("INSERT INTO work (id, status) VALUES (?, ?)", ("w2", "completed"))
    conn.commit()
    conn.close()
    return db_path


def test_finalization_without_optimization_or_sensitivity(tmp_path):
    db_path = create_db_no_optional(tmp_path)
    work_store = create_mock_work_store(db_path)
    output_dir = tmp_path / "out"
    cfg = FinalizationConfig(
        work_id="w2",
        output_dir=output_dir,
        tool_version="test-version",
    )
    FinalizationService(cfg, work_store=work_store).run()

    manifest = json.loads((output_dir / "manifest.json").read_text())
    assert manifest["artefacts"]["optimization"] == []
    assert manifest["artefacts"]["sensitivity"] == []

    full_results = json.loads((output_dir / "full_results.json").read_text())
    assert full_results["optimization"]["studies"] == []
    assert full_results["sensitivity"]["analyses"] == []
