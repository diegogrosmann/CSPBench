import json
from pathlib import Path
from src.infrastructure.export.finalization_service import FinalizationConfig, FinalizationService


def create_db_no_optional(tmp_path: Path):
    import sqlite3
    db_path = tmp_path / "work.db"
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("CREATE TABLE executions (id INTEGER PRIMARY KEY, combination_id INTEGER, status TEXT, params_json TEXT, result_json TEXT, objective REAL, unit_id TEXT, sequencia INTEGER)")
    cur.execute("CREATE TABLE events (id INTEGER PRIMARY KEY, event_type TEXT, entity_data_json TEXT, timestamp TEXT)")
    cur.execute("CREATE TABLE work (id TEXT, status TEXT)")
    cur.execute("INSERT INTO work (id, status) VALUES (?, ?)", ("w2", "completed"))
    conn.commit()
    conn.close()
    return db_path


def test_finalization_without_optimization_or_sensitivity(tmp_path):
    db_path = create_db_no_optional(tmp_path)
    output_dir = tmp_path / "out"
    cfg = FinalizationConfig(work_id="w2", db_path=db_path, output_dir=output_dir, tool_version="test-version")
    FinalizationService(cfg).run()

    manifest = json.loads((output_dir / "manifest.json").read_text())
    assert manifest["artefacts"]["optimization"] == []
    assert manifest["artefacts"]["sensitivity"] == []

    full_results = json.loads((output_dir / "full_results.json").read_text())
    assert full_results["optimization"]["studies"] == []
    assert full_results["sensitivity"]["analyses"] == []

