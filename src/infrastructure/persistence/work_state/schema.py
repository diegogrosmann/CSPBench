"""Database schema definitions for work state persistence."""

SCHEMA_STATEMENTS = [
    """
    CREATE TABLE IF NOT EXISTS work(
        id TEXT PRIMARY KEY,
        config_json TEXT,
        status TEXT CHECK(status IN ('queued', 'running', 'paused', 'canceled', 'completed', 'failed', 'error')),
        created_at REAL,
        updated_at REAL,
        output_path TEXT,
        error TEXT,
        extra_json TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS datasets(
        id TEXT PRIMARY KEY,
        name TEXT,
        meta_json TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS dataset_sequences(
        dataset_id TEXT,
        seq_index INTEGER,
        sequence TEXT,
        PRIMARY KEY(dataset_id, seq_index),
        FOREIGN KEY(dataset_id) REFERENCES datasets(id) ON DELETE CASCADE
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS combinations(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        work_id TEXT,
        task_id TEXT,
        dataset_id TEXT,
        preset_id TEXT,
        algorithm_id TEXT,
        mode TEXT,
        status TEXT CHECK(status IN ('queued', 'running', 'paused', 'canceled', 'completed', 'failed', 'error')),
        total_sequences INTEGER,
        created_at REAL,
        started_at REAL,
        finished_at REAL,
        FOREIGN KEY(work_id) REFERENCES work(id) ON DELETE CASCADE,
        FOREIGN KEY(dataset_id) REFERENCES datasets(id) ON DELETE SET NULL,
        UNIQUE(work_id, task_id, dataset_id, preset_id, algorithm_id, mode)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS executions(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        unit_id TEXT UNIQUE NOT NULL,
        combination_id INTEGER,
        sequencia INTEGER,
        status TEXT CHECK(status IN ('queued', 'running', 'paused', 'canceled', 'completed', 'failed', 'error')),
        started_at REAL,
        finished_at REAL,
        params_json TEXT,
        result_json TEXT,
        objective REAL,
        FOREIGN KEY(combination_id) REFERENCES combinations(id) ON DELETE CASCADE
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS execution_progress(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        execution_id INTEGER,
        progress REAL CHECK(progress >= 0.0 AND progress <= 1.0),
        message TEXT,
        timestamp REAL,
        FOREIGN KEY(execution_id) REFERENCES executions(id) ON DELETE CASCADE
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS events(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        work_id TEXT,
        event_type TEXT CHECK(event_type IN ('error', 'progress', 'warning')),
        event_category TEXT CHECK(event_category IN ('work', 'task', 'dataset', 'preset', 'combination', 'unit', 'other')),
        entity_data_json TEXT,
        timestamp REAL,
        FOREIGN KEY(work_id) REFERENCES work(id) ON DELETE CASCADE
    )
    """,
    """
    CREATE INDEX IF NOT EXISTS idx_events_work_category ON events(work_id, event_category)
    """,
    """
    CREATE INDEX IF NOT EXISTS idx_events_timestamp ON events(timestamp)
    """,
    """
    CREATE INDEX IF NOT EXISTS idx_execution_progress_execution_id ON execution_progress(execution_id)
    """,
]
