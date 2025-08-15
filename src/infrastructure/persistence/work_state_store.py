"""WorkStateStore: persistência por work em arquivo SQLite outputs/<work_id>/state.db.

Responsável por:
 - Schema inicial (tabelas: work, config, datasets, dataset_sequences, executions, progress)
 - Armazenar config (JSON)
 - Registrar datasets e sequências (quando disponíveis)
 - Registrar início/fim de cada unidade de execução (experiment / optimization / sensitivity)
 - Atualizar status do work

Este é um MVP para possibilitar futura retomada do pipeline.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
import json
import sqlite3
import threading
import time
from typing import Any, Iterable, Sequence

_SCHEMA = [
    """
    CREATE TABLE IF NOT EXISTS work(
        id TEXT PRIMARY KEY,
        status TEXT,
        created_at REAL,
        updated_at REAL,
        output_path TEXT,
        error TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS config(
        id TEXT PRIMARY KEY CHECK (id='config'),
        json TEXT NOT NULL
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS datasets(
        id TEXT PRIMARY KEY,
        name TEXT,
        type TEXT,
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
    CREATE TABLE IF NOT EXISTS executions(
        unit_id TEXT PRIMARY KEY,
        task_id TEXT,
        dataset_id TEXT,
        algorithm TEXT,
        mode TEXT,
        repetition INTEGER,
        trial INTEGER,
        sample INTEGER,
        status TEXT,
        started_at REAL,
        finished_at REAL,
        params_json TEXT,
        result_json TEXT,
        objective REAL
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS progress(
        key TEXT PRIMARY KEY,
        value REAL,
        updated_at REAL
    )
    """,
]


class WorkStateStore:
    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._lock = threading.RLock()
        self._init_done = False

    # -- internal helpers --
    def _execute(self, sql: str, params: Sequence[Any] | None = None) -> None:
        with self._lock:
            self._conn.execute(sql, params or [])
            self._conn.commit()

    def _executemany(self, sql: str, seq_of_params: Iterable[Sequence[Any]]):
        with self._lock:
            self._conn.executemany(sql, seq_of_params)
            self._conn.commit()

    # -- public API --
    def init_schema(self) -> None:
        if self._init_done:
            return
        with self._lock:
            cur = self._conn.cursor()
            for stmt in _SCHEMA:
                cur.execute(stmt)
            self._conn.commit()
            self._init_done = True

    def save_config(self, config_obj: Any) -> None:
        if is_dataclass(config_obj):
            data = asdict(config_obj)
        else:
            # Pode ser dict já
            data = config_obj
        js = json.dumps(data, ensure_ascii=False)
        self._execute(
            "INSERT OR REPLACE INTO config(id, json) VALUES('config', ?)", (js,)
        )

    def upsert_work(self, work_dict: dict[str, Any]) -> None:
        # work_dict vem de WorkManager.get()
        self._execute(
            """
            INSERT INTO work(id, status, created_at, updated_at, output_path, error)
            VALUES(?,?,?,?,?,?)
            ON CONFLICT(id) DO UPDATE SET
              status=excluded.status,
              updated_at=excluded.updated_at,
              output_path=excluded.output_path,
              error=excluded.error
            """,
            (
                work_dict.get("id"),
                work_dict.get("status"),
                work_dict.get("created_at"),
                work_dict.get("updated_at"),
                work_dict.get("output_path"),
                work_dict.get("error"),
            ),
        )

    def update_work_status(self, work_id: str, status: str, **fields: Any) -> None:
        cols = ["status=?", "updated_at=?"]
        params: list[Any] = [status, time.time()]
        for k, v in fields.items():
            cols.append(f"{k}=?")
            params.append(v)
        params.append(work_id)
        sql = f"UPDATE work SET {', '.join(cols)} WHERE id=?"
        self._execute(sql, params)

    def ensure_dataset(self, dataset_obj: Any) -> None:
        ds_id = getattr(dataset_obj, "id", None)
        if not ds_id:
            return
        name = getattr(dataset_obj, "name", ds_id)
        dtype = getattr(dataset_obj, "type", dataset_obj.__class__.__name__)
        meta = {}
        # Tentar capturar campos simples
        for attr in [
            "n",
            "L",
            "alphabet",
            "noise",
            "filename",
            "query",
            "db",
            "retmax",
        ]:
            if hasattr(dataset_obj, attr):
                meta[attr] = getattr(dataset_obj, attr)
        js = json.dumps(meta, ensure_ascii=False)
        self._execute(
            "INSERT OR IGNORE INTO datasets(id, name, type, meta_json) VALUES(?,?,?,?)",
            (ds_id, name, dtype, js),
        )

    def add_sequences(self, dataset_id: str, sequences: list[str]):
        if not sequences:
            return
        rows = [(dataset_id, idx, seq) for idx, seq in enumerate(sequences)]
        self._executemany(
            "INSERT OR IGNORE INTO dataset_sequences(dataset_id, seq_index, sequence) VALUES(?,?,?)",
            rows,
        )

    def record_execution_start(
        self,
        *,
        unit_id: str,
        task_id: str,
        dataset_id: str,
        algorithm: str,
        mode: str,
        repetition: int | None = None,
        trial: int | None = None,
        sample: int | None = None,
        params: dict[str, Any] | None = None,
    ) -> None:
        params_js = json.dumps(params or {}, ensure_ascii=False)
        self._execute(
            """
            INSERT OR REPLACE INTO executions(
              unit_id, task_id, dataset_id, algorithm, mode, repetition, trial, sample,
              status, started_at, finished_at, params_json, result_json, objective
            ) VALUES(?,?,?,?,?,?,?,?,?, ?, NULL, ?, NULL, NULL)
            """,
            (
                unit_id,
                task_id,
                dataset_id,
                algorithm,
                mode,
                repetition,
                trial,
                sample,
                "running",
                time.time(),
                params_js,
            ),
        )

    def record_execution_end(
        self,
        unit_id: str,
        status: str,
        result: dict[str, Any] | None,
        objective: float | None,
    ) -> None:
        res_js = json.dumps(result or {}, ensure_ascii=False)
        self._execute(
            """
            UPDATE executions SET status=?, finished_at=?, result_json=?, objective=?
            WHERE unit_id=?
            """,
            (status, time.time(), res_js, objective, unit_id),
        )

    def close(self):
        with self._lock:
            try:
                self._conn.close()
            except Exception:  # noqa: BLE001
                pass


# --- Registro global simples (work_id -> store) ---
_STORES: dict[str, WorkStateStore] = {}
_STORES_LOCK = threading.RLock()


def register_work_store(work_id: str, store: WorkStateStore) -> None:
    with _STORES_LOCK:
        _STORES[work_id] = store


def get_work_store(work_id: str) -> WorkStateStore | None:
    with _STORES_LOCK:
        return _STORES.get(work_id)
