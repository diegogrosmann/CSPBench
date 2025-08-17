"""DEPRECATED: SQLiteWorkRepository removido - use WorkStateStore.

Este arquivo está obsoleto. A funcionalidade foi consolidada no WorkStateStore.
"""

from __future__ import annotations

# Este módulo foi removido - use WorkStateStore para persistência


class SQLiteWorkRepository(WorkRepository):
    """Persistência mínima em SQLite apenas para WorkItem (fase 1).

    Extensível depois para passos, execuções etc. Mantém thread-safety via lock + check_same_thread=False.
    """

    def __init__(self, db_path: str | Path = "cspbench.db"):
        self._db_path = str(db_path)
        Path(self._db_path).parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self._db_path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        self._lock = threading.RLock()
        self._init_schema()

    def _init_schema(self) -> None:
        with self._conn:  # autotransaction
            cur = self._conn.cursor()
            for stmt in _SCHEMA:
                cur.execute(stmt)

    # --- helpers ---
    def _serialize_config(self, cfg: CSPBenchConfig) -> str:
        # CSPBenchConfig provavelmente dataclass / pydantic-like com .model_dump ou .dict
        if hasattr(cfg, "model_dump"):
            data = cfg.model_dump()
        elif hasattr(cfg, "to_dict"):
            data = cfg.to_dict()
        elif hasattr(cfg, "__dict__"):
            data = {k: v for k, v in cfg.__dict__.items() if not k.startswith("_")}
        else:
            raise TypeError("Unsupported config object for serialization")
        return json.dumps(data, ensure_ascii=False)

    def _deserialize_config(self, raw: str) -> CSPBenchConfig:
        from src.domain.config import CSPBenchConfig  # import tardio para evitar ciclos

        data = json.loads(raw)
        return CSPBenchConfig(**data)  # type: ignore[arg-type]

    def add(self, item: WorkItem) -> None:
        with self._lock, self._conn:
            self._conn.execute(
                "INSERT OR REPLACE INTO work_items (id,status,created_at,updated_at,output_path,error,directory,config_json,extra_json) VALUES (?,?,?,?,?,?,?,?,?)",
                (
                    item.id,
                    item.status.value,
                    item.created_at,
                    item.updated_at,
                    item.output_path,
                    item.error,
                    item.directory,
                    self._serialize_config(item.config),
                    json.dumps(item.extra or {}, ensure_ascii=False),
                ),
            )

    def get(self, work_id: str) -> Optional[WorkItem]:
        with self._lock:
            cur = self._conn.execute("SELECT * FROM work_items WHERE id=?", (work_id,))
            row = cur.fetchone()
            if not row:
                return None
            return self._row_to_item(row)

    def list(self):
        with self._lock:
            cur = self._conn.execute(
                "SELECT * FROM work_items ORDER BY created_at DESC"
            )
            return [self._row_to_item(r) for r in cur.fetchall()]

    def update(self, item: WorkItem) -> None:
        with self._lock, self._conn:
            self._conn.execute(
                "UPDATE work_items SET status=?, updated_at=?, output_path=?, error=?, directory=?, config_json=?, extra_json=? WHERE id=?",
                (
                    item.status.value,
                    item.updated_at,
                    item.output_path,
                    item.error,
                    item.directory,
                    self._serialize_config(item.config),
                    json.dumps(item.extra or {}, ensure_ascii=False),
                    item.id,
                ),
            )

    def remove(self, work_id: str) -> None:
        with self._lock, self._conn:
            self._conn.execute("DELETE FROM work_items WHERE id=?", (work_id,))

    # --- internal ---
    def _row_to_item(self, row: sqlite3.Row) -> WorkItem:
        cfg = self._deserialize_config(row["config_json"])
        extra = json.loads(row["extra_json"]) if row["extra_json"] else {}
        return WorkItem(
            id=row["id"],
            config=cfg,
            status=WorkStatus(row["status"]),
            created_at=row["created_at"],
            updated_at=row["updated_at"],
            output_path=row["output_path"],
            error=row["error"],
            directory=row["directory"],
            extra=extra,
        )

    def close(self) -> None:
        with self._lock:
            self._conn.close()
