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
import logging
from typing import Any, Iterable, Sequence, Optional

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
        work_id TEXT PRIMARY KEY,
        json TEXT NOT NULL,
        FOREIGN KEY(work_id) REFERENCES work(id) ON DELETE CASCADE
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
        sequencia INTEGER,
        status TEXT,
        started_at REAL,
        finished_at REAL,
        params_json TEXT,
        result_json TEXT,
        objective REAL
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS pipeline_combinations(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        work_id TEXT,
        task_index INTEGER,
        dataset_index INTEGER,
        preset_index INTEGER,
        algorithm_index INTEGER,
        status TEXT CHECK(status IN ('pending', 'running', 'completed', 'failed')),
        total_sequences INTEGER,
        created_at REAL,
        started_at REAL,
        finished_at REAL,
        FOREIGN KEY(work_id) REFERENCES work(id) ON DELETE CASCADE,
        UNIQUE(work_id, task_index, dataset_index, preset_index, algorithm_index)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS logs(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        work_id TEXT,
        level TEXT,
        message TEXT,
        context_json TEXT,
        timestamp REAL,
        component TEXT,
        FOREIGN KEY(work_id) REFERENCES work(id) ON DELETE CASCADE
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS progress_events(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        work_id TEXT,
        event_type TEXT,
        entity_id TEXT,
        entity_data_json TEXT,
        timestamp REAL,
        FOREIGN KEY(work_id) REFERENCES work(id) ON DELETE CASCADE
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS unit_events(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        work_id TEXT,
        unit_id TEXT,
        event_type TEXT CHECK(event_type IN ('started', 'finished', 'error', 'callback')),
        event_data_json TEXT,
        timestamp REAL,
        FOREIGN KEY(work_id) REFERENCES work(id) ON DELETE CASCADE
    )
    """,
]

class WorkStatePersistence:
    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._lock = threading.RLock()
        self._init_done = False
        self._logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

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

    def save_config(self, work_id: str, config_obj: Any) -> None:
        """Salva configuração associada ao work."""
        if is_dataclass(config_obj):
            data = asdict(config_obj)
        else:
            data = config_obj
        js = json.dumps(data, ensure_ascii=False)
        self._execute(
            "INSERT OR REPLACE INTO config(work_id, json) VALUES(?, ?)", (work_id, js)
        )

    def get_config(self, work_id: str) -> dict[str, Any] | None:
        """Recupera configuração do work."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute("SELECT json FROM config WHERE work_id=?", (work_id,))
            row = cursor.fetchone()
            if row:
                try:
                    return json.loads(row[0])
                except json.JSONDecodeError:
                    return None
            return None

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
        sequencia: int | None = None,
        params: dict[str, Any] | None = None,
    ) -> None:
        params_js = json.dumps(params or {}, ensure_ascii=False)
        self._execute(
            """
            INSERT OR REPLACE INTO executions(
              unit_id, task_id, dataset_id, algorithm, mode, sequencia,
              status, started_at, finished_at, params_json, result_json, objective
            ) VALUES(?,?,?,?,?,?,?, ?, NULL, ?, NULL, NULL)
            """,
            (
                unit_id,
                task_id,
                dataset_id,
                algorithm,
                mode,
                sequencia,
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

    def has_dataset(self, dataset_id: str) -> bool:
        """Verifica se dataset já está persistido com sequências."""
        if not dataset_id:
            return False
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                "SELECT COUNT(*) FROM dataset_sequences WHERE dataset_id=?",
                (dataset_id,),
            )
            count = cursor.fetchone()[0]
            return count > 0

    def get_completed_repetitions(
        self, task_id: str, dataset_id: str, algorithm: str
    ) -> list[dict[str, Any]]:
        """Obtém sequencias já completadas para funcionalidade de resume."""
        if not all([task_id, dataset_id, algorithm]):
            return []

        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                """
                SELECT sequencia, result_json, objective, params_json 
                FROM executions 
                WHERE task_id=? AND dataset_id=? AND algorithm=? 
                  AND mode='experiment' AND status='ok' AND sequencia IS NOT NULL
                ORDER BY sequencia
                """,
                (task_id, dataset_id, algorithm),
            )
            rows = cursor.fetchall()

            results = []
            for row in rows:
                sequencia, result_json, objective, params_json = row
                try:
                    result = json.loads(result_json) if result_json else {}
                    params = json.loads(params_json) if params_json else {}
                    results.append(
                        {
                            "sequencia": sequencia,
                            "result": result,
                            "objective": objective,
                            "params": params,
                        }
                    )
                except json.JSONDecodeError:
                    continue  # Skip malformed entries

            return results

    def init_pipeline_combinations(
        self, 
        work_id: str, 
        tasks_combinations: list[dict[str, Any]]
    ) -> None:
        """Inicializa todas as combinações do pipeline."""
        if not tasks_combinations:
            return

        combinations = []
        timestamp = time.time()
        
        for combo in tasks_combinations:
            combinations.append((
                work_id,
                combo['task_index'],
                combo['dataset_index'], 
                combo['preset_index'],
                combo['algorithm_index'],
                'pending',
                combo.get('total_sequences', 1),
                timestamp,
                None,  # started_at
                None   # finished_at
            ))

        self._executemany(
            """
            INSERT OR IGNORE INTO pipeline_combinations(
                work_id, task_index, dataset_index, preset_index, algorithm_index,
                status, total_sequences, created_at, started_at, finished_at
            ) VALUES(?,?,?,?,?,?,?,?,?,?)
            """,
            combinations
        )

    def update_combination_status(
        self,
        work_id: str,
        task_index: int,
        dataset_index: int,
        preset_index: int,
        algorithm_index: int,
        status: str
    ) -> None:
        """Atualiza status de uma combinação específica."""
        timestamp = time.time()
        
        if status == 'running':
            self._execute(
                """
                UPDATE pipeline_combinations 
                SET status=?, started_at=? 
                WHERE work_id=? AND task_index=? AND dataset_index=? 
                  AND preset_index=? AND algorithm_index=?
                """,
                (status, timestamp, work_id, task_index, dataset_index, preset_index, algorithm_index)
            )
        elif status in ('completed', 'failed'):
            self._execute(
                """
                UPDATE pipeline_combinations 
                SET status=?, finished_at=? 
                WHERE work_id=? AND task_index=? AND dataset_index=? 
                  AND preset_index=? AND algorithm_index=?
                """,
                (status, timestamp, work_id, task_index, dataset_index, preset_index, algorithm_index)
            )
        else:
            self._execute(
                """
                UPDATE pipeline_combinations 
                SET status=? 
                WHERE work_id=? AND task_index=? AND dataset_index=? 
                  AND preset_index=? AND algorithm_index=?
                """,
                (status, work_id, task_index, dataset_index, preset_index, algorithm_index)
            )

    def get_next_pending_combination(self, work_id: str) -> dict[str, Any] | None:
        """Obtém próxima combinação pendente para execução."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                """
                SELECT task_index, dataset_index, preset_index, algorithm_index, total_sequences
                FROM pipeline_combinations 
                WHERE work_id=? AND status='pending'
                ORDER BY task_index, dataset_index, preset_index, algorithm_index
                LIMIT 1
                """,
                (work_id,)
            )
            row = cursor.fetchone()
            if row:
                return {
                    "task_index": row[0],
                    "dataset_index": row[1], 
                    "preset_index": row[2],
                    "algorithm_index": row[3],
                    "total_sequences": row[4]
                }
            return None

    def get_pipeline_progress(self, work_id: str) -> dict[str, Any]:
        """Obtém progresso atual do pipeline."""
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                """
                SELECT status, COUNT(*) 
                FROM pipeline_combinations 
                WHERE work_id=? 
                GROUP BY status
                """,
                (work_id,)
            )
            
            progress = {'pending': 0, 'running': 0, 'completed': 0, 'failed': 0}
            for status, count in cursor.fetchall():
                progress[status] = count
            
            total = sum(progress.values())
            completion_rate = (progress['completed'] / total * 100) if total > 0 else 0
            
            return {
                **progress,
                'total': total,
                'completion_rate': completion_rate
            }

    def unit_started(self, work_id: str, unit_id: str, unit_data: dict[str, Any]) -> None:
        """Registra início de uma unidade de execução."""
        self._execute(
            "INSERT INTO unit_events(work_id, unit_id, event_type, event_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, unit_id, "started", json.dumps(unit_data, ensure_ascii=False), time.time())
        )

    def unit_finished(self, work_id: str, unit_id: str, result: dict[str, Any]) -> None:
        """Registra fim de uma unidade de execução."""
        self._execute(
            "INSERT INTO unit_events(work_id, unit_id, event_type, event_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, unit_id, "finished", json.dumps(result, ensure_ascii=False), time.time())
        )

    def error(self, work_id: str, unit_id: str, error: Exception) -> None:
        """Registra erro em uma unidade."""
        error_data = {
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": __import__('traceback').format_exc()
        }
        self._execute(
            "INSERT INTO unit_events(work_id, unit_id, event_type, event_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, unit_id, "error", json.dumps(error_data, ensure_ascii=False), time.time())
        )

    def unit_callback(self, work_id: str, unit_id: str, callback_data: dict[str, Any]) -> None:
        """Registra callback de uma unidade."""
        self._execute(
            "INSERT INTO unit_events(work_id, unit_id, event_type, event_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, unit_id, "callback", json.dumps(callback_data, ensure_ascii=False), time.time())
        )

    def get_unit_events(self, work_id: str, unit_id: str | None = None, event_type: str | None = None, limit: int = 100) -> list[dict[str, Any]]:
        """Obtém eventos de unidades."""
        conditions = ["work_id = ?"]
        params = [work_id]
        
        if unit_id:
            conditions.append("unit_id = ?")
            params.append(unit_id)
        
        if event_type:
            conditions.append("event_type = ?")
            params.append(event_type)
        
        where_clause = " AND ".join(conditions)
        params.append(limit)
        
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(f"""
                SELECT unit_id, event_type, event_data_json, timestamp
                FROM unit_events 
                WHERE {where_clause}
                ORDER BY timestamp DESC 
                LIMIT ?
            """, params)
            
            results = []
            for row in cursor.fetchall():
                try:
                    event_data = json.loads(row[2]) if row[2] else {}
                    results.append({
                        "unit_id": row[0],
                        "event_type": row[1],
                        "event_data": event_data,
                        "timestamp": row[3],
                    })
                except json.JSONDecodeError:
                    continue
            
            return results

    def get_execution_history(
        self, 
        task_id: Optional[str] = None,
        dataset_id: Optional[str] = None,
        algorithm: Optional[str] = None,
        limit: int = 100
    ) -> list[dict[str, Any]]:
        """Obtém histórico de execuções com filtros opcionais."""
        conditions = []
        params = []
        
        if task_id:
            conditions.append("task_id = ?")
            params.append(task_id)
        if dataset_id:
            conditions.append("dataset_id = ?")
            params.append(dataset_id)
        if algorithm:
            conditions.append("algorithm = ?")
            params.append(algorithm)
        
        where_clause = " AND ".join(conditions) if conditions else "1=1"
        params.append(limit)
        
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(f"""
                SELECT unit_id, task_id, dataset_id, algorithm, mode, sequencia,
                       status, started_at, finished_at, 
                       params_json, result_json, objective
                FROM executions 
                WHERE {where_clause}
                ORDER BY started_at DESC 
                LIMIT ?
            """, params)
            
            results = []
            for row in cursor.fetchall():
                try:
                    params_data = json.loads(row[9]) if row[9] else {}
                    result_data = json.loads(row[10]) if row[10] else {}
                    
                    results.append({
                        "unit_id": row[0],
                        "task_id": row[1],
                        "dataset_id": row[2],
                        "algorithm": row[3],
                        "mode": row[4],
                        "sequencia": row[5],
                        "status": row[6],
                        "started_at": row[7],
                        "finished_at": row[8],
                        "params": params_data,
                        "result": result_data,
                        "objective": row[11],
                        "duration": (row[8] - row[7]) if row[8] and row[7] else None,
                    })
                except (json.JSONDecodeError, TypeError) as e:
                    self._logger.warning(f"Failed to parse execution data for {row[0]}: {e}")
                    continue
            
            return results

    def get_work_statistics(self, work_id: str) -> dict[str, Any]:
        """Obtém estatísticas completas do work."""
        if not work_id:
            return {}
        
        with self._lock:
            cursor = self._conn.cursor()
            
            # Estatísticas básicas
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_executions,
                    COUNT(CASE WHEN status = 'ok' THEN 1 END) as completed,
                    COUNT(CASE WHEN status = 'error' THEN 1 END) as failed,
                    COUNT(CASE WHEN status = 'running' THEN 1 END) as running,
                    AVG(CASE WHEN finished_at IS NOT NULL AND started_at IS NOT NULL 
                             THEN finished_at - started_at END) as avg_duration,
                    COUNT(DISTINCT dataset_id) as datasets_count,
                    COUNT(DISTINCT algorithm) as algorithms_count
                FROM executions
                WHERE unit_id LIKE ?
            """, (f"{work_id}:%",))
            stats = cursor.fetchone()
            
            # Progresso por modo
            cursor.execute("""
                SELECT mode, status, COUNT(*) 
                FROM executions 
                WHERE unit_id LIKE ?
                GROUP BY mode, status
            """, (f"{work_id}:%",))
            mode_progress = {}
            for mode, status, count in cursor.fetchall():
                if mode not in mode_progress:
                    mode_progress[mode] = {}
                mode_progress[mode][status] = count
            
            # Progresso do pipeline
            pipeline_progress = self.get_pipeline_progress(work_id)
            
            return {
                "total_executions": stats[0] or 0,
                "completed": stats[1] or 0,
                "failed": stats[2] or 0,
                "running": stats[3] or 0,
                "avg_duration_seconds": stats[4],
                "datasets_count": stats[5] or 0,
                "algorithms_count": stats[6] or 0,
                "mode_progress": mode_progress,
                "completion_rate": (stats[1] or 0) / max(stats[0] or 1, 1) * 100,
                "pipeline_progress": pipeline_progress,
            }

    # Remove métodos obsoletos
    def save_pipeline_state(self, *args, **kwargs) -> None:
        """Método obsoleto - usar init_pipeline_combinations e update_combination_status."""
        pass

    def load_pipeline_state(self, work_id: str) -> dict[str, Any] | None:
        """Método obsoleto - usar get_next_pending_combination."""
        return self.get_next_pending_combination(work_id)

    def update_pipeline_status(self, work_id: str, status: str) -> None:
        """Método obsoleto - usar update_work_status."""
        self.update_work_status(work_id, status)

    def clear_pipeline_state(self, work_id: str) -> None:
        """Método obsoleto - combinações são mantidas para histórico."""
        pass

    def has_pending_executions(self, work_id: str) -> bool:
        """Verifica se há execuções pendentes/em andamento."""
        if not work_id:
            return False
        
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                "SELECT COUNT(*) FROM pipeline_combinations WHERE work_id=? AND status IN ('pending', 'running')",
                (work_id,)
            )
            count = cursor.fetchone()[0]
            return count > 0

    def cleanup_old_executions(self, days_old: int = 30) -> int:
        """Remove execuções antigas baseado na data de criação."""
        cutoff_time = time.time() - (days_old * 24 * 3600)
        
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                "DELETE FROM executions WHERE started_at < ? AND status IN ('ok', 'error')",
                (cutoff_time,)
            )
            deleted_count = cursor.rowcount
            self._conn.commit()
            
            self._logger.info(f"Cleaned up {deleted_count} old executions older than {days_old} days")
            return deleted_count

    def backup_database(self, backup_path: Optional[Path] = None) -> Path:
        """Cria backup do banco de dados."""
        if backup_path is None:
            timestamp = int(time.time())
            backup_path = self.db_path.parent / f"backup_{timestamp}.db"
        
        backup_path = Path(backup_path)
        backup_path.parent.mkdir(parents=True, exist_ok=True)
        
        with self._lock:
            try:
                backup_conn = sqlite3.connect(backup_path)
                self._conn.backup(backup_conn)
                backup_conn.close()
                self._logger.info(f"Database backed up to {backup_path}")
                return backup_path
            except Exception as e:
                self._logger.error(f"Failed to backup database: {e}")
                raise

    def get_execution_history(
        self, 
        task_id: Optional[str] = None,
        dataset_id: Optional[str] = None,
        algorithm: Optional[str] = None,
        limit: int = 100
    ) -> list[dict[str, Any]]:
        """Obtém histórico de execuções com filtros opcionais."""
        conditions = []
        params = []
        
        if task_id:
            conditions.append("task_id = ?")
            params.append(task_id)
        if dataset_id:
            conditions.append("dataset_id = ?")
            params.append(dataset_id)
        if algorithm:
            conditions.append("algorithm = ?")
            params.append(algorithm)
        
        where_clause = " AND ".join(conditions) if conditions else "1=1"
        params.append(limit)
        
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(f"""
                SELECT unit_id, task_id, dataset_id, algorithm, mode, sequencia,
                       status, started_at, finished_at, 
                       params_json, result_json, objective
                FROM executions 
                WHERE {where_clause}
                ORDER BY started_at DESC 
                LIMIT ?
            """, params)
            
            results = []
            for row in cursor.fetchall():
                try:
                    params_data = json.loads(row[9]) if row[9] else {}
                    result_data = json.loads(row[10]) if row[10] else {}
                    
                    results.append({
                        "unit_id": row[0],
                        "task_id": row[1],
                        "dataset_id": row[2],
                        "algorithm": row[3],
                        "mode": row[4],
                        "sequencia": row[5],
                        "status": row[6],
                        "started_at": row[7],
                        "finished_at": row[8],
                        "params": params_data,
                        "result": result_data,
                        "objective": row[11],
                        "duration": (row[8] - row[7]) if row[8] and row[7] else None,
                    })
                except (json.JSONDecodeError, TypeError) as e:
                    self._logger.warning(f"Failed to parse execution data for {row[0]}: {e}")
                    continue
            
            return results

    def vacuum_database(self) -> None:
        """Otimiza o banco de dados removendo espaço não utilizado."""
        with self._lock:
            try:
                self._conn.execute("VACUUM")
                self._conn.commit()
                self._logger.info("Database vacuum completed")
            except Exception as e:
                self._logger.error(f"Failed to vacuum database: {e}")
                raise

    def get_database_info(self) -> dict[str, Any]:
        """Obtém informações sobre o banco de dados."""
        with self._lock:
            cursor = self._conn.cursor()
            
            # Tamanho das tabelas
            cursor.execute("""
                SELECT name, COUNT(*) as row_count
                FROM sqlite_master 
                WHERE type='table' AND name NOT LIKE 'sqlite_%'
            """)
            tables = {}
            for table_name, _ in cursor.fetchall():
                cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
                count = cursor.fetchone()[0]
                tables[table_name] = count
            
            # Tamanho do arquivo
            file_size = self.db_path.stat().st_size if self.db_path.exists() else 0
            
            return {
                "database_path": str(self.db_path),
                "file_size_bytes": file_size,
                "file_size_mb": round(file_size / 1024 / 1024, 2),
                "table_counts": tables,
                "wal_mode": True,  # Sempre WAL baseado na inicialização
            }

    def log(self, work_id: str, level: str, message: str, context: dict[str, Any] | None = None, component: str = "CSPBench") -> None:
        """Registra uma mensagem de log."""
        context_json = json.dumps(context or {}, ensure_ascii=False)
        self._execute(
            "INSERT INTO logs(work_id, level, message, context_json, timestamp, component) VALUES(?,?,?,?,?,?)",
            (work_id, level, message, context_json, time.time(), component)
        )

    def pipeline_started(self, work_id: str, config: Any) -> None:
        """Registra início do pipeline."""
        config_data = asdict(config) if is_dataclass(config) else config
        self._execute(
            "INSERT INTO progress_events(work_id, event_type, entity_id, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, "pipeline_started", "pipeline", json.dumps(config_data, ensure_ascii=False), time.time())
        )

    def pipeline_finished(self, work_id: str, success: bool) -> None:
        """Registra fim do pipeline."""
        self._execute(
            "INSERT INTO progress_events(work_id, event_type, entity_id, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, "pipeline_finished", "pipeline", json.dumps({"success": success}, ensure_ascii=False), time.time())
        )

    def task_started(self, work_id: str, task_id: str, task_data: dict[str, Any]) -> None:
        """Registra início de uma tarefa."""
        self._execute(
            "INSERT INTO progress_events(work_id, event_type, entity_id, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, "task_started", task_id, json.dumps(task_data, ensure_ascii=False), time.time())
        )

    def task_finished(self, work_id: str, task_id: str, task_data: dict[str, Any]) -> None:
        """Registra fim de uma tarefa."""
        self._execute(
            "INSERT INTO progress_events(work_id, event_type, entity_id, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, "task_finished", task_id, json.dumps(task_data, ensure_ascii=False), time.time())
        )

    def unit_started(self, work_id: str, unit_id: str, unit_data: dict[str, Any]) -> None:
        """Registra início de uma unidade de execução."""
        self._execute(
            "INSERT INTO progress_events(work_id, event_type, entity_id, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, "unit_started", unit_id, json.dumps(unit_data, ensure_ascii=False), time.time())
        )

    def unit_finished(self, work_id: str, unit_id: str, result: dict[str, Any]) -> None:
        """Registra fim de uma unidade de execução."""
        self._execute(
            "INSERT INTO progress_events(work_id, event_type, entity_id, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, "unit_finished", unit_id, json.dumps(result, ensure_ascii=False), time.time())
        )

    def error(self, work_id: str, unit_id: str, error: Exception) -> None:
        """Registra erro em uma unidade."""
        error_data = {
            "error_type": error.__class__.__name__,
            "error_message": str(error),
            "traceback": __import__('traceback').format_exc()
        }
        self._execute(
            "INSERT INTO progress_events(work_id, event_type, entity_id, entity_data_json, timestamp) VALUES(?,?,?,?,?)",
            (work_id, "unit_error", unit_id, json.dumps(error_data, ensure_ascii=False), time.time())
        )

    def get_logs(self, work_id: str, level: str | None = None, limit: int = 100) -> list[dict[str, Any]]:
        """Obtém logs do work."""
        conditions = ["work_id = ?"]
        params = [work_id]
        
        if level:
            conditions.append("level = ?")
            params.append(level)
        
        where_clause = " AND ".join(conditions)
        params.append(limit)
        
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(f"""
                SELECT level, message, context_json, timestamp, component
                FROM logs 
                WHERE {where_clause}
                ORDER BY timestamp DESC 
                LIMIT ?
            """, params)
            
            results = []
            for row in cursor.fetchall():
                try:
                    context = json.loads(row[2]) if row[2] else {}
                    results.append({
                        "level": row[0],
                        "message": row[1],
                        "context": context,
                        "timestamp": row[3],
                        "component": row[4],
                    })
                except json.JSONDecodeError:
                    continue
            
            return results

    def get_progress_events(self, work_id: str, event_type: str | None = None, limit: int = 100) -> list[dict[str, Any]]:
        """Obtém eventos de progresso."""
        conditions = ["work_id = ?"]
        params = [work_id]
        
        if event_type:
            conditions.append("event_type = ?")
            params.append(event_type)
        
        where_clause = " AND ".join(conditions)
        params.append(limit)
        
        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(f"""
                SELECT event_type, entity_id, entity_data_json, timestamp
                FROM progress_events 
                WHERE {where_clause}
                ORDER BY timestamp DESC 
                LIMIT ?
            """, params)
            
            results = []
            for row in cursor.fetchall():
                try:
                    entity_data = json.loads(row[2]) if row[2] else {}
                    results.append({
                        "event_type": row[0],
                        "entity_id": row[1],
                        "entity_data": entity_data,
                        "timestamp": row[3],
                    })
                except json.JSONDecodeError:
                    continue
            
            return results
