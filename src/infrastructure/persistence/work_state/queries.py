"""
Generic query service for work state database.
Provides reusable, type-safe queries for monitoring and reporting.
"""

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, Tuple


@dataclass
class ExecutionDetail:
    """Detailed execution information for monitoring."""
    unit_id: str
    combination_id: int
    sequencia: int
    status: str
    progress: float
    progress_message: Optional[str]
    started_at: Optional[float]
    finished_at: Optional[float]
    objective: Optional[float]
    task_id: str
    dataset_id: str
    preset_id: str
    algorithm_id: str
    mode: str
    total_sequences: int


@dataclass
class ProgressSummary:
    """Summary of progress for a work item."""
    work_id: str
    tasks: Dict[str, Union[int, str]]
    datasets: Dict[str, Union[int, str]]
    configs: Dict[str, Union[int, str]]
    algorithms: Dict[str, Union[int, str]]
    execution: Dict[str, int]
    global_execution: Dict[str, int]
    global_progress: float
    current_combination_details: Optional[Dict[str, Any]] = None


@dataclass
class ErrorSummary:
    """Summary of errors for a work item."""
    unit_id: str
    error_type: str
    error_message: str
    timestamp: float


class WorkStateQueries:
    """Generic query service for work state database."""
    
    def __init__(self, db_path: Union[str, Path]):
        """Initialize with database path."""
        self.db_path = Path(db_path)
        self._conn = None
    
    def _get_connection(self) -> sqlite3.Connection:
        """Get database connection with row factory."""
        if self._conn is None:
            self._conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
            self._conn.row_factory = sqlite3.Row
        return self._conn
    
    def close(self):
        """Close database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    # === PROGRESS QUERIES ===

    def get_running_combination(self, work_id: str) -> Optional[Dict[str, Any]]:
        """Retorna dict da primeira combinação em execução ou None."""
        cursor = self._get_connection().cursor()
        cursor.execute("""
            SELECT c.id as combination_id,
                   c.task_id,
                   c.dataset_id,
                   c.preset_id,
                   c.algorithm_id,
                   c.total_sequences
            FROM combinations c
            WHERE c.work_id = ? AND c.status = 'running'
            ORDER BY c.task_id, c.dataset_id, c.preset_id, c.algorithm_id
            LIMIT 1
        """, (work_id,))
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_task_status_lists(self, work_id: str) -> Tuple[List[str], List[str], List[str]]:
        """
        Retorna (finished_ids, running_ids, queued_ids) para tasks.
        finished: todas as combinações da task = completed
        queued: todas as combinações da task = queued
        running: qualquer outro caso (mistura / em execução / falhas / etc.)
        """
        cursor = self._get_connection().cursor()
        cursor.execute("""
            SELECT task_id,
                   -- considerar completed + failed + error como finalizados
                   SUM(CASE WHEN status IN ('completed','failed','error') THEN 1 ELSE 0 END) as completed_cnt,
                   SUM(CASE WHEN status='queued' THEN 1 ELSE 0 END) as queued_cnt,
                   COUNT(*) as total_cnt
            FROM combinations
            WHERE work_id = ?
            GROUP BY task_id
        """, (work_id,))
        finished, queued, running = [], [], []
        for row in cursor.fetchall():
            if row['completed_cnt'] == row['total_cnt']:
                finished.append(row['task_id'])
            elif row['queued_cnt'] == row['total_cnt']:
                queued.append(row['task_id'])
            else:
                running.append(row['task_id'])
        return finished, running, queued

    def get_dataset_status_lists(self, work_id: str, task_id: str) -> Tuple[List[str], List[str], List[str]]:
        """Mesmo critério de classificação aplicado a datasets dentro da task."""
        cursor = self._get_connection().cursor()
        cursor.execute("""
            SELECT dataset_id,
                   SUM(CASE WHEN status IN ('completed','failed','error') THEN 1 ELSE 0 END) as completed_cnt,
                   SUM(CASE WHEN status='queued' THEN 1 ELSE 0 END) as queued_cnt,
                   COUNT(*) as total_cnt
            FROM combinations
            WHERE work_id = ? AND task_id = ?
            GROUP BY dataset_id
        """, (work_id, task_id))
        finished, queued, running = [], [], []
        for row in cursor.fetchall():
            if row['completed_cnt'] == row['total_cnt']:
                finished.append(row['dataset_id'])
            elif row['queued_cnt'] == row['total_cnt']:
                queued.append(row['dataset_id'])
            else:
                running.append(row['dataset_id'])
        return finished, running, queued

    def get_config_status_lists(self, work_id: str, task_id: str, dataset_id: str) -> Tuple[List[str], List[str], List[str]]:
        """Classificação de presets (configs) dentro do dataset."""
        cursor = self._get_connection().cursor()
        cursor.execute("""
            SELECT preset_id,
                   SUM(CASE WHEN status IN ('completed','failed','error') THEN 1 ELSE 0 END) as completed_cnt,
                   SUM(CASE WHEN status='queued' THEN 1 ELSE 0 END) as queued_cnt,
                   COUNT(*) as total_cnt
            FROM combinations
            WHERE work_id = ? AND task_id = ? AND dataset_id = ?
            GROUP BY preset_id
        """, (work_id, task_id, dataset_id))
        finished, queued, running = [], [], []
        for row in cursor.fetchall():
            if row['completed_cnt'] == row['total_cnt']:
                finished.append(row['preset_id'])
            elif row['queued_cnt'] == row['total_cnt']:
                queued.append(row['preset_id'])
            else:
                running.append(row['preset_id'])
        return finished, running, queued

    def get_algorithm_status_lists(self, work_id: str, task_id: str, dataset_id: str, preset_id: str) -> Tuple[List[str], List[str], List[str]]:
        """Classificação de algorithms dentro da config."""
        cursor = self._get_connection().cursor()
        cursor.execute("""
            SELECT algorithm_id,
                   SUM(CASE WHEN status IN ('completed','failed','error') THEN 1 ELSE 0 END) as completed_cnt,
                   SUM(CASE WHEN status='queued' THEN 1 ELSE 0 END) as queued_cnt,
                   COUNT(*) as total_cnt
            FROM combinations
            WHERE work_id = ? AND task_id = ? AND dataset_id = ? AND preset_id = ?
            GROUP BY algorithm_id
        """, (work_id, task_id, dataset_id, preset_id))
        finished, queued, running = [], [], []
        for row in cursor.fetchall():
            if row['completed_cnt'] == row['total_cnt']:
                finished.append(row['algorithm_id'])
            elif row['queued_cnt'] == row['total_cnt']:
                queued.append(row['algorithm_id'])
            else:
                running.append(row['algorithm_id'])
        return finished, running, queued

    def get_execution_status_lists(self, combination_id: int) -> Tuple[List[int], List[int], List[int]]:
        """Retorna listas de execution IDs por status (finished, running, queued) para uma combinação."""
        cursor = self._get_connection().cursor()
        cursor.execute("""
            SELECT id, status FROM executions
            WHERE combination_id = ?
        """, (combination_id,))
        finished, running, queued = [], [], []
        for row in cursor.fetchall():
            if row['status'] in ('completed','failed','error'):
                finished.append(row['id'])
            elif row['status'] == 'running':
                running.append(row['id'])
            elif row['status'] == 'queued':
                queued.append(row['id'])
        return finished, running, queued

    def get_global_execution_stats(self, work_id: str) -> Dict[str, Union[int, float]]:
        """Estatísticas globais: execuções concluídas vs total (soma de total_sequences)."""
        cursor = self._get_connection().cursor()
        cursor.execute("""
            SELECT COALESCE(SUM(total_sequences),0) as total_sequences
            FROM combinations
            WHERE work_id = ?
        """, (work_id,))
        total_sequences = cursor.fetchone()["total_sequences"] or 0

        cursor.execute("""
            SELECT COUNT(*) as finished
            FROM executions e
            JOIN combinations c ON e.combination_id = c.id
            WHERE c.work_id = ? AND e.status IN ('completed','failed','error')
        """, (work_id,))
        finished = cursor.fetchone()["finished"] or 0

        progress = finished / total_sequences if total_sequences else 0.0
        return {"Finished": finished, "Total": total_sequences, "Progress": progress}

    def get_work_progress_summary(self, work_id: str) -> Optional[ProgressSummary]:
        """Resumo hierárquico usando métodos reutilizáveis."""
        if not self.work_exists(work_id):
            return None

        current_combo = self.get_running_combination(work_id)
        # TASKS
        task_finished, task_running, task_queued = self.get_task_status_lists(work_id)
        running_task_id = current_combo['task_id'] if current_combo else ""
        # DATASETS (somente da task corrente)
        if current_combo:
            dataset_finished, dataset_running, dataset_queued = self.get_dataset_status_lists(
                work_id, current_combo['task_id']
            )
            running_dataset_id = current_combo['dataset_id']
        else:
            dataset_finished = dataset_running = dataset_queued = []
            running_dataset_id = ""
        # CONFIGS
        if current_combo and current_combo['dataset_id']:
            config_finished, config_running, config_queued = self.get_config_status_lists(
                work_id, current_combo['task_id'], current_combo['dataset_id']
            )
            running_config_id = current_combo['preset_id']
        else:
            config_finished = config_running = config_queued = []
            running_config_id = ""
        # ALGORITHMS
        if current_combo and current_combo['preset_id']:
            alg_finished, alg_running, alg_queued = self.get_algorithm_status_lists(
                work_id, current_combo['task_id'], current_combo['dataset_id'], current_combo['preset_id']
            )
            running_alg_id = current_combo['algorithm_id']
        else:
            alg_finished = alg_running = alg_queued = []
            running_alg_id = ""
        # EXECUTIONS
        if current_combo:
            exec_finished_list, exec_running_list, exec_queued_list = self.get_execution_status_lists(
                current_combo['combination_id']
            )
            execution_total = current_combo.get('total_sequences') or 0
        else:
            exec_finished_list = exec_running_list = exec_queued_list = []
            execution_total = 0

        # GLOBAL EXECUTIONS (finished via executions; total via soma de combinations.total_sequences)
        global_stats = self.get_global_execution_stats(work_id)
        global_finished = global_stats["Finished"]
        global_total = global_stats["Total"]
        global_progress = global_stats["Progress"]

        return ProgressSummary(
            work_id=work_id,
            tasks={
                "Finished": len(task_finished),
                "Running": running_task_id,
                "Queued": len(task_queued),
                "Total": len(task_finished) + len(task_running) + len(task_queued),
            },
            datasets={
                "Finished": len(dataset_finished),
                "Running": running_dataset_id,
                "Queued": len(dataset_queued),
                "Total": len(dataset_finished) + len(dataset_running) + len(dataset_queued),
            },
            configs={
                "Finished": len(config_finished),
                "Running": running_config_id,
                "Queued": len(config_queued),
                "Total": len(config_finished) + len(config_running) + len(config_queued),
            },
            algorithms={
                "Finished": len(alg_finished),
                "Running": running_alg_id,
                "Queued": len(alg_queued),
                "Total": len(alg_finished) + len(alg_running) + len(alg_queued),
            },
            execution={
                "Finished": len(exec_finished_list),
                "Running": len(exec_running_list),
                "Queued": len(exec_queued_list),
                "Total": execution_total,
            },
            global_execution={
                "Finished": global_finished,
                "Total": global_total,
            },
            global_progress=global_progress,
            current_combination_details=dict(current_combo) if current_combo else None
        )
    
    def get_running_executions_detail(self, work_id: str, limit: int = 20) -> List[ExecutionDetail]:
        """Get detailed information about running executions."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        query = """
            SELECT 
                e.unit_id,
                e.combination_id,
                e.sequencia,
                e.status,
                e.started_at,
                e.finished_at,
                e.objective,
                c.task_id,
                c.dataset_id,
                c.preset_id,
                c.algorithm_id,
                c.mode,
                c.total_sequences,
                -- Get latest progress
                (SELECT ep.progress 
                 FROM execution_progress ep 
                 WHERE ep.execution_id = e.id 
                 ORDER BY ep.timestamp DESC LIMIT 1) as progress,
                (SELECT ep.message 
                 FROM execution_progress ep 
                 WHERE ep.execution_id = e.id 
                 ORDER BY ep.timestamp DESC LIMIT 1) as progress_message
            FROM executions e
            JOIN combinations c ON e.combination_id = c.id
            WHERE c.work_id = ? AND e.status = 'running'
            ORDER BY e.started_at DESC
            LIMIT ?
        """
        
        cursor.execute(query, (work_id, limit))
        rows = cursor.fetchall()
        
        results = []
        for row in rows:
            results.append(ExecutionDetail(
                unit_id=row['unit_id'],
                combination_id=row['combination_id'],
                sequencia=row['sequencia'] or 0,
                status=row['status'],
                progress=row['progress'] or 0.0,
                progress_message=row['progress_message'],
                started_at=row['started_at'],
                finished_at=row['finished_at'],
                objective=row['objective'],
                task_id=row['task_id'],
                dataset_id=row['dataset_id'],
                preset_id=row['preset_id'],
                algorithm_id=row['algorithm_id'],
                mode=row['mode'],
                total_sequences=row['total_sequences'] or 0
            ))
        
        return results
    
    def get_error_summary(self, work_id: str, limit: int = 10) -> List[ErrorSummary]:
        """Get recent errors for a work item."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        query = """
            SELECT 
                e.unit_id,
                e.status as error_type,
                COALESCE(json_extract(e.result_json, '$.error_message'), 'Unknown error') as error_message,
                e.finished_at as timestamp
            FROM executions e
            JOIN combinations c ON e.combination_id = c.id
            WHERE c.work_id = ? AND e.status IN ('failed', 'error')
            ORDER BY e.finished_at DESC
            LIMIT ?
        """
        
        cursor.execute(query, (work_id, limit))
        rows = cursor.fetchall()
        
        results = []
        for row in rows:
            results.append(ErrorSummary(
                unit_id=row['unit_id'],
                error_type=row['error_type'],
                error_message=row['error_message'],
                timestamp=row['timestamp'] or 0.0
            ))
        
        return results
    
    def get_execution_warnings(self, work_id: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent warnings from events table."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        query = """
            SELECT 
                event_type,
                event_category,
                entity_data_json,
                timestamp
            FROM events
            WHERE work_id = ? AND event_type = 'warning'
            ORDER BY timestamp DESC
            LIMIT ?
        """
        
        cursor.execute(query, (work_id, limit))
        rows = cursor.fetchall()
        
        results = []
        for row in rows:
            entity_data = json.loads(row['entity_data_json']) if row['entity_data_json'] else {}
            results.append({
                'event_type': row['event_type'],
                'event_category': row['event_category'],
                'message': entity_data.get('message', 'Unknown warning'),
                'unit_id': entity_data.get('unit_id'),
                'combination_id': entity_data.get('combination_id'),
                'timestamp': row['timestamp']
            })
        
        return results
    
    # === STATUS QUERIES ===
    
    def get_work_status(self, work_id: str) -> Optional[str]:
        """Get current work status."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("SELECT status FROM work WHERE id = ?", (work_id,))
        row = cursor.fetchone()
        return row['status'] if row else None
    
    def get_combination_status_counts(self, work_id: str) -> Dict[str, int]:
        """Get counts of combinations by status."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT status, COUNT(*) as count
            FROM combinations
            WHERE work_id = ?
            GROUP BY status
        """, (work_id,))
        
        return {row['status']: row['count'] for row in cursor.fetchall()}
    
    def get_execution_status_counts(self, work_id: str) -> Dict[str, int]:
        """Get counts of executions by status."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT e.status, COUNT(*) as count
            FROM executions e
            JOIN combinations c ON e.combination_id = c.id
            WHERE c.work_id = ?
            GROUP BY e.status
        """, (work_id,))
        
        return {row['status']: row['count'] for row in cursor.fetchall()}
    
    # === UTILITY QUERIES ===
    
    def work_exists(self, work_id: str) -> bool:
        """Check if work exists."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("SELECT 1 FROM work WHERE id = ? LIMIT 1", (work_id,))
        return cursor.fetchone() is not None
    
    def get_work_info(self, work_id: str) -> Optional[Dict[str, Any]]:
        """Get basic work information."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT id, status, created_at, updated_at, output_path, error
            FROM work WHERE id = ?
        """, (work_id,))
        
        row = cursor.fetchone()
        return dict(row) if row else None
