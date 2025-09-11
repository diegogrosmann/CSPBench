"""Advanced queries mixin for complex database operations."""

import time
from typing import Any, Dict, List, Optional, Tuple, Union
from sqlalchemy import and_, desc, func, case, text
from dataclasses import dataclass


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


class QueriesMixin:
    """Advanced queries for complex database operations."""

    def get_running_combination(self, work_id: str) -> Optional[Dict[str, Any]]:
        """Get running combination."""
        from ..models import Combination
        
        with self.session_scope() as session:
            combination = session.query(Combination).filter(
                Combination.work_id == work_id,
                Combination.status == 'running'
            ).order_by(
                Combination.task_id,
                Combination.dataset_id,
                Combination.preset_id,
                Combination.algorithm_id
            ).first()
            
            if combination:
                return {
                    'id': combination.id,
                    'combination_id': combination.id,
                    'task_id': combination.task_id,
                    'dataset_id': combination.dataset_id,
                    'preset_id': combination.preset_id,
                    'algorithm_id': combination.algorithm_id,
                    'total_sequences': combination.total_sequences
                }
            return None

    def get_last_paused_or_incomplete_combination(self, work_id: str) -> Optional[Dict[str, Any]]:
        """Get the last paused or incomplete combination for context when no running combination exists."""
        from ..models import Combination
        
        with self.session_scope() as session:
            # Try paused first
            combination = session.query(Combination).filter(
                Combination.work_id == work_id,
                Combination.status == 'paused'
            ).order_by(desc(Combination.id)).first()
            
            if combination:
                return self._combination_to_dict(combination)
            
            # Try queued
            combination = session.query(Combination).filter(
                Combination.work_id == work_id,
                Combination.status == 'queued'
            ).order_by(Combination.id).first()
            
            if combination:
                return self._combination_to_dict(combination)
            
            # Try failed/error
            combination = session.query(Combination).filter(
                Combination.work_id == work_id,
                Combination.status.in_(['failed', 'error'])
            ).order_by(desc(Combination.id)).first()
            
            if combination:
                return self._combination_to_dict(combination)
            
            # Try cancelled
            combination = session.query(Combination).filter(
                Combination.work_id == work_id,
                Combination.status == 'canceled'
            ).order_by(desc(Combination.id)).first()
            
            if combination:
                return self._combination_to_dict(combination)
                
            return None

    def _combination_to_dict(self, combination) -> Dict[str, Any]:
        """Convert combination to dict format expected by queries."""
        return {
            'id': combination.id,
            'combination_id': combination.id,
            'task_id': combination.task_id,
            'dataset_id': combination.dataset_id,
            'preset_id': combination.preset_id,
            'algorithm_id': combination.algorithm_id,
            'total_sequences': combination.total_sequences
        }

    def get_task_status_lists(self, work_id: str) -> Tuple[List[str], List[str], List[str]]:
        """Get task status lists (finished, running, queued)."""
        from ..models import Combination
        
        with self.session_scope() as session:
            # Get task status aggregation
            result = session.query(
                Combination.task_id,
                func.sum(case((Combination.status.in_(['completed', 'failed', 'error']), 1), else_=0)).label('completed_cnt'),
                func.sum(case((Combination.status == 'queued', 1), else_=0)).label('queued_cnt'),
                func.count().label('total_cnt')
            ).filter(
                Combination.work_id == work_id
            ).group_by(Combination.task_id).all()
            
            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.task_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.task_id)
                else:
                    running.append(row.task_id)
                    
            return finished, running, queued

    def get_dataset_status_lists(self, work_id: str, task_id: str) -> Tuple[List[str], List[str], List[str]]:
        """Get dataset status lists within a task."""
        from ..models import Combination
        
        with self.session_scope() as session:
            result = session.query(
                Combination.dataset_id,
                func.sum(case((Combination.status.in_(['completed', 'failed', 'error']), 1), else_=0)).label('completed_cnt'),
                func.sum(case((Combination.status == 'queued', 1), else_=0)).label('queued_cnt'),
                func.count().label('total_cnt')
            ).filter(
                Combination.work_id == work_id,
                Combination.task_id == task_id
            ).group_by(Combination.dataset_id).all()
            
            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.dataset_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.dataset_id)
                else:
                    running.append(row.dataset_id)
                    
            return finished, running, queued

    def get_config_status_lists(self, work_id: str, task_id: str, dataset_id: str) -> Tuple[List[str], List[str], List[str]]:
        """Get config (preset) status lists within a dataset."""
        from ..models import Combination
        
        with self.session_scope() as session:
            result = session.query(
                Combination.preset_id,
                func.sum(case((Combination.status.in_(['completed', 'failed', 'error']), 1), else_=0)).label('completed_cnt'),
                func.sum(case((Combination.status == 'queued', 1), else_=0)).label('queued_cnt'),
                func.count().label('total_cnt')
            ).filter(
                Combination.work_id == work_id,
                Combination.task_id == task_id,
                Combination.dataset_id == dataset_id
            ).group_by(Combination.preset_id).all()
            
            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.preset_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.preset_id)
                else:
                    running.append(row.preset_id)
                    
            return finished, running, queued

    def get_algorithm_status_lists(self, work_id: str, task_id: str, dataset_id: str, preset_id: str) -> Tuple[List[str], List[str], List[str]]:
        """Get algorithm status lists within a config."""
        from ..models import Combination
        
        with self.session_scope() as session:
            result = session.query(
                Combination.algorithm_id,
                func.sum(case((Combination.status.in_(['completed', 'failed', 'error']), 1), else_=0)).label('completed_cnt'),
                func.sum(case((Combination.status == 'queued', 1), else_=0)).label('queued_cnt'),
                func.count().label('total_cnt')
            ).filter(
                Combination.work_id == work_id,
                Combination.task_id == task_id,
                Combination.dataset_id == dataset_id,
                Combination.preset_id == preset_id
            ).group_by(Combination.algorithm_id).all()
            
            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.algorithm_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.algorithm_id)
                else:
                    running.append(row.algorithm_id)
                    
            return finished, running, queued

    def get_execution_status_lists(self, combination_id: int) -> Tuple[List[int], List[int], List[int]]:
        """Get execution status lists for a combination."""
        from ..models import Execution
        
        with self.session_scope() as session:
            executions = session.query(Execution).filter(
                Execution.combination_id == combination_id
            ).all()
            
            finished, running, queued = [], [], []
            for execution in executions:
                if execution.status in ('completed', 'failed', 'error'):
                    finished.append(execution.id)
                elif execution.status == 'running':
                    running.append(execution.id)
                elif execution.status == 'queued':
                    queued.append(execution.id)
                    
            return finished, running, queued

    def get_global_execution_stats(self, work_id: str) -> Dict[str, Union[int, float]]:
        """Get global execution statistics."""
        from ..models import Combination, Execution
        
        with self.session_scope() as session:
            # Total sequences from combinations
            total_sequences = session.query(
                func.coalesce(func.sum(Combination.total_sequences), 0)
            ).filter(Combination.work_id == work_id).scalar() or 0
            
            # Finished executions
            finished = session.query(func.count(Execution.id)).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Combination.work_id == work_id,
                Execution.status.in_(['completed', 'failed', 'error'])
            ).scalar() or 0
            
            progress = finished / total_sequences if total_sequences else 0.0
            
            return {
                "Finished": finished,
                "Total": total_sequences,
                "Progress": progress
            }

    def get_work_progress_summary(self, work_id: str) -> Optional[ProgressSummary]:
        """Get comprehensive progress summary for a work item."""
        work = self.work_get(work_id)
        if not work:
            return None

        current_combo = self.get_running_combination(work_id)
        if not current_combo:
            current_combo = self.get_last_paused_or_incomplete_combination(work_id)

        # TASKS
        task_finished, task_running, task_queued = self.get_task_status_lists(work_id)
        running_task_id = current_combo['task_id'] if current_combo else ""

        # DATASETS (only for current task)
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

        # GLOBAL EXECUTIONS
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

    def get_combination_executions_detail(self, combination_id: int) -> List[ExecutionDetail]:
        """Get detailed information about all executions for a specific combination."""
        from ..models import Execution, Combination, ExecutionProgress
        
        with self.session_scope() as session:
            # Query executions with combination data
            query = session.query(Execution, Combination).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Execution.combination_id == combination_id
            ).order_by(Execution.sequencia)
            
            results = []
            for execution, combination in query.all():
                # Get latest progress
                latest_progress = session.query(ExecutionProgress).filter(
                    ExecutionProgress.execution_id == execution.id
                ).order_by(desc(ExecutionProgress.timestamp)).first()
                
                progress = latest_progress.progress if latest_progress else 0.0
                progress_message = latest_progress.message if latest_progress else None
                
                # Format progress message based on status
                if execution.status == 'completed':
                    progress_message = f"Objective: {execution.objective}"
                elif execution.status in ('failed', 'error') and execution.result_json:
                    progress_message = execution.result_json.get('error_message', 'Unknown error')
                
                results.append(ExecutionDetail(
                    unit_id=execution.unit_id,
                    combination_id=execution.combination_id,
                    sequencia=execution.sequencia,
                    status=execution.status,
                    progress=progress,
                    progress_message=progress_message,
                    started_at=execution.started_at,
                    finished_at=execution.finished_at,
                    objective=execution.objective,
                    task_id=combination.task_id,
                    dataset_id=combination.dataset_id,
                    preset_id=combination.preset_id,
                    algorithm_id=combination.algorithm_id,
                    mode=combination.mode,
                    total_sequences=combination.total_sequences,
                ))
        
        return results

    def get_running_executions_detail(self, work_id: str, limit: int = 20) -> List[ExecutionDetail]:
        """Get detailed information about running executions."""
        from ..models import Execution, Combination, ExecutionProgress
        
        with self.session_scope() as session:
            query = session.query(Execution, Combination).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Combination.work_id == work_id,
                Execution.status == 'running'
            ).order_by(desc(Execution.started_at)).limit(limit)
            
            results = []
            for execution, combination in query.all():
                # Get latest progress
                latest_progress = session.query(ExecutionProgress).filter(
                    ExecutionProgress.execution_id == execution.id
                ).order_by(desc(ExecutionProgress.timestamp)).first()
                
                progress = latest_progress.progress if latest_progress else 0.0
                progress_message = latest_progress.message if latest_progress else None
                
                results.append(ExecutionDetail(
                    unit_id=execution.unit_id,
                    combination_id=execution.combination_id,
                    sequencia=execution.sequencia or 0,
                    status=execution.status,
                    progress=progress,
                    progress_message=progress_message,
                    started_at=execution.started_at,
                    finished_at=execution.finished_at,
                    objective=execution.objective,
                    task_id=combination.task_id,
                    dataset_id=combination.dataset_id,
                    preset_id=combination.preset_id,
                    algorithm_id=combination.algorithm_id,
                    mode=combination.mode,
                    total_sequences=combination.total_sequences or 0
                ))
        
        return results

    def get_error_summary(self, work_id: str, limit: int = 10) -> List[ErrorSummary]:
        """Get recent errors for a work item."""
        from ..models import Execution, Combination
        
        with self.session_scope() as session:
            query = session.query(Execution).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Combination.work_id == work_id,
                Execution.status.in_(['failed', 'error'])
            ).order_by(desc(Execution.finished_at)).limit(limit)
            
            results = []
            for execution in query.all():
                error_message = "Unknown error"
                if execution.result_json and 'error_message' in execution.result_json:
                    error_message = execution.result_json['error_message']
                
                results.append(ErrorSummary(
                    unit_id=execution.unit_id,
                    error_type=execution.status,
                    error_message=error_message,
                    timestamp=execution.finished_at or 0.0,
                ))
        
        return results

    def get_execution_warnings(self, work_id: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent warnings from events table."""
        from ..models import Event
        
        with self.session_scope() as session:
            query = session.query(Event).filter(
                Event.work_id == work_id,
                Event.event_type == 'warning'
            ).order_by(desc(Event.timestamp)).limit(limit)
            
            results = []
            for event in query.all():
                entity_data = event.entity_data_json or {}
                results.append({
                    'event_type': event.event_type,
                    'event_category': event.event_category,
                    'message': entity_data.get('message', 'Unknown warning'),
                    'unit_id': entity_data.get('unit_id'),
                    'combination_id': entity_data.get('combination_id'),
                    'timestamp': event.timestamp
                })
        
        return results

    def get_events(self, work_id: str, limit: int = 100, event_types: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """Get events for a work item."""
        from ..models import Event
        
        results = []
        with self.session_scope() as session:
            query = session.query(Event).filter(
                Event.work_id == work_id
            )
            
            if event_types:
                query = query.filter(Event.event_type.in_(event_types))
                
            query = query.order_by(desc(Event.timestamp)).limit(limit)
            
            for event in query.all():
                entity_data = event.entity_data_json or {}
                results.append({
                    'id': event.id,
                    'work_id': event.work_id,
                    'event_type': event.event_type,
                    'event_category': event.event_category,
                    'entity_data': entity_data,
                    'timestamp': event.timestamp,
                    'message': entity_data.get('message', '')
                })
        
        return results

    def get_combination_status_counts(self, work_id: str) -> Dict[str, int]:
        """Get counts of combinations by status."""
        from ..models import Combination
        
        with self.session_scope() as session:
            result = session.query(
                Combination.status,
                func.count().label('count')
            ).filter(
                Combination.work_id == work_id
            ).group_by(Combination.status).all()
            
            return {row.status: row.count for row in result}

    def get_execution_status_counts(self, work_id: str) -> Dict[str, int]:
        """Get counts of executions by status."""
        from ..models import Execution, Combination
        
        with self.session_scope() as session:
            result = session.query(
                Execution.status,
                func.count().label('count')
            ).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Combination.work_id == work_id
            ).group_by(Execution.status).all()
            
            return {row.status: row.count for row in result}

    def get_execution_stats(self, work_id: str) -> Dict[str, Any]:
        """Return aggregated execution stats across all combinations."""
        from ..models import Execution, Combination
        
        with self.session_scope() as session:
            # Status counts at execution level
            status_counts = session.query(
                Execution.status,
                func.count().label('cnt')
            ).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Combination.work_id == work_id
            ).group_by(Execution.status).all()
            
            status_dict = {row.status: row.cnt for row in status_counts}
            
            # Total sequences across combinations
            total_sequences = session.query(
                func.coalesce(func.sum(Combination.total_sequences), 0)
            ).filter(Combination.work_id == work_id).scalar() or 0
            
            completed = status_dict.get('completed', 0)
            running = status_dict.get('running', 0)
            queued = status_dict.get('queued', 0)
            failed = status_dict.get('failed', 0) + status_dict.get('error', 0)
            
            return {
                'completed': completed,
                'running': running,
                'queued': queued,
                'failed': failed,
                'total_sequences': total_sequences,
                'raw_status_counts': status_dict,
            }

    # -------------- EXPORT SPECIFIC QUERIES --------------
    def get_work_export_data(self, work_id: str) -> Dict[str, Any]:
        """Get work data for export."""
        from ..models import Work
        
        with self.session_scope() as session:
            work = session.query(Work).filter(Work.id == work_id).first()
            return work.to_dict() if work else {}

    def list_combinations(self, work_id: str) -> List[Dict[str, Any]]:
        """List all combinations for a work with combination_id field."""
        from ..models import Combination
        
        with self.session_scope() as session:
            combinations = session.query(Combination).filter(
                Combination.work_id == work_id
            ).order_by(Combination.id).all()
            
            result = []
            for combo in combinations:
                combo_dict = combo.to_dict()
                # Ensure combination_id field is present (some code expects this field name)
                combo_dict['combination_id'] = combo_dict.get('id')
                result.append(combo_dict)
            
            return result

    def get_combinations_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all combinations for a work for export."""
        from ..models import Combination
        
        with self.session_scope() as session:
            combinations = session.query(Combination).filter(
                Combination.work_id == work_id
            ).all()
            return [combo.to_dict() for combo in combinations]

    def get_executions_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all executions for a work for export."""
        from ..models import Execution, Combination
        
        with self.session_scope() as session:
            executions = session.query(Execution).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(Combination.work_id == work_id).all()
            return [exec.to_dict() for exec in executions]

    def get_execution_progress_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all execution progress for a work for export."""
        from ..models import ExecutionProgress, Execution, Combination
        
        with self.session_scope() as session:
            progress_records = session.query(ExecutionProgress).join(
                Execution, ExecutionProgress.execution_id == Execution.id
            ).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(Combination.work_id == work_id).all()
            return [progress.to_dict() for progress in progress_records]

    def get_events_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all events for a work for export."""
        from ..models import Event
        
        with self.session_scope() as session:
            # Check if work_id column exists
            try:
                events = session.query(Event).filter(Event.work_id == work_id).all()
                return [event.to_dict() for event in events]
            except Exception:
                # Fallback: return all events if work_id column doesn't exist
                events = session.query(Event).all()
                return [event.to_dict() for event in events]

    def get_datasets_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get datasets used by a work for export."""
        from ..models import Dataset
        
        with self.session_scope() as session:
            datasets = session.query(Dataset).filter(
                Dataset.work_id == work_id
            ).all()
            return [dataset.to_dict() for dataset in datasets]

    def get_dataset_sequences_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get dataset sequences used by a work for export."""
        from ..models import DatasetSequence, Dataset
        
        with self.session_scope() as session:
            sequences = session.query(DatasetSequence).join(
                Dataset, DatasetSequence.dataset_id == Dataset.id
            ).filter(Dataset.work_id == work_id).all()
            return [seq.to_dict() for seq in sequences]

    def get_optimization_executions_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get optimization executions for a work for export."""
        from ..models import Execution, Combination
        
        with self.session_scope() as session:
            executions = session.query(Execution).join(
                Combination, Execution.combination_id == Combination.id
            ).filter(
                Combination.work_id == work_id,
                Execution.unit_id.like('optimization:%')
            ).order_by(Execution.unit_id, Execution.sequencia).all()
            return [exec.to_dict() for exec in executions]

    def get_sensitivity_events_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get sensitivity analysis events for a work for export."""
        from ..models import Event
        
        with self.session_scope() as session:
            try:
                # Try with work_id filter first - look for progress events with sensitivity_analysis unit_id
                events = session.query(Event).filter(
                    Event.event_type == 'progress',
                    Event.work_id == work_id
                ).filter(
                    Event.entity_data_json.like('%"unit_id": "sensitivity_analysis"%')
                ).order_by(Event.timestamp).all()
                return [event.to_dict() for event in events]
            except Exception:
                # Fallback: return all sensitivity events if work_id column doesn't exist
                events = session.query(Event).filter(
                    Event.event_type == 'progress'
                ).filter(
                    Event.entity_data_json.like('%"unit_id": "sensitivity_analysis"%')
                ).order_by(Event.timestamp).all()
                return [event.to_dict() for event in events]