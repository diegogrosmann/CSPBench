"""
API routes for real-time monitoring and execution management
"""

from fastapi import APIRouter, HTTPException, Depends
from typing import List, Dict, Any, Optional
from datetime import datetime, timedelta
import os
import glob
from pathlib import Path

from src.application.work.global_manager import get_global_work_manager
from src.application.work.manager import WorkManager
from src.presentation.web.websocket_manager import connection_manager
from src.infrastructure.logging_config import get_logger

logger = get_logger(__name__)
router = APIRouter()


@router.get("/api/monitoring/active-executions")
async def get_active_executions(work_manager: WorkManager = Depends(get_global_work_manager)):
    """Get list of active executions"""
    try:
        # Get active work items from work manager
        active_works = []
        
        # For now, we'll check the outputs directory for recent executions
        outputs_dir = Path("outputs")
        if outputs_dir.exists():
            # Look for recent directories (created in last 24 hours)
            recent_threshold = datetime.now() - timedelta(hours=24)
            
            for work_dir in outputs_dir.iterdir():
                if work_dir.is_dir():
                    try:
                        # Get directory creation time
                        stat = work_dir.stat()
                        created_time = datetime.fromtimestamp(stat.st_ctime)
                        
                        if created_time > recent_threshold:
                            # Check for various status indicators
                            log_files = list(work_dir.glob("*.log"))
                            state_files = list(work_dir.glob("state.db"))
                            output_files = list(work_dir.glob("*.json")) + list(work_dir.glob("*.csv"))
                            
                            # Determine status based on files and timing
                            status = "unknown"
                            now = datetime.now()
                            time_since_created = now - created_time
                            time_since_modified = now - datetime.fromtimestamp(stat.st_mtime)
                            
                            # Check if execution completed successfully
                            if state_files and output_files:
                                status = "finished"
                            elif state_files:
                                # Check creation/modification time to determine if it just finished
                                if time_since_created < timedelta(minutes=2):
                                    # Very recent execution, likely just finished
                                    status = "finished" 
                                elif time_since_modified < timedelta(minutes=5):
                                    # Recently modified, might be running
                                    status = "running"
                                else:
                                    # Older execution
                                    status = "finished"
                            
                            # Check main log file for more detailed status
                            latest_log_content = ""
                            main_log_path = Path("logs/cspbench.log")
                            if main_log_path.exists():
                                try:
                                    with open(main_log_path, 'r') as f:
                                        lines = f.readlines()
                                        # Filter lines containing this work_id
                                        relevant_lines = [line for line in lines if work_dir.name in line]
                                        latest_log_content = ''.join(relevant_lines[-10:]) if relevant_lines else ""
                                except Exception:
                                    pass
                            
                            # Check log content for completion indicators
                            if latest_log_content:
                                log_lower = latest_log_content.lower()
                                
                                # Check for definitive completion messages first
                                if any(phrase in log_lower for phrase in ["completed successfully", "pipeline concluído", "finalizado", "finished", "marcado como finalizado"]):
                                    status = "finished"
                                elif any(phrase in log_lower for phrase in ["error", "failed", "falhou", "exception"]):
                                    status = "failed"
                                elif any(phrase in log_lower for phrase in ["iniciando", "starting", "running"]):
                                    # Only consider it running if very recent
                                    if time_since_created < timedelta(minutes=2):
                                        # Check if there are any completion messages after the start
                                        if any(phrase in log_lower for phrase in ["concluído", "completed", "finalizado", "marcado como finalizado"]):
                                            status = "finished"
                                        else:
                                            status = "running"
                                    else:
                                        status = "finished"  # Old execution, assume finished
                            
                            # Check log files for more detailed status
                            if log_files:
                                latest_log = max(log_files, key=lambda f: f.stat().st_mtime)
                                modified_time = datetime.fromtimestamp(latest_log.stat().st_mtime)
                                
                                # If log was modified very recently, consider it running
                                if datetime.now() - modified_time < timedelta(minutes=2):
                                    status = "running"
                                else:
                                    # Read log content to check completion
                                    try:
                                        with open(latest_log, 'r') as f:
                                            log_content = f.read().lower()
                                            if "completed" in log_content or "finished" in log_content:
                                                status = "finished"
                                            elif "error" in log_content or "failed" in log_content:
                                                status = "failed"
                                    except Exception:
                                        pass
                            
                            active_works.append({
                                "work_id": work_dir.name,
                                "status": status,
                                "created_at": created_time.isoformat(),
                                "last_modified": datetime.fromtimestamp(stat.st_mtime).isoformat(),
                                "output_path": str(work_dir),
                                "log_files": [str(f) for f in log_files]
                            })
                            
                    except Exception as e:
                        logger.warning(f"Error processing work directory {work_dir}: {e}")
                        continue
        
        # Sort by creation time (newest first)
        active_works.sort(key=lambda x: x["created_at"], reverse=True)
        
        return {
            "executions": active_works,
            "total": len(active_works),
            "timestamp": datetime.now().isoformat()
        }
        
    except Exception as e:
        logger.error(f"Error getting active executions: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/execution/{work_id}/detailed")
async def get_execution_detailed(work_id: str, work_manager: WorkManager = Depends(get_global_work_manager)):
    """Get detailed execution information including task progress and trials"""
    try:
        work_dir = Path("outputs") / work_id
        if not work_dir.exists():
            raise HTTPException(status_code=404, detail=f"Execution {work_id} not found")
        
        # Get basic execution info
        basic_info = await get_execution_details(work_id, work_manager)
        
        # Get pipeline state from database
        state_db_path = work_dir / "state.db"
        pipeline_state = {}
        trial_data = []
        
        if state_db_path.exists():
            try:
                import sqlite3
                conn = sqlite3.connect(str(state_db_path))
                conn.row_factory = sqlite3.Row
                cursor = conn.cursor()
                
                # Get pipeline information
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
                tables = [row[0] for row in cursor.fetchall()]
                
                if 'pipeline_state' in tables:
                    cursor.execute("SELECT * FROM pipeline_state ORDER BY id DESC LIMIT 1")
                    state_row = cursor.fetchone()
                    if state_row:
                        pipeline_state = dict(state_row)
                
                # Get task progress from work table
                if 'work' in tables:
                    cursor.execute("SELECT * FROM work ORDER BY created_at")
                    work_rows = [dict(row) for row in cursor.fetchall()]
                    pipeline_state['work_history'] = work_rows
                    
                    # Count tasks by status
                    status_counts = {}
                    for work in work_rows:
                        status = work.get('status', 'unknown')
                        status_counts[status] = status_counts.get(status, 0) + 1
                    pipeline_state['status_counts'] = status_counts
                
                # Get configuration data
                if 'config' in tables:
                    cursor.execute("SELECT * FROM config WHERE id = 'config'")
                    config_row = cursor.fetchone()
                    if config_row:
                        import json
                        try:
                            config_data = json.loads(config_row['json'])
                            pipeline_state['config'] = config_data
                            
                            # Extract task information from config
                            tasks = config_data.get('tasks', {}).get('items', [])
                            pipeline_state['tasks'] = tasks
                        except Exception as e:
                            logger.warning(f"Error parsing config JSON: {e}")
                
                # Get dataset information
                if 'datasets' in tables:
                    cursor.execute("SELECT * FROM datasets")
                    datasets = [dict(row) for row in cursor.fetchall()]
                    pipeline_state['datasets'] = datasets
                
                # Create trial/repetition data from work table and config
                trial_data = []
                if 'work' in tables and pipeline_state.get('config'):
                    config = pipeline_state['config']
                    tasks = config.get('tasks', {}).get('items', [])
                    
                    for i, task in enumerate(tasks):
                        repetitions = task.get('repetitions', 1)
                        algorithms = task.get('algorithms', [])
                        datasets = task.get('datasets', [])
                        
                        # Create trial entries for each combination
                        trial_id_counter = 1
                        for rep in range(repetitions):
                            for algo in algorithms if algorithms else [{'id': 'default', 'name': 'Default'}]:
                                for dataset in datasets:
                                    trial_data.append({
                                        'trial_id': f"trial_{trial_id_counter:03d}",
                                        'task_id': task.get('id', f'task_{i+1}'),
                                        'task_name': task.get('name', f'Task {i+1}'),
                                        'algorithm_name': algo.get('name', 'Default'),
                                        'algorithm_id': algo.get('id', 'default'),
                                        'dataset_id': dataset.get('id', 'unknown'),
                                        'dataset_name': dataset.get('name', 'Unknown'),
                                        'repetition_number': rep + 1,
                                        'status': 'completed' if work_rows and work_rows[-1].get('status') == 'finished' else 'queued',
                                        'start_time': work_rows[0].get('created_at') if work_rows else None,
                                        'end_time': work_rows[-1].get('updated_at') if work_rows and work_rows[-1].get('status') == 'finished' else None,
                                        'duration': None,
                                        'result_path': None,
                                        'error_message': work_rows[-1].get('error') if work_rows else None,
                                        'config': algo
                                    })
                                    trial_id_counter += 1
                
                # If no trials from config, create from executions table
                if not trial_data and 'executions' in tables:
                    cursor.execute("""
                        SELECT trial_id, algorithm_name, dataset_id, status, 
                               start_time, end_time, duration, 
                               result_path, error_message, config
                        FROM executions 
                        ORDER BY start_time DESC
                    """)
                    trial_data = [dict(row) for row in cursor.fetchall()]
                
                conn.close()
                
            except Exception as e:
                logger.warning(f"Error reading state database: {e}")
        
        # Get log progress indicators
        main_log_path = Path("logs/cspbench.log")
        progress_info = {}
        if main_log_path.exists():
            try:
                with open(main_log_path, 'r') as f:
                    lines = f.readlines()
                    # Filter lines for this work_id
                    relevant_lines = [line for line in lines if work_id in line]
                    
                    # Extract progress information
                    progress_info = {
                        'current_task': None,
                        'current_algorithm': None,
                        'current_dataset': None,
                        'current_trial': None,
                        'last_activity': None
                    }
                    
                    for line in reversed(relevant_lines[-50:]):  # Check last 50 relevant lines
                        if 'Iniciando experiment:' in line:
                            progress_info['current_task'] = line.split('Iniciando experiment:')[1].strip()
                        elif 'algoritmo:' in line or 'Algorithm:' in line:
                            progress_info['current_algorithm'] = line.split(':')[-1].strip()
                        elif 'dataset:' in line or 'Dataset:' in line:
                            progress_info['current_dataset'] = line.split(':')[-1].strip()
                        elif 'trial' in line.lower() or 'repetição' in line.lower():
                            progress_info['current_trial'] = line.strip()
                        
                        if not progress_info['last_activity']:
                            progress_info['last_activity'] = line.strip()
                            
            except Exception as e:
                logger.warning(f"Error reading log file: {e}")
        
        # Debug: log summary of detailed execution info being returned
        logger.debug(
            "get_execution_detailed: work_id=%s, trials=%d, pipeline_tasks=%d, progress=%s",
            work_id,
            len(trial_data),
            len(pipeline_state.get('tasks', [])),
            {k: v for k, v in progress_info.items()}
        )

        return {
            **basic_info,
            "pipeline_state": pipeline_state,
            "trials": trial_data,
            "progress": progress_info,
            "detailed_metrics": {
                "total_trials": len(trial_data),
                "completed_trials": len([t for t in trial_data if t.get('status') == 'completed']),
                "running_trials": len([t for t in trial_data if t.get('status') == 'running']),
                "failed_trials": len([t for t in trial_data if t.get('status') == 'failed']),
                "total_tasks": len(pipeline_state.get('tasks', [])),
                "completed_tasks": len([t for t in pipeline_state.get('tasks', []) if t.get('status') == 'completed']),
                "running_tasks": len([t for t in pipeline_state.get('tasks', []) if t.get('status') == 'running']),
                "work_status_counts": pipeline_state.get('status_counts', {}),
                "config_summary": {
                    "batch_name": pipeline_state.get('config', {}).get('metadata', {}).get('name', 'Unknown'),
                    "total_datasets": len(pipeline_state.get('datasets', [])),
                    "algorithms_configured": len(pipeline_state.get('config', {}).get('tasks', {}).get('items', [{}])[0].get('algorithms', [])) if pipeline_state.get('config', {}).get('tasks', {}).get('items') else 0
                }
            }
        }
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting detailed execution info: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/execution/{work_id}/executions")
async def get_execution_executions(work_id: str, status: Optional[str] = None):
    """Get execution/repetition data for a work with optional status filter"""
    try:
        work_dir = Path("outputs") / work_id
        state_db_path = work_dir / "state.db"
        
        if not state_db_path.exists():
            return {"executions": [], "total": 0}
        
        import sqlite3
        conn = sqlite3.connect(str(state_db_path))
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        
        # Get execution data
        execution_data = []
        
        # Get configuration and create execution data
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]
        
        if 'config' in tables and 'work' in tables:
            # Get work history
            cursor.execute("SELECT * FROM work ORDER BY created_at")
            work_rows = [dict(row) for row in cursor.fetchall()]
            
            # Get config
            cursor.execute("SELECT * FROM config WHERE id = 'config'")
            config_row = cursor.fetchone()
            if config_row:
                import json
                try:
                    config_data = json.loads(config_row['json'])
                    tasks = config_data.get('tasks', {}).get('items', [])
                    
                    for i, task in enumerate(tasks):
                        repetitions = task.get('repetitions', 1)
                        algorithms = task.get('algorithms', [])
                        datasets = task.get('datasets', [])
                        
                        # Create execution entries for each combination
                        execution_id_counter = 1
                        for rep in range(repetitions):
                            for algo in algorithms if algorithms else [{'id': 'default', 'name': 'Default'}]:
                                for dataset in datasets:
                                    execution_status = 'completed' if work_rows and work_rows[-1].get('status') == 'finished' else 'queued'
                                    
                                    execution = {
                                        'execution_id': f"exec_{execution_id_counter:03d}",
                                        'task_id': task.get('id', f'task_{i+1}'),
                                        'task_name': task.get('name', f'Task {i+1}'),
                                        'algorithm_name': algo.get('name', 'Default'),
                                        'algorithm_id': algo.get('id', 'default'),
                                        'dataset_id': dataset.get('id', 'unknown'),
                                        'dataset_name': dataset.get('name', 'Unknown'),
                                        'repetition_number': rep + 1,
                                        'status': execution_status,
                                        'start_time': work_rows[0].get('created_at') if work_rows else None,
                                        'end_time': work_rows[-1].get('updated_at') if work_rows and work_rows[-1].get('status') == 'finished' else None,
                                        'duration': None,
                                        'result_path': None,
                                        'error_message': work_rows[-1].get('error') if work_rows else None,
                                        'config': algo
                                    }
                                    
                                    # Apply status filter
                                    if not status or execution['status'] == status:
                                        execution_data.append(execution)
                                    
                                    execution_id_counter += 1
                                    
                except Exception as e:
                    logger.warning(f"Error parsing config: {e}")
        
        # Fallback to executions table if it exists
        if not execution_data and 'executions' in tables:
            query = """
                SELECT rowid as execution_id, algorithm, dataset_id, status, 
                       started_at as start_time, finished_at as end_time, 
                       repetition, trial, sample
                FROM executions
            """
            params = []
            
            if status:
                query += " WHERE status = ?"
                params.append(status)
            
            query += " ORDER BY started_at DESC"
            
            cursor.execute(query, params)
            rows = cursor.fetchall()
            for row in rows:
                execution_dict = dict(row)
                execution_dict['execution_id'] = f"exec_{execution_dict['execution_id']:03d}"
                execution_dict['algorithm_name'] = execution_dict.get('algorithm', 'Unknown')
                execution_dict['dataset_name'] = execution_dict.get('dataset_id', 'Unknown')
                execution_dict['algorithm_id'] = execution_dict.get('algorithm', 'unknown')
                execution_dict['repetition_number'] = execution_dict.get('repetition', 1)
                execution_dict['duration'] = None
                execution_dict['result_path'] = None
                execution_dict['error_message'] = None
                execution_data.append(execution_dict)
        
        conn.close()
        
        # Debug: log executions list summary
        logger.debug(
            "get_execution_executions: work_id=%s, executions_returned=%d, status_filter=%s",
            work_id,
            len(execution_data),
            status
        )

        return {
            "executions": execution_data,
            "total": len(execution_data),
            "status_filter": status
        }
        
    except Exception as e:
        logger.error(f"Error getting executions: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/execution/{work_id}/execution/{execution_id}")
async def get_execution_details(work_id: str, execution_id: str):
    """Get detailed information about a specific execution"""
    try:
        work_dir = Path("outputs") / work_id
        state_db_path = work_dir / "state.db"
        
        if not state_db_path.exists():
            raise HTTPException(status_code=404, detail="Execution state not found")
        
        import sqlite3
        conn = sqlite3.connect(str(state_db_path))
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        
        # Check what tables exist
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]
        
        execution_dict = None
        
        # Try to find execution in executions table
        if 'executions' in tables:
            # Extract execution number from execution_id (e.g., "exec_001" -> 1)
            exec_num = execution_id.replace('exec_', '').lstrip('0') or '1'
            cursor.execute("""
                SELECT rowid, algorithm, dataset_id, status, 
                       started_at, finished_at, repetition, trial, sample
                FROM executions 
                WHERE rowid = ?
            """, (exec_num,))
            execution = cursor.fetchone()
            
            if execution:
                execution_dict = dict(execution)
                execution_dict['execution_id'] = execution_id
                execution_dict['algorithm_name'] = execution_dict.get('algorithm', 'Unknown')
                execution_dict['dataset_name'] = execution_dict.get('dataset_id', 'Unknown')
                execution_dict['start_time'] = execution_dict.get('started_at')
                execution_dict['end_time'] = execution_dict.get('finished_at')
                execution_dict['error_message'] = None
        
        # If no execution found, create mock data
        if not execution_dict:
            execution_dict = {
                "execution_id": execution_id,
                "algorithm_name": "Baseline",
                "dataset_name": "Simple Test",
                "status": "completed",
                "start_time": "2025-01-18 17:30:00",
                "end_time": "2025-01-18 17:35:00",
                "error_message": None,
                "repetition": 1
            }
        
        # Get basic callback data (mock for now)
        callback_data = [
            {
                "callback_type": "progress",
                "timestamp": execution_dict.get('start_time', '2025-01-18 17:30:00'),
                "data": {"progress": 0},
                "message": "Execution started"
            }
        ]
        
        if execution_dict.get('end_time'):
            callback_data.append({
                "callback_type": "completion",
                "timestamp": execution_dict['end_time'],
                "data": {"progress": 100},
                "message": "Execution completed"
            })
        
        conn.close()
        
        return {
            "execution": execution_dict,
            "callbacks": callback_data,
            "callback_count": len(callback_data)
        }
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting execution details: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/execution/{work_id}")
async def get_execution_details(work_id: str, work_manager: WorkManager = Depends(get_global_work_manager)):
    """Get detailed information about a specific execution"""
    try:
        # Check if work directory exists
        work_dir = Path("outputs") / work_id
        if not work_dir.exists():
            raise HTTPException(status_code=404, detail=f"Execution {work_id} not found")
        
        # Get basic info
        stat = work_dir.stat()
        created_time = datetime.fromtimestamp(stat.st_ctime)
        modified_time = datetime.fromtimestamp(stat.st_mtime)
        
        # Get log files and output files
        log_files = list(work_dir.glob("*.log"))
        output_files = list(work_dir.glob("*.json")) + list(work_dir.glob("*.csv"))
        state_files = list(work_dir.glob("state.db"))
        
        # Get latest log content (also check main log file)
        latest_log_content = ""
        main_log_path = Path("logs/cspbench.log")
        
        if log_files:
            latest_log = max(log_files, key=lambda f: f.stat().st_mtime)
            try:
                with open(latest_log, 'r') as f:
                    lines = f.readlines()
                    latest_log_content = ''.join(lines[-50:])
            except Exception as e:
                logger.warning(f"Error reading log file {latest_log}: {e}")
        elif main_log_path.exists():
            # If no local logs, check main log for this work_id
            try:
                with open(main_log_path, 'r') as f:
                    lines = f.readlines()
                    # Filter lines containing this work_id
                    relevant_lines = [line for line in lines if work_id in line]
                    latest_log_content = ''.join(relevant_lines[-20:]) if relevant_lines else "No relevant logs found"
            except Exception as e:
                logger.warning(f"Error reading main log file: {e}")
        
        # Determine status with improved logic
        status = "unknown"
        now = datetime.now()
        time_since_created = now - created_time
        time_since_modified = now - modified_time
        
        # Check if execution completed successfully
        if state_files and output_files:
            status = "finished"
        elif state_files:
            # Check creation/modification time
            if time_since_created < timedelta(minutes=2):
                status = "finished"  # Very recent, likely just completed
            elif time_since_modified < timedelta(minutes=5):
                status = "running"   # Recently modified
            else:
                status = "finished"  # Older execution
        
        # Check log content for completion indicators
        if latest_log_content:
            log_lower = latest_log_content.lower()
            
            # Check for definitive completion messages first
            if any(phrase in log_lower for phrase in ["completed successfully", "pipeline concluído", "finalizado", "finished"]):
                status = "finished"
            elif any(phrase in log_lower for phrase in ["error", "failed", "falhou", "exception"]):
                status = "failed"
            elif any(phrase in log_lower for phrase in ["iniciando", "starting", "running"]):
                # Only consider it running if very recent
                if time_since_created < timedelta(minutes=2):
                    # Check if there are any completion messages after the start
                    if any(phrase in log_lower for phrase in ["concluído", "completed", "finalizado"]):
                        status = "finished"
                    else:
                        status = "running"
                else:
                    status = "finished"  # Old execution, assume finished
        
        # Debug: log execution overview being returned
        logger.debug(
            "get_execution_overview: work_id=%s, status=%s, files=%d, metrics_present=%s",
            work_id,
            status,
            len(log_files) + len(output_files),
            (work_dir.exists())
        )

        return {
            "work_id": work_id,
            "status": status,
            "created_at": created_time.isoformat(),
            "last_modified": modified_time.isoformat(),
            "output_path": str(work_dir),
            "files": {
                "logs": [str(f) for f in log_files],
                "outputs": [str(f) for f in output_files],
                "all": [str(f) for f in work_dir.iterdir() if f.is_file()]
            },
            "latest_log": latest_log_content,
            "metrics": await _get_execution_metrics(work_dir)
        }
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting execution details for {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/execution/{work_id}/logs")
async def get_execution_logs(work_id: str, lines: int = 100):
    """Get log content for a specific execution"""
    try:
        work_dir = Path("outputs") / work_id
        if not work_dir.exists():
            raise HTTPException(status_code=404, detail=f"Execution {work_id} not found")
        
        log_files = list(work_dir.glob("*.log"))
        if not log_files:
            return {"logs": "", "message": "No log files found"}
        
        # Get the most recent log file
        latest_log = max(log_files, key=lambda f: f.stat().st_mtime)
        
        try:
            with open(latest_log, 'r') as f:
                all_lines = f.readlines()
                # Get last N lines
                log_lines = all_lines[-lines:] if len(all_lines) > lines else all_lines
                
            return {
                "work_id": work_id,
                "log_file": str(latest_log),
                "logs": ''.join(log_lines),
                "total_lines": len(all_lines),
                "showing_lines": len(log_lines),
                "timestamp": datetime.now().isoformat()
            }
            
        except Exception as e:
            logger.error(f"Error reading log file {latest_log}: {e}")
            raise HTTPException(status_code=500, detail=f"Error reading log file: {e}")
            
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting logs for {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/stats")
async def get_monitoring_stats():
    """Get monitoring system statistics"""
    try:
        # Get WebSocket connection stats
        ws_stats = await connection_manager.get_connection_stats()
        
        # Get execution stats
        outputs_dir = Path("outputs")
        total_executions = 0
        recent_executions = 0
        
        if outputs_dir.exists():
            all_dirs = [d for d in outputs_dir.iterdir() if d.is_dir()]
            total_executions = len(all_dirs)
            
            # Count recent executions (last 24 hours)
            recent_threshold = datetime.now() - timedelta(hours=24)
            for work_dir in all_dirs:
                try:
                    created_time = datetime.fromtimestamp(work_dir.stat().st_ctime)
                    if created_time > recent_threshold:
                        recent_executions += 1
                except Exception:
                    continue
        
        return {
            "websocket": ws_stats,
            "executions": {
                "total": total_executions,
                "recent_24h": recent_executions
            },
            "timestamp": datetime.now().isoformat()
        }
        
    except Exception as e:
        logger.error(f"Error getting monitoring stats: {e}")
        raise HTTPException(status_code=500, detail=str(e))


async def _get_execution_metrics(work_dir: Path) -> Dict[str, Any]:
    """Extract metrics from execution directory"""
    metrics = {
        "duration": None,
        "file_count": 0,
        "total_size": 0,
        "algorithms_run": [],
        "datasets_processed": []
    }
    
    try:
        # Count files and total size
        for file_path in work_dir.rglob("*"):
            if file_path.is_file():
                metrics["file_count"] += 1
                metrics["total_size"] += file_path.stat().st_size
        
        # Try to extract duration from timestamps and logs
        if work_dir.exists():
            created = datetime.fromtimestamp(work_dir.stat().st_ctime)
            modified = datetime.fromtimestamp(work_dir.stat().st_mtime)
            
            # Try to find actual completion time from logs
            completion_time = None
            
            # Check main log file for completion info
            main_log_path = Path("logs/cspbench.log")
            if main_log_path.exists():
                try:
                    with open(main_log_path, 'r') as f:
                        lines = f.readlines()
                        # Look for completion messages for this work_id
                        for line in reversed(lines[-200:]):  # Check last 200 lines
                            if work_dir.name in line and ("concluído" in line.lower() or "completed" in line.lower()):
                                # Try to extract timestamp from log line
                                # Format: 2025-08-18 16:30:01 - INFO - message
                                if " - " in line:
                                    timestamp_str = line.split(" - ")[0].strip()
                                    try:
                                        # Try different timestamp formats
                                        for fmt in ["%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%S.%f"]:
                                            try:
                                                completion_time = datetime.strptime(timestamp_str, fmt)
                                                break
                                            except ValueError:
                                                continue
                                        if completion_time:
                                            break
                                    except Exception:
                                        pass
                except Exception as e:
                    logger.debug(f"Error reading main log for duration: {e}")
            
            # Calculate duration
            if completion_time and completion_time > created:
                duration = (completion_time - created).total_seconds()
            else:
                # Fallback to file timestamps
                duration = (modified - created).total_seconds()
                # If duration is 0 or very small, estimate based on file count and size
                if duration < 1:
                    # Estimate: 1-5 seconds for small batches
                    if metrics["file_count"] > 1:
                        duration = max(1, metrics["file_count"] * 0.5)
                    else:
                        duration = 1
                        
            metrics["duration"] = round(duration, 1)
        
        # Format file size
        size_mb = metrics["total_size"] / (1024 * 1024)
        metrics["total_size_mb"] = round(size_mb, 2)
        
        # Try to extract algorithm and dataset info from state.db or logs
        state_file = work_dir / "state.db"
        if state_file.exists():
            # This could be enhanced to read SQLite database
            # For now, just indicate that execution has state
            metrics["has_state"] = True
        
    except Exception as e:
        logger.warning(f"Error calculating metrics for {work_dir}: {e}")
    
    return metrics



@router.get("/api/monitoring/execution/{work_id}/monitor")
async def get_work_monitor(work_id: str):
    """Retorna o snapshot do monitor associado a um work_id.

    Estratégia:
    - tenta recuperar monitor registrado no WorkManager (in-memory)
    - se presente e possuir método `get_state`, retorna esse snapshot
    - se ausente, tenta ler `WorkStateStore` (state.db) e retorna metadados
    """
    try:
        from src.application.work.global_manager import get_global_work_manager

        wm = get_global_work_manager()
        mon = None
        try:
            mon = wm.get_registered_monitor(work_id)
        except Exception:
            mon = None

        if mon is not None:
            get_state = getattr(mon, "get_state", None)
            if callable(get_state):
                try:
                    return {"monitor": get_state(), "source": "live"}
                except Exception:
                    logger.debug("Falha ao obter estado do monitor vivo para %s", work_id)

        # Fallback: read state.db
        work_dir = Path("outputs") / work_id
        state_db_path = work_dir / "state.db"
        if state_db_path.exists():
            try:
                import sqlite3

                conn = sqlite3.connect(str(state_db_path))
                conn.row_factory = sqlite3.Row
                cursor = conn.cursor()

                cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
                tables = [r[0] for r in cursor.fetchall()]
                result: dict = {"work_id": work_id, "source": "state_db"}

                if 'work' in tables:
                    cursor.execute("SELECT * FROM work ORDER BY created_at DESC LIMIT 1")
                    row = cursor.fetchone()
                    if row:
                        result['work'] = dict(row)

                if 'pipeline_state' in tables:
                    cursor.execute("SELECT * FROM pipeline_state ORDER BY id DESC LIMIT 1")
                    row = cursor.fetchone()
                    if row:
                        result['pipeline_state'] = dict(row)

                conn.close()
                return result
            except Exception as e:
                logger.warning(f"Error reading state.db for monitor fallback: {e}")

        raise HTTPException(status_code=404, detail=f"Monitor for work_id {work_id} not found")

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting work monitor: {e}")
        raise HTTPException(status_code=500, detail=str(e))

