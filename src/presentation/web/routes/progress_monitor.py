"""
Real-time monitoring API endpoints for algorithm progress tracking.
Updated to use the new work management system.
"""

from datetime import datetime
from fastapi import APIRouter, HTTPException, Depends
from pathlib import Path
from typing import Dict, Any, List, Optional

from src.application.services.work_service import get_work_service
from src.application.work.work_manager import WorkManager

try:
    from src.infrastructure.persistence.work_state.queries import WorkStateQueries
except ImportError:
    WorkStateQueries = None

from src.infrastructure.logging_config import get_logger

logger = get_logger(__name__)
router = APIRouter(prefix="/api/monitoring", tags=["monitoring"])


def format_progress_bar(progress: float, width: int = 20) -> str:
    """Create ASCII progress bar."""
    filled = int(progress * width)
    empty = width - filled
    return f"[{'█' * filled}{'-' * empty}]"


def format_timestamp(timestamp: Optional[float]) -> str:
    """Format timestamp for display."""
    if not timestamp:
        return "N/A"
    return datetime.fromtimestamp(timestamp).strftime("%H:%M:%S")


def format_duration(start_time: Optional[float], end_time: Optional[float] = None) -> str:
    """Format duration for display."""
    if not start_time:
        return "N/A"
    
    end = end_time or datetime.now().timestamp()
    duration = end - start_time
    
    if duration < 60:
        return f"{int(duration)}s"
    elif duration < 3600:
        return f"{int(duration // 60)}m {int(duration % 60)}s"
    else:
        hours = int(duration // 3600)
        minutes = int((duration % 3600) // 60)
        return f"{hours}h {minutes}m"


@router.get("/progress/{work_id}")
async def get_work_progress(work_id: str, work_manager: WorkManager = Depends(get_work_service)) -> Dict[str, Any]:
    """Get comprehensive progress information for a work item."""
    print(f"DEBUG: get_work_progress called for work_id: {work_id}")
    try:
        # Get work item from work manager
        work_item = work_manager.get(work_id)
        if not work_item:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        # Check if work state database exists
        work_output_path = work_item.get("output_path", "")
        state_db_path = Path(work_output_path) / "state.db" if work_output_path else None
        
        print(f"DEBUG (progress_monitor) Work output path: {work_output_path}")
        print(f"DEBUG (progress_monitor) State DB path: {state_db_path}")
        print(f"DEBUG (progress_monitor) State DB exists: {state_db_path.exists() if state_db_path else False}")
        print(f"DEBUG (progress_monitor) WorkStateQueries available: {WorkStateQueries is not None}")
        
        if not state_db_path or not state_db_path.exists() or not WorkStateQueries:
            # Return basic work info if no state database or queries not available
            return {
                "work_id": work_id,
                "work_status": work_item.get("status", "unknown"),
                "progress_summary": {
                    "task": {"current": 0, "total": 0},
                    "dataset": {"current": 0, "total": 0},
                    "config": {"current": 0, "total": 0},
                    "algorithm": {"current": 0, "total": 0},
                    "sequence": {"current": 0, "total": 0}
                },
                "running_executions": [],
                "errors": [],
                "warnings": [],
                "terminal_output": ""
            }
        
        with WorkStateQueries(state_db_path) as queries:
            # Check if work exists in state database
            if not queries.work_exists(work_id):
                raise HTTPException(status_code=404, detail=f"Work {work_id} not found in state database")
            
            # Get progress summary
            progress_summary = queries.get_work_progress_summary(work_id)
            if not progress_summary:
                raise HTTPException(status_code=404, detail=f"No progress data for work {work_id}")
            
            # Get running executions detail
            running_executions = queries.get_running_executions_detail(work_id, limit=20)
            
            # Get error summary
            errors = queries.get_error_summary(work_id, limit=10)
            
            # Get warnings
            warnings = queries.get_execution_warnings(work_id, limit=10)
            
            # Get work info
            work_info = queries.get_work_info(work_id)
            
            return {
                "work_id": work_id,
                "work_status": work_info.get("status") if work_info else work_item.get("status", "unknown"),
                "progress_summary": {
                    "task": {
                        "current": progress_summary.tasks["Finished"],
                        "total": progress_summary.tasks["Total"]
                    },
                    "dataset": {
                        "current": progress_summary.datasets["Finished"],
                        "total": progress_summary.datasets["Total"]
                    },
                    "config": {
                        "current": progress_summary.configs["Finished"],
                        "total": progress_summary.configs["Total"]
                    },
                    "algorithm": {
                        "current": progress_summary.algorithms["Finished"],
                        "total": progress_summary.algorithms["Total"]
                    },
                    "sequence": {
                        "current": progress_summary.execution["Finished"],
                        "total": progress_summary.execution["Total"]
                    }
                },
                "global_progress": progress_summary.global_progress,
                "current_combination": progress_summary.current_combination_details,
                "running_executions": [
                    {
                        "unit_id": exec.unit_id,
                        "combination_id": exec.combination_id,
                        "sequencia": exec.sequencia,
                        "task_id": exec.task_id,
                        "dataset_id": exec.dataset_id,
                        "preset_id": exec.preset_id,
                        "algorithm_id": exec.algorithm_id,
                        "progress": exec.progress,
                        "progress_message": exec.progress_message,
                        "start_time": exec.started_at,
                        "duration": format_duration(exec.started_at) if exec.started_at else "Unknown",
                        "formatted_time": format_timestamp(exec.started_at) if exec.started_at else "Unknown"
                    }
                    for exec in running_executions
                ],
                "errors": [
                    {
                        "unit_id": error.unit_id,
                        "error_type": error.error_type,
                        "error_message": error.error_message,
                        "timestamp": error.timestamp,
                        "formatted_time": format_timestamp(error.timestamp)
                    }
                    for error in errors
                ],
                "warnings": [
                    {
                        "message": warning.get("message", "Unknown warning"),
                        "unit_id": warning.get("unit_id"),
                        "timestamp": warning.get("timestamp"),
                        "formatted_time": format_timestamp(warning.get("timestamp")) if warning.get("timestamp") else "Unknown"
                    }
                    for warning in warnings
                ],
                "terminal_output": ""  # Could be enhanced to read from log files
            }

    except Exception as e:
        logger.error(f"Error getting work progress: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/terminal/{work_id}")
async def get_terminal_output(work_id: str, work_manager: WorkManager = Depends(get_work_service)):
    """Get terminal output for a work item."""
    try:
        # Get work item from work manager
        work_item = work_manager.get(work_id)
        if not work_item:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        # Try to read log files from the work output directory
        work_output_path = work_item.get("output_path", "")
        terminal_output = []
        
        if work_output_path:
            output_dir = Path(work_output_path)
            log_files = list(output_dir.glob("*.log"))
            
            # Read recent log content
            for log_file in log_files:
                try:
                    with open(log_file, "r") as f:
                        lines = f.readlines()
                        # Get last 100 lines
                        recent_lines = lines[-100:] if len(lines) > 100 else lines
                        terminal_output.extend([line.rstrip() for line in recent_lines])
                except Exception as e:
                    logger.warning(f"Error reading log file {log_file}: {e}")
        
        return {
            "work_id": work_id,
            "terminal_output": "\n".join(terminal_output) if terminal_output else "No log output available",
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        logger.error(f"Error getting terminal output: {e}")
        raise HTTPException(status_code=500, detail=str(e))
