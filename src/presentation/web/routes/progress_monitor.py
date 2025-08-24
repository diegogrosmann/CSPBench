"""
Real-time monitoring API endpoints for algorithm progress tracking.
"""

from datetime import datetime
from fastapi import APIRouter, HTTPException
from pathlib import Path
from typing import Dict, Any, List, Optional

from src.infrastructure.persistence.work_state.queries import WorkStateQueries
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
async def get_work_progress(work_id: str) -> Dict[str, Any]:
    """Get comprehensive progress information for a work item."""
    try:
        # Check if work state database exists
        state_db_path = Path("data/outputs") / work_id / "state.db"
        if not state_db_path.exists():
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        with WorkStateQueries(state_db_path) as queries:
            # Check if work exists
            if not queries.work_exists(work_id):
                raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
            
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
                "work_status": work_info.get("status") if work_info else "unknown",
                "progress_summary": {
                    "task": {
                        "current": progress_summary.current_task,
                        "total": progress_summary.total_tasks
                    },
                    "dataset": {
                        "current": progress_summary.current_dataset,
                        "total": progress_summary.total_datasets
                    },
                    "config": {
                        "current": progress_summary.current_config,
                        "total": progress_summary.total_configs
                    },
                    "algorithm": {
                        "current": progress_summary.current_algorithm,
                        "total": progress_summary.total_algorithms
                    },
                    "sequence": {
                        "current": progress_summary.current_sequence,
                        "total": progress_summary.total_sequences
                    },
                    "executions": {
                        "completed": progress_summary.completed_executions,
                        "total": progress_summary.total_executions
                    },
                    "overall_progress": progress_summary.overall_progress
                },
                "running_executions": [
                    {
                        "unit_id": exec_detail.unit_id,
                        "sequencia": exec_detail.sequencia,
                        "progress": exec_detail.progress,
                        "progress_message": exec_detail.progress_message,
                        "task_id": exec_detail.task_id,
                        "dataset_id": exec_detail.dataset_id,
                        "algorithm_id": exec_detail.algorithm_id,
                        "started_at": exec_detail.started_at,
                        "duration": format_duration(exec_detail.started_at),
                        "total_sequences": exec_detail.total_sequences
                    }
                    for exec_detail in running_executions
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
                        "message": warning["message"],
                        "unit_id": warning.get("unit_id"),
                        "timestamp": warning["timestamp"],
                        "formatted_time": format_timestamp(warning["timestamp"])
                    }
                    for warning in warnings
                ]
            }
    
    except Exception as e:
        logger.error(f"Error getting progress for work {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/progress/{work_id}/console")
async def get_console_view(work_id: str) -> Dict[str, Any]:
    """Get console-formatted progress view for terminal display."""
    try:
        # Get the progress data
        progress_data = await get_work_progress(work_id)
        
        # Build console output
        summary = progress_data["progress_summary"]
        running = progress_data["running_executions"]
        errors = progress_data["errors"]
        warnings = progress_data["warnings"]
        
        # Header line
        header = (
            f"Task: {summary['task']['current']}/{summary['task']['total']} | "
            f"DataSet: {summary['dataset']['current']}/{summary['dataset']['total']} | "
            f"Config: {summary['config']['current']}/{summary['config']['total']} | "
            f"Alg: {summary['algorithm']['current']}/{summary['algorithm']['total']} | "
            f"Run Seq ({summary['sequence']['current']}/{summary['sequence']['total']})"
        )
        
        # Table header
        table_lines = [
            "+-----------------------------------------------------------+",
            "| SEQ   Unit ID      Progress                                |",
            "|-----------------------------------------------------------|"
        ]
        
        # Running executions
        for exec_data in running[:10]:  # Limit to 10 for display
            unit_id = exec_data["unit_id"][:10]  # Truncate unit ID
            progress = exec_data["progress"]
            sequencia = exec_data["sequencia"]
            total_seq = exec_data["total_sequences"]
            
            # Progress bar
            progress_bar = format_progress_bar(progress, 15)
            progress_text = f"{int(progress * total_seq)}/{total_seq}  {int(progress * 100)}%"
            
            # Main line
            main_line = f"| {sequencia:<4} {unit_id:<10} {progress_bar} {progress_text:<15} |"
            table_lines.append(main_line)
            
            # Progress message
            if exec_data["progress_message"]:
                msg = exec_data["progress_message"][:50]  # Truncate message
                msg_line = f"|        Progress Message: {msg:<35} |"
                table_lines.append(msg_line)
            
            # Check for warnings for this unit
            unit_warnings = [w for w in warnings if w.get("unit_id") == exec_data["unit_id"]]
            if unit_warnings:
                warn_msg = unit_warnings[0]["message"][:40]  # Truncate warning
                warn_line = f"|        Warning: {warn_msg:<43} |"
                table_lines.append(warn_line)
            
            table_lines.append("|                                                            |")
        
        table_lines.append("+-----------------------------------------------------------+")
        
        # Errors section
        error_lines = []
        if errors:
            error_lines.append("")
            error_lines.append("Errors:")
            for error in errors[:5]:  # Limit to 5 errors
                error_line = f"Unit ID: {error['unit_id']} - {error['error_message'][:50]}"
                error_lines.append(error_line)
        
        # Combine all output
        console_output = [header, ""] + table_lines + error_lines
        
        return {
            "work_id": work_id,
            "console_output": "\n".join(console_output),
            "raw_data": progress_data
        }
    
    except Exception as e:
        logger.error(f"Error getting console view for work {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/status/{work_id}")
async def get_work_status(work_id: str) -> Dict[str, Any]:
    """Get basic status information for a work item."""
    try:
        state_db_path = Path("data/outputs") / work_id / "state.db"
        if not state_db_path.exists():
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        with WorkStateQueries(state_db_path) as queries:
            work_info = queries.get_work_info(work_id)
            if not work_info:
                raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
            
            combination_counts = queries.get_combination_status_counts(work_id)
            execution_counts = queries.get_execution_status_counts(work_id)
            
            return {
                "work_id": work_id,
                "status": work_info["status"],
                "created_at": work_info["created_at"],
                "updated_at": work_info["updated_at"],
                "output_path": work_info["output_path"],
                "error": work_info["error"],
                "combination_counts": combination_counts,
                "execution_counts": execution_counts
            }
    
    except Exception as e:
        logger.error(f"Error getting status for work {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/list")
async def list_active_works() -> List[Dict[str, Any]]:
    """List all active work items with basic status."""
    try:
        outputs_dir = Path("data/outputs")
        if not outputs_dir.exists():
            return []
        
        active_works = []
        
        for work_dir in outputs_dir.iterdir():
            if work_dir.is_dir():
                work_id = work_dir.name
                state_db_path = work_dir / "state.db"
                
                if state_db_path.exists():
                    try:
                        with WorkStateQueries(state_db_path) as queries:
                            work_info = queries.get_work_info(work_id)
                            if work_info:
                                active_works.append({
                                    "work_id": work_id,
                                    "status": work_info["status"],
                                    "created_at": work_info["created_at"],
                                    "updated_at": work_info["updated_at"]
                                })
                    except Exception as e:
                        logger.warning(f"Error reading work {work_id}: {e}")
                        continue
        
        # Sort by update time, most recent first
        active_works.sort(key=lambda x: x["updated_at"] or 0, reverse=True)
        
        return active_works
    
    except Exception as e:
        logger.error(f"Error listing active works: {e}")
        raise HTTPException(status_code=500, detail=str(e))
