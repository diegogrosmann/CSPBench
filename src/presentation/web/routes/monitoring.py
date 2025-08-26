"""
API routes for real-time monitoring and execution management
Updated to work with the refactored system without monitor components.
"""

from fastapi import APIRouter, HTTPException, Depends
from typing import List, Dict, Any, Optional
from datetime import datetime, timedelta
import os
import glob
from pathlib import Path

from src.application.services.work_service import get_work_service
from src.application.work.work_manager import WorkManager
from src.domain.status import BaseStatus
from src.presentation.web.websocket_manager import connection_manager
from src.infrastructure.logging_config import get_logger

# Get logger
logger = get_logger(__name__)

try:
    from src.infrastructure.persistence.work_state.queries import WorkStateQueries
except ImportError:
    # Fallback if queries not available
    WorkStateQueries = None

logger = get_logger(__name__)
router = APIRouter()


@router.get("/api/monitoring/list")
async def list_works(work_manager: WorkManager = Depends(get_work_service)):
    """Get list of all work items from database"""
    try:
        # Get all work items from the work manager
        works = work_manager.list()
        
        # Format for frontend
        formatted_works = []
        for work in works:
            formatted_works.append({
                "work_id": work.get("id", "unknown"),
                "status": work.get("status", "unknown"),
                "created_at": work.get("created_at", 0),
                "updated_at": work.get("updated_at", 0),
                "output_path": work.get("output_path", ""),
                "config_name": work.get("config", {}).get("metadata", {}).get("name", "Unknown"),
                "description": work.get("extra", {}).get("description", ""),
                "batch_file": work.get("extra", {}).get("batch_file", ""),
                "origin": work.get("extra", {}).get("origin", "unknown")
            })
        
        # Sort by creation time (newest first)
        formatted_works.sort(key=lambda x: x["created_at"], reverse=True)
        
        return formatted_works

    except Exception as e:
        logger.error(f"Error listing works: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/work/{work_id}/status")
async def get_work_status(work_id: str, work_manager: WorkManager = Depends(get_work_service)):
    """Get work status and basic information"""
    try:
        work_item = work_manager.get(work_id)
        if not work_item:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        return {
            "work_id": work_id,
            "status": work_item.get("status", "unknown"),
            "created_at": work_item.get("created_at", 0),
            "updated_at": work_item.get("updated_at", 0),
            "error": work_item.get("error"),
            "config_name": work_item.get("config", {}).get("metadata", {}).get("name", "Unknown"),
            "description": work_item.get("extra", {}).get("description", ""),
        }

    except Exception as e:
        logger.error(f"Error getting work status: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# @router.get("/api/monitoring/work/{work_id}/progress")
# async def get_work_progress(work_id: str, work_manager: WorkManager = Depends(get_work_service)):
#     """Get detailed progress information for a work item"""
#     try:
#         # Get work item from database
#         work_item = work_manager.get(work_id)
#         if not work_item:
#             raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
#         
#         # Check if work state database exists for detailed progress
#         work_output_path = work_item.get("output_path", "")
#         state_db_path = Path(work_output_path) / "state.db" if work_output_path else None
#         
#         print(f"DEBUG State DB path: {state_db_path}")
#         print(f"DEBUG State DB exists: {state_db_path.exists() if state_db_path else False}")
#         print(f"DEBUG WorkStateQueries available: {WorkStateQueries is not None}")
#         
#         # Default progress data structure
#         progress_data = {
#             "work_id": work_id,
#             "status": work_item.get("status", "unknown"),
#             "updated_at": work_item.get("updated_at", 0),
#             "error": work_item.get("error"),
#             "progress_summary": {
#                 "task": {"current": 0, "total": 1},
#                 "dataset": {"current": 0, "total": 1}, 
#                 "config": {"current": 0, "total": 1},
#                 "algorithm": {"current": 0, "total": 1},
#                 "sequence": {"current": 0, "total": 1}
#             },
#             "running_executions": [],
#             "errors": [],
#             "warnings": [],
#             "console_output": []
#         }
#         
#         # Simulate progress based on status for demo purposes
#         status = work_item.get("status", "unknown")
#         if status == "completed":
#             progress_data["progress_summary"] = {
#                 "task": {"current": 1, "total": 1},
#                 "dataset": {"current": 1, "total": 1}, 
#                 "config": {"current": 1, "total": 1},
#                 "algorithm": {"current": 1, "total": 1},
#                 "sequence": {"current": 1, "total": 1}
#             }
#         elif status == "running":
#             progress_data["progress_summary"] = {
#                 "task": {"current": 1, "total": 2},
#                 "dataset": {"current": 1, "total": 1}, 
#                 "config": {"current": 1, "total": 1},
#                 "algorithm": {"current": 0, "total": 1},
#                 "sequence": {"current": 50, "total": 100}
#             }
#         
#         # If work state database exists and WorkStateQueries available, get detailed progress
#         print(f"DEBUG Entering query section...")
#         print(f"DEBUG WorkStateQueries check: {WorkStateQueries}")
#         print(f"DEBUG state_db_path check: {state_db_path}")
#         print(f"DEBUG state_db_path exists: {state_db_path.exists() if state_db_path else 'Path is None'}")
#         
#         if WorkStateQueries and state_db_path and state_db_path.exists():
#             print(f"DEBUG Entering WorkStateQueries section...")
#             try:
#                 with WorkStateQueries(state_db_path) as queries:
#                     if queries.work_exists(work_id):
#                         # Get progress summary
#                         summary = queries.get_work_progress_summary(work_id)
#                         if summary:
#                             progress_data["progress_summary"] = {
#                                 "task": {"current": summary.tasks["Finished"], "total": summary.tasks["Total"]},
#                                 "dataset": {"current": summary.datasets["Finished"], "total": summary.datasets["Total"]},
#                                 "config": {"current": summary.configs["Finished"], "total": summary.configs["Total"]},
#                                 "algorithm": {"current": summary.algorithms["Finished"], "total": summary.algorithms["Total"]},
#                                 "sequence": {"current": summary.execution["Finished"], "total": summary.execution["Total"]}
#                             }
#                             progress_data["global_progress"] = summary.global_progress
#                             progress_data["current_combination"] = summary.current_combination_details
#                         
#                         # Get running executions
#                         running = queries.get_running_executions_detail(work_id, limit=20)
#                         progress_data["running_executions"] = [
#                             {
#                                 "unit_id": exec.unit_id,
#                                 "combination_id": exec.combination_id,
#                                 "sequencia": exec.sequencia,
#                                 "task_id": exec.task_id, 
#                                 "dataset_id": exec.dataset_id,
#                                 "preset_id": exec.preset_id,
#                                 "algorithm_id": exec.algorithm_id,
#                                 "progress": exec.progress,
#                                 "progress_message": exec.progress_message,
#                                 "started_at": exec.started_at,
#                                 "formatted_time": datetime.fromtimestamp(exec.started_at).strftime("%H:%M:%S") if exec.started_at else "N/A"
#                             }
#                             for exec in running
#                         ]
#                         
#                         # Get errors
#                         errors = queries.get_error_summary(work_id, limit=10)
#                         progress_data["errors"] = [
#                             {
#                                 "unit_id": error.unit_id,
#                                 "error_type": error.error_type,
#                                 "error_message": error.error_message,
#                                 "timestamp": error.timestamp,
#                                 "formatted_time": datetime.fromtimestamp(error.timestamp).strftime("%H:%M:%S") if error.timestamp else "N/A"
#                             }
#                             for error in errors
#                         ]
#                         
#                         # Get warnings  
#                         warnings = queries.get_execution_warnings(work_id, limit=10)
#                         progress_data["warnings"] = [
#                             {
#                                 "message": warning.get("message", "Unknown warning"),
#                                 "unit_id": warning.get("unit_id"),
#                                 "timestamp": warning.get("timestamp"),
#                                 "formatted_time": datetime.fromtimestamp(warning.get("timestamp")).strftime("%H:%M:%S") if warning.get("timestamp") else "N/A"
#                             }
#                             for warning in warnings
#                         ]
#                         
#             except Exception as e:
#                 logger.warning(f"Error reading state database for {work_id}: {e}")
#         
#         return progress_data
# 
#     except Exception as e:
#         logger.error(f"Error getting work progress: {e}")
#         raise HTTPException(status_code=500, detail=str(e))


from pydantic import BaseModel

# Add control request model
class WorkControlRequest(BaseModel):
    action: str


@router.post("/api/monitoring/work/{work_id}/control")
async def control_work(
    work_id: str, 
    request: WorkControlRequest,
    work_manager: WorkManager = Depends(get_work_service)
):
    """Control work execution (pause, resume, cancel)"""
    try:
        action = request.action
        if action == "pause":
            success = work_manager.pause(work_id)
        elif action == "resume":
            success = work_manager.resume(work_id)
        elif action == "cancel":
            success = work_manager.cancel(work_id)
        else:
            raise HTTPException(status_code=400, detail=f"Invalid action: {action}")
        
        if not success:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found or action failed")
        
        return {"success": True, "action": action, "work_id": work_id}

    except Exception as e:
        logger.error(f"Error controlling work {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/monitoring/terminal/{work_id}")
async def get_terminal_output(work_id: str, work_manager: WorkManager = Depends(get_work_service)):
    """Get terminal output for a work item."""
    try:
        from pathlib import Path
        
        # Get work item from work manager
        work_item = work_manager.get(work_id)
        if not work_item:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        # Try to read log files from the work output directory
        work_output_path = work_item.get("output_path", "")
        terminal_output = []
        
        if work_output_path:
            output_dir = Path(work_output_path)
            if output_dir.exists():
                log_files = list(output_dir.glob("*.log"))
                
                # Read recent log content
                for log_file in log_files:
                    try:
                        with open(log_file, "r") as f:
                            lines = f.readlines()
                            # Get last 100 lines
                            recent_lines = lines[-100:] if len(lines) > 100 else lines
                            terminal_output.extend([line.rstrip() for line in recent_lines])
                    except Exception as file_error:
                        logger.warning(f"Could not read log file {log_file}: {file_error}")
        
        return {
            "work_id": work_id,
            "output": "\n".join(terminal_output) if terminal_output else "No log output available",
            "timestamp": datetime.now()
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting terminal output: {e}")
        raise HTTPException(status_code=500, detail=str(e))
