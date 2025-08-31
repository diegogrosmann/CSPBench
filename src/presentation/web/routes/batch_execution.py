"""
Batch Execution API Routes
Provides REST endpoints for batch execution management using unified ExecutionManager.
"""

import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse

from src.application.services.execution_manager import ExecutionManager
from src.application.services.work_service import get_work_service
from src.domain.config import load_cspbench_config
from src.domain.errors import BatchConfigurationError
from src.infrastructure.utils.path_utils import get_batch_directory

import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException, BackgroundTasks, Depends
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from src.application.services.execution_manager import ExecutionManager
from src.application.work.work_manager import WorkManager
from src.domain import WorkStatus
from src.domain.status import BaseStatus
from src.application.services.work_service import get_work_service
from src.domain.config import load_cspbench_config

from ..core.batch_execution_models import (
    BatchExecutionRequest,
    BatchExecutionResponse,
    BatchStatusResponse,
    BatchResultsResponse,
    BatchListResponse,
    BatchControlResponse,
)

# Get logger
logger = logging.getLogger(__name__)

# Create router
router = APIRouter(prefix="/api/batch", tags=["batch-execution"])


@router.post("/execute", response_model=BatchExecutionResponse)
async def execute_batch(
    request: BatchExecutionRequest, background_tasks: BackgroundTasks
) -> BatchExecutionResponse:
    """
    Execute batch configuration using unified ExecutionManager.

    Args:
        request: Batch execution request with batch file path
        background_tasks: FastAPI background tasks

    Returns:
        Batch execution response with work_id and status
    """
    try:
        logger.info(f"Received batch execution request for file: {request.batch_file}")

        # Get batch directory and construct full path
        batch_dir = get_batch_directory()
        batch_path = batch_dir / request.batch_file

        # Validate file exists and is readable
        if not batch_path.exists():
            raise HTTPException(
                status_code=404, detail=f"Batch file not found: {request.batch_file}"
            )

        if not batch_path.is_file():
            raise HTTPException(
                status_code=400, detail=f"Path is not a file: {request.batch_file}"
            )

        # Load configuration using the same function as CLI
        try:
            config = load_cspbench_config(batch_path)
            logger.info(f"Loaded batch configuration: {config.metadata.name}")
        except Exception as e:
            logger.error(f"Config loading error: {e}")
            raise HTTPException(
                status_code=400, detail=f"Failed to load configuration: {e}"
            )

        # Use ExecutionManager without monitor
        execution_manager = ExecutionManager()
        extra_data = {
            "description": f"Web execution of {request.batch_file}",
            "batch_file": request.batch_file,
            "origin": "web",
        }

        work_id = execution_manager.execute(
            config=config, extra=extra_data
        )

        logger.info(f"Work submitted with ID: {work_id}")

        # Get work service to get work details
        work_service = get_work_service()
        work_details = work_service.get(work_id)

        return BatchExecutionResponse(
            work_id=work_id,
            status=BaseStatus.QUEUED.value,
            message="Batch execution started successfully",
            submitted_at=datetime.now(),
            output_path=work_details.get("output_path") if work_details else None,
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Unexpected error in batch execution: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Internal server error: {e}")


@router.get("/{work_id}/status", response_model=BatchStatusResponse)
async def get_batch_status(work_id: str, work_manager: WorkManager = Depends(get_work_service)) -> BatchStatusResponse:
    """
    Get status of a batch execution.

    Args:
        work_id: Work item ID

    Returns:
        Batch status response
    """
    try:
        work_item = work_manager.get(work_id)
        if not work_item:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")

        return BatchStatusResponse(
            work_id=work_id,
            status=work_item.get("status", "unknown"),
            created_at=datetime.fromtimestamp(work_item.get("created_at", 0)),
            updated_at=datetime.fromtimestamp(work_item.get("updated_at", 0)),
            output_path=work_item.get("output_path"),
            error_message=work_item.get("error"),
            config_name=work_item.get("config", {}).get("metadata", {}).get("name", "Unknown"),
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting batch status: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/list", response_model=BatchListResponse)
async def list_batches(work_manager: WorkManager = Depends(get_work_service)) -> BatchListResponse:
    """
    List all batch executions.

    Returns:
        List of batch executions
    """
    try:
        works = work_manager.list()
        
        batches = []
        for work in works:
            batches.append({
                "work_id": work.get("id", "unknown"),
                "status": work.get("status", "unknown"),
                "created_at": datetime.fromtimestamp(work.get("created_at", 0)),
                "updated_at": datetime.fromtimestamp(work.get("updated_at", 0)),
                "config_name": work.get("config", {}).get("metadata", {}).get("name", "Unknown"),
                "batch_file": work.get("extra", {}).get("batch_file", ""),
                "description": work.get("extra", {}).get("description", ""),
                "origin": work.get("extra", {}).get("origin", "unknown"),
            })

        return BatchListResponse(
            batches=batches,
            total=len(batches),
            timestamp=datetime.now()
        )

    except Exception as e:
        logger.error(f"Error listing batches: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{work_id}/control", response_model=BatchControlResponse)
async def control_batch(
    work_id: str, 
    action: str,
    work_manager: WorkManager = Depends(get_work_service)
) -> BatchControlResponse:
    """
    Control batch execution (pause, resume, cancel).

    Args:
        work_id: Work item ID
        action: Control action (pause, resume, cancel)

    Returns:
        Control response
    """
    try:
        if action == "pause":
            success = work_manager.pause(work_id)
        elif action == "resume":
            success = work_manager.resume(work_id)
        elif action == "cancel":
            success = work_manager.cancel(work_id)
        else:
            raise HTTPException(status_code=400, detail=f"Invalid action: {action}")

        if not success:
            raise HTTPException(
                status_code=404, 
                detail=f"Work {work_id} not found or action failed"
            )

        return BatchControlResponse(
            work_id=work_id,
            action=action,
            success=True,
            message=f"Action '{action}' executed successfully",
            timestamp=datetime.now()
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error controlling batch {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))
