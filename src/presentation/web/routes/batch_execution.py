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
from pydantic import BaseModel

from src.application.services.execution_manager import ExecutionManager
from src.domain import WorkStatus
from src.domain.status import BaseStatus
from src.presentation.display.web_monitor import WebMonitor
from src.application.services.work_service import (
    get_work_service,
    get_work_details,
    list_all_work,
    control_work,
)
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
        batch_dir = Path(os.getenv("BATCH_DIRECTORY", "./batches"))
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

        # Use ExecutionManager with async mode
        execution_manager = ExecutionManager()
        web_monitor = WebMonitor()  # TODO: Pass session info if needed
        extra_data = {
            "description": f"Web execution of {request.batch_file}",
            "batch_file": request.batch_file,
            "origin": "web",
        }

        work_id = execution_manager.execute(
            config=config, monitor=web_monitor, mode="async", extra=extra_data
        )

        logger.info(f"Work submitted with ID: {work_id}")

        # Get work details for response
        work_details = get_work_details(work_id)

        return BatchExecutionResponse(
            work_id=work_id,
            status=BaseStatus.QUEUED.value,
            message="Batch execution started successfully",
            submitted_at=datetime.now(),
            output_path=work_details.get("output_path"),
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Unexpected error in batch execution: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Internal server error: {e}")


import logging
import os
import threading
import yaml
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional

from fastapi import APIRouter, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from src.application.services.execution_manager import ExecutionManager
from src.presentation.display.web_monitor import WebMonitor
from src.application.work.global_manager import (
    get_global_work_manager,
    get_work_status,
    get_work_details,
    list_all_work,
    control_work,
)
from src.domain.config import CSPBenchConfig, load_cspbench_config

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
    Execute batch configuration with unified ExecutionManager.

    Args:
        request: Batch execution request with batch file path
        background_tasks: FastAPI background tasks

    Returns:
        Batch execution response with work_id and status
    """
    try:
        logger.info(f"Received batch execution request for file: {request.batch_file}")

        # Get batch directory and construct full path
        batch_dir = Path(os.getenv("BATCH_DIRECTORY", "./batches"))
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

        # Create WebMonitor for this execution
        web_monitor = WebMonitor()

        # Use ExecutionManager for unified orchestration
        execution_manager = ExecutionManager()
        extra_data = {
            "description": f"Web execution of {request.batch_file}",
            "batch_file": request.batch_file,
            "origin": "web",
        }

        # Execute in async mode (returns work_id immediately)
        work_id = execution_manager.execute(
            config=config, monitor=web_monitor, mode="async", extra=extra_data
        )

        logger.info(f"Work submitted with ID: {work_id}")

        # Get work details for response
        work_details = get_work_details(work_id)

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
async def get_batch_status(work_id: str) -> BatchStatusResponse:
    """
    Get batch execution status from global WorkManager.

    Args:
        work_id: Work item identifier

    Returns:
        Work status response with current state
    """
    try:
        logger.info(f"Getting status for work_id: {work_id}")

        # Get work details from global WorkManager
        work_details = get_work_details(work_id)
        if not work_details:
            raise HTTPException(
                status_code=404, detail=f"Work item '{work_id}' not found"
            )

        return BatchStatusResponse(
            work_id=work_id,
            status=work_details["status"],
            created_at=work_details["created_at"],
            updated_at=work_details["updated_at"],
            output_path=work_details.get("output_path"),
            progress=work_details.get("extra", {}).get("progress"),
            error=work_details.get("error"),
            config_name=work_details.get("extra", {}).get("config_name"),
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting status for work_id {work_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Failed to get work status: {e}")


@router.get("/{work_id}/results", response_model=BatchResultsResponse)
async def get_batch_results(work_id: str) -> BatchResultsResponse:
    """
    Get batch execution results from global WorkManager.

    Args:
        work_id: Work item identifier

    Returns:
        Work results response with execution output
    """
    try:
        logger.info(f"Getting results for work_id: {work_id}")

        # Get work details from global WorkManager
        work_details = get_work_details(work_id)
        if not work_details:
            raise HTTPException(
                status_code=404, detail=f"Work item '{work_id}' not found"
            )

        # Check if work is completed
        if work_details["status"] not in [
            BaseStatus.COMPLETED.value,
            BaseStatus.FAILED.value,
        ]:
            raise HTTPException(
                status_code=409,
                detail=f"Work item is not completed (status: {work_details['status']})",
            )

        return BatchResultsResponse(
            work_id=work_id,
            status=work_details["status"],
            output_path=work_details.get("output_path"),
            error=work_details.get("error"),
            output_files=[],  # TODO: Implement file listing
            summary={},  # TODO: Implement execution summary
            details=work_details,
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting results for work_id {work_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Failed to get work results: {e}")


@router.get("/list", response_model=BatchListResponse)
async def list_batch_work(status: Optional[str] = None) -> BatchListResponse:
    """
    List all batch work items from global WorkManager.

    Args:
        status: Optional status filter

    Returns:
        List of all work items with their status
    """
    try:
        logger.info(f"Listing batch work items with status filter: {status}")

        # Get all work items from global WorkManager
        work_items = list_all_work()

        # Transform to expected format
        executions = []
        for item in work_items:
            # Convert Unix timestamps to datetime
            from datetime import datetime

            # Handle both old and new format
            work_id = item.get("work_id") or item.get("id")
            created_at_value = item.get("created_at")
            updated_at_value = item.get("updated_at")

            # Convert timestamps safely
            if isinstance(created_at_value, str):
                created_at = datetime.fromisoformat(
                    created_at_value.replace("Z", "+00:00")
                )
            else:
                created_at = datetime.fromtimestamp(created_at_value)

            if isinstance(updated_at_value, str):
                updated_at = datetime.fromisoformat(
                    updated_at_value.replace("Z", "+00:00")
                )
            else:
                updated_at = datetime.fromtimestamp(updated_at_value)

            execution = BatchStatusResponse(
                work_id=work_id,
                status=item["status"],
                created_at=created_at,
                updated_at=updated_at,
                output_path=item.get("output_path"),
                progress=item.get("progress"),
                error=item.get("error"),
                config_name=item.get("config_name"),
            )
            executions.append(execution)

        # Apply status filter if provided
        total_count = len(executions)
        if status:
            executions = [exec for exec in executions if exec.status == status]

        filtered_count = len(executions)

        return BatchListResponse(
            total=total_count, executions=executions, filtered=filtered_count
        )

    except Exception as e:
        logger.error(f"Error listing batch work: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Failed to list batch work: {e}")


# Simple request model for control actions
class BatchControlRequest(BaseModel):
    """Request model for batch control actions."""

    action: str


@router.post("/{work_id}/control", response_model=BatchControlResponse)
async def control_batch_execution(
    work_id: str, request: BatchControlRequest
) -> BatchControlResponse:
    """
    Control batch execution (pause, resume, cancel, restart).

    Args:
        work_id: Work item identifier
        request: Control request with action

    Returns:
        Control response with operation result
    """
    try:
        logger.info(f"Control action '{request.action}' for work_id: {work_id}")

        # Validate action
        valid_actions = ["pause", "resume", "cancel", "restart"]
        if request.action not in valid_actions:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid action '{request.action}'. Valid actions: {valid_actions}",
            )

        # Execute control action using global WorkManager
        success = control_work(work_id, request.action)

        if not success:
            # Check if work item exists
            work_details = get_work_details(work_id)
            if not work_details:
                raise HTTPException(
                    status_code=404, detail=f"Work item '{work_id}' not found"
                )

            # Work exists but action failed
            raise HTTPException(
                status_code=409,
                detail=f"Cannot {request.action} work item in current state: {work_details['status']}",
            )

        return BatchControlResponse(
            work_id=work_id,
            operation=request.action,
            success=True,
            message=f"Successfully {request.action}{'ed' if request.action.endswith('e') else 'ed'} work item",
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error controlling work_id {work_id}: {e}", exc_info=True)
        raise HTTPException(
            status_code=500, detail=f"Failed to control work execution: {e}"
        )


@router.delete("/{work_id}")
async def delete_batch_work(work_id: str) -> JSONResponse:
    """
    Delete batch work item from global WorkManager.

    Args:
        work_id: Work item identifier

    Returns:
        JSON response confirming deletion
    """
    try:
        logger.info(f"Deleting work_id: {work_id}")

        # Check if work item exists
        work_details = get_work_details(work_id)
        if not work_details:
            raise HTTPException(
                status_code=404, detail=f"Work item '{work_id}' not found"
            )

        # TODO: Implement actual deletion in WorkManager/Repository
        # For now, we'll cancel the work if it's running
        if work_details["status"] in ["queued", "running"]:
            control_work(work_id, "cancel")

        return JSONResponse(
            status_code=200,
            content={
                "message": f"Work item '{work_id}' deleted successfully",
                "work_id": work_id,
            },
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error deleting work_id {work_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Failed to delete work item: {e}")
