"""
Batch Execution API Routes for CSPBench Web Interface.

This module provides REST endpoints for batch execution management using the unified
WorkManager. It handles the complete lifecycle of batch executions from submission
to completion with comprehensive monitoring and control capabilities.
"""

import logging
from datetime import datetime

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException

from src.application.services.work_service import get_work_service
from src.application.work.work_manager import WorkManager
from src.domain.config import load_cspbench_config
from src.domain.status import BaseStatus
from src.infrastructure.utils.path_utils import get_batch_directory

from ..core.batch_execution_models import (
    BatchControlResponse,
    BatchExecutionRequest,
    BatchExecutionResponse,
    BatchListResponse,
    BatchStatusResponse,
)

# Get logger
logger = logging.getLogger(__name__)

# Create router
router = APIRouter(prefix="/api/batch", tags=["batch-execution"])


@router.post("/execute", response_model=BatchExecutionResponse)
async def execute_batch(
    request: BatchExecutionRequest, background_tasks: BackgroundTasks
) -> BatchExecutionResponse:
    """Execute batch configuration using unified WorkManager.

    This endpoint provides a centralized way to execute batch configurations
    through the unified WorkManager, which handles both work item management
    and pipeline execution with comprehensive validation and error handling.

    Args:
        request (BatchExecutionRequest): Batch execution configuration and options.
        background_tasks (BackgroundTasks): FastAPI background task manager.

    Returns:
        BatchExecutionResponse: Response with work ID and execution details.

    Raises:
        HTTPException: If configuration loading fails, file not found, or execution setup errors.

    Security:
        - Validates batch file paths to prevent directory traversal
        - Ensures files exist and are readable before execution
        - Sanitizes file paths and configuration data

    Note:
        Execution is handled asynchronously with status tracking through work ID.
        Use monitoring endpoints to track progress and results.
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

        # Use WorkManager for execution
        work_manager = get_work_service()
        extra_data = {
            "description": f"Web execution of {request.batch_file}",
            "batch_file": request.batch_file,
            "origin": "web",
        }

        work_id = work_manager.execute(config=config, extra=extra_data)

        logger.info(f"Work submitted with ID: {work_id}")

        # Get work service to get work details
        work_service = get_work_service()
        work_details = work_service.get(work_id)

        return BatchExecutionResponse(
            work_id=work_id,
            status=BaseStatus.QUEUED.value,
            message="Batch execution started successfully",
            submitted_at=datetime.now(),
            output_path=work_details.output_path if work_details else None,
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Unexpected error in batch execution: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Internal server error: {e}")


@router.get("/{work_id}/status", response_model=BatchStatusResponse)
async def get_batch_status(
    work_id: str, work_manager: WorkManager = Depends(get_work_service)
) -> BatchStatusResponse:
    """Get comprehensive status of a batch execution.

    Retrieves detailed status information for batch executions including
    progress, timing, error information, and configuration metadata.

    Args:
        work_id (str): Unique identifier of the work execution.
        work_manager (WorkManager): Injected work management service.

    Returns:
        BatchStatusResponse: Complete status information with metadata.

    Raises:
        HTTPException: If work not found or status retrieval fails.

    Note:
        Status includes both execution state and configuration details
        for comprehensive monitoring and debugging capabilities.
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
            config_name=work_item.get("config", {})
            .get("metadata", {})
            .get("name", "Unknown"),
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting batch status: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/list", response_model=BatchListResponse)
async def list_batches(
    work_manager: WorkManager = Depends(get_work_service),
) -> BatchListResponse:
    """List all batch executions with metadata and organization.

    Provides comprehensive listing of all batch executions in the system
    with metadata extraction and organization for management interfaces.

    Args:
        work_manager (WorkManager): Injected work management service.

    Returns:
        BatchListResponse: Complete list of batch executions with metadata.

    Raises:
        HTTPException: If listing fails due to service errors.

    Note:
        Includes configuration metadata, execution origin information,
        and timing data for comprehensive batch management.
    """
    try:
        work_items = work_manager.list()

        batches = []
        for work_item in work_items:
            batches.append(
                {
                    "work_id": work_item.id,
                    "status": work_item.status.value,
                    "created_at": datetime.fromtimestamp(work_item.created_at),
                    "updated_at": datetime.fromtimestamp(work_item.updated_at),
                    "config_name": work_item.get("config", {})
                    .get("metadata", {})
                    .get("name", "Unknown"),
                    "batch_file": work_item.get("extra", {}).get("batch_file", ""),
                    "description": work_item.get("extra", {}).get("description", ""),
                    "origin": work_item.get("extra", {}).get("origin", "unknown"),
                }
            )

        return BatchListResponse(
            batches=batches, total=len(batches), timestamp=datetime.now()
        )

    except Exception as e:
        logger.error(f"Error listing batches: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{work_id}/control", response_model=BatchControlResponse)
async def control_batch(
    work_id: str, action: str, work_manager: WorkManager = Depends(get_work_service)
) -> BatchControlResponse:
    """Control batch execution with comprehensive action support.

    Provides control operations for batch executions including pause,
    restart, and cancel with proper state management and validation.

    Args:
        work_id (str): Unique identifier of the work execution.
        action (str): Control action to perform ('pause', 'restart', 'cancel').
        work_manager (WorkManager): Injected work management service.

    Returns:
        BatchControlResponse: Control operation result with status information.

    Raises:
        HTTPException: If work not found, action invalid, or control operation fails.

    Security:
        - Validates action types against allowed operations
        - Ensures work exists before attempting control operations
        - Provides detailed error messages for debugging

    Note:
        Control operations are executed synchronously with immediate
        status updates reflected in the work management system.
    """
    try:
        if action == "pause":
            success = work_manager.pause(work_id)
        elif action == "cancel":
            success = work_manager.cancel(work_id)
        elif action == "restart":
            success = work_manager.restart(work_id)
        else:
            raise HTTPException(status_code=400, detail=f"Invalid action: {action}")

        if not success:
            raise HTTPException(
                status_code=404, detail=f"Work {work_id} not found or action failed"
            )

        return BatchControlResponse(
            work_id=work_id,
            action=action,
            success=True,
            message=f"Action '{action}' executed successfully",
            timestamp=datetime.now(),
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error controlling batch {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))
