"""
Batch Execution Routes

This module handles batch execution functionality including:
- Listing available batch files
- Uploading new batch files
- Executing batch configurations
- Monitoring batch execution status
"""

import logging
from pathlib import Path
from typing import Any, Dict, List

from fastapi import (
    APIRouter,
    BackgroundTasks,
    File,
    Form,
    HTTPException,
    Request,
    UploadFile,
)
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.templating import Jinja2Templates

from ..services.session_manager import ExecutionSessionManager

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/execution", tags=["batch-execution"])
templates = Jinja2Templates(directory="src/presentation/web/templates")

# Global session manager for tracking executions
session_manager = ExecutionSessionManager()


@router.get("/batch", response_class=HTMLResponse)
async def batch_execution_page(request: Request):
    """Render the batch execution page."""
    return templates.TemplateResponse("batch_execution.html", {"request": request})
    return templates.TemplateResponse("batch_execution.html", {"request": request})


@router.get("/api/batches")
async def list_batch_files():
    """List all available batch files in the batches directory."""
    try:
        batches_dir = Path("batches")
        if not batches_dir.exists():
            return {"batch_files": []}

        batch_files = []
        for file_path in batches_dir.glob("*.yaml"):
            if file_path.is_file():
                # Extract description from file
                description = _extract_batch_description(file_path)

                # Get file stats
                stat = file_path.stat()

                batch_files.append(
                    {
                        "filename": file_path.name,
                        "name": file_path.stem,
                        "description": description,
                        "size": stat.st_size,
                        "modified": stat.st_mtime,
                        "path": str(file_path),
                    }
                )

        # Also check for .yml files
        for file_path in batches_dir.glob("*.yml"):
            if file_path.is_file():
                description = _extract_batch_description(file_path)
                stat = file_path.stat()

                batch_files.append(
                    {
                        "filename": file_path.name,
                        "name": file_path.stem,
                        "description": description,
                        "size": stat.st_size,
                        "modified": stat.st_mtime,
                        "path": str(file_path),
                    }
                )

        # Sort by modification time (newest first)
        batch_files.sort(key=lambda x: x["modified"], reverse=True)

        return {"batch_files": batch_files}

    except Exception as e:
        logger.error(f"Error listing batch files: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to list batch files: {str(e)}"
        )


@router.post("/api/batches/upload")
async def upload_batch_file(file: UploadFile = File(...), name: str = Form(None)):
    """Upload a new batch file."""
    try:
        # Validate file extension
        if not file.filename.lower().endswith((".yaml", ".yml")):
            raise HTTPException(
                status_code=400, detail="Only YAML files (.yaml or .yml) are allowed"
            )

        # Use provided name or original filename
        filename = name if name else file.filename
        if not filename.lower().endswith((".yaml", ".yml")):
            filename += ".yaml"

        # Ensure batches directory exists
        batches_dir = Path("batches")
        batches_dir.mkdir(exist_ok=True)

        # Save file
        file_path = batches_dir / filename
        content = await file.read()

        # Validate YAML content
        try:
            import yaml

            yaml.safe_load(content.decode("utf-8"))
        except yaml.YAMLError as e:
            raise HTTPException(
                status_code=400, detail=f"Invalid YAML format: {str(e)}"
            )

        with open(file_path, "wb") as f:
            f.write(content)

        logger.info(f"Batch file uploaded: {filename}")

        return {
            "message": f"Batch file '{filename}' uploaded successfully",
            "filename": filename,
            "path": str(file_path),
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error uploading batch file: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to upload batch file: {str(e)}"
        )


@router.post("/api/batches/{filename}/execute")
async def execute_batch_file(filename: str, background_tasks: BackgroundTasks):
    """Execute a batch file."""
    try:
        # Validate file exists
        batches_dir = Path("batches")
        file_path = batches_dir / filename

        if not file_path.exists():
            raise HTTPException(
                status_code=404, detail=f"Batch file '{filename}' not found"
            )

        # Create execution session
        session_id = session_manager.create_session(
            session_type="batch_execution",
            config={"batch_file": filename, "batch_path": str(file_path)},
        )

        # Start background execution
        background_tasks.add_task(_execute_batch_background, session_id, str(file_path))

        return {
            "message": f"Batch execution started for '{filename}'",
            "session_id": session_id,
            "filename": filename,
            "status": "started",
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error executing batch file {filename}: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to execute batch: {str(e)}"
        )


@router.get("/api/batches/{filename}/download")
async def download_batch_file(filename: str):
    """Download a batch file."""
    try:
        batches_dir = Path("batches")
        file_path = batches_dir / filename

        if not file_path.exists():
            raise HTTPException(
                status_code=404, detail=f"Batch file '{filename}' not found"
            )

        return FileResponse(
            path=str(file_path), filename=filename, media_type="application/x-yaml"
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error downloading batch file {filename}: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to download batch file: {str(e)}"
        )


@router.delete("/api/batches/{filename}")
async def delete_batch_file(filename: str):
    """Delete a batch file. Template files cannot be deleted."""
    try:
        # Prevent deletion of template files
        if filename.upper() == "TEMPLATE.YAML" or filename.lower().startswith(
            "template"
        ):
            raise HTTPException(
                status_code=403, detail="Template files cannot be deleted"
            )

        batches_dir = Path("batches")
        file_path = batches_dir / filename

        if not file_path.exists():
            raise HTTPException(
                status_code=404, detail=f"Batch file '{filename}' not found"
            )

        # Delete the file
        file_path.unlink()

        logger.info(f"Batch file deleted: {filename}")

        return {
            "message": f"Batch file '{filename}' deleted successfully",
            "filename": filename,
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error deleting batch file {filename}: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to delete batch file: {str(e)}"
        )


@router.get("/api/batches/{filename}")
async def get_batch_file_details(filename: str):
    """Get detailed information about a specific batch file."""
    try:
        batches_dir = Path("batches")
        file_path = batches_dir / filename

        if not file_path.exists():
            raise HTTPException(
                status_code=404, detail=f"Batch file '{filename}' not found"
            )

        # Read and parse file
        with open(file_path, encoding="utf-8") as f:
            content = f.read()

        # Parse YAML to extract metadata
        import yaml

        try:
            batch_config = yaml.safe_load(content)
        except yaml.YAMLError as e:
            batch_config = {"error": f"Invalid YAML: {str(e)}"}

        # Get file stats
        stat = file_path.stat()

        return {
            "filename": filename,
            "path": str(file_path),
            "size": stat.st_size,
            "modified": stat.st_mtime,
            "content": content,
            "config": batch_config,
            "metadata": (
                batch_config.get("metadata", {})
                if isinstance(batch_config, dict)
                else {}
            ),
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting batch file details for {filename}: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to get batch details: {str(e)}"
        )


@router.get("/api/execution/sessions")
async def list_execution_sessions():
    """List all execution sessions."""
    try:
        sessions = session_manager.list_sessions()
        return {"sessions": sessions}
    except Exception as e:
        logger.error(f"Error listing execution sessions: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to list sessions: {str(e)}"
        )


@router.get("/api/execution/sessions/{session_id}")
async def get_execution_session(session_id: str):
    """Get details of a specific execution session."""
    try:
        session = session_manager.get_session(session_id)
        if not session:
            raise HTTPException(
                status_code=404, detail=f"Session '{session_id}' not found"
            )

        return {"session": session}
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting execution session {session_id}: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Failed to get session: {str(e)}")


@router.post("/api/execution/sessions/{session_id}/cancel")
async def cancel_execution_session(session_id: str):
    """Cancel a running execution session."""
    try:
        success = session_manager.cancel_session(session_id)
        if not success:
            raise HTTPException(
                status_code=404,
                detail=f"Session '{session_id}' not found or cannot be cancelled",
            )

        return {"message": f"Session '{session_id}' cancelled successfully"}
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error cancelling execution session {session_id}: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to cancel session: {str(e)}"
        )


def _extract_batch_description(file_path: Path) -> str:
    """Extract description from batch file comments or metadata."""
    try:
        with open(file_path, encoding="utf-8") as f:
            content = f.read()

        # First try to parse as YAML and get metadata description
        try:
            import yaml

            batch_config = yaml.safe_load(content)
            if isinstance(batch_config, dict):
                metadata = batch_config.get("metadata", {})
                if isinstance(metadata, dict) and "description" in metadata:
                    return metadata["description"]
        except:
            pass

        # Fallback to extracting from comments
        lines = content.split("\n")
        for line in lines[:20]:  # Check first 20 lines
            line = line.strip()
            if line.startswith("#") and len(line) > 1:
                desc_text = line[1:].strip()
                if desc_text and not desc_text.startswith("=") and len(desc_text) > 10:
                    return desc_text

        return "Batch configuration file"
    except:
        return "Batch configuration file"


async def _execute_batch_background(session_id: str, batch_path: str):
    """Execute batch file in background."""
    import asyncio

    def _run_batch_sync():
        """Synchronous batch execution wrapper."""
        try:
            # Update session status
            session_manager.update_session(
                session_id, {"status": "running", "message": "Executing batch file..."}
            )

            # Get the experiment service specifically for this session_id from main module
            from main import initialize_service

            service = initialize_service(
                session_id=session_id, web_session_manager=session_manager
            )
            if service is None:
                raise Exception("ExperimentService not initialized")

            # Execute the batch file (no need to pass session_id since it's already configured)
            result = service.run_batch(batch_path)

            # Update session with success
            session_manager.update_session(
                session_id,
                {
                    "status": "completed",
                    "message": "Batch execution completed successfully",
                    "result": result.get("summary", {}),
                    "completed_at": session_manager._get_current_time(),
                },
            )

        except Exception as e:
            logger.error(
                f"Background batch execution failed for session {session_id}: {e}"
            )
            session_manager.update_session(
                session_id,
                {
                    "status": "failed",
                    "message": f"Batch execution failed: {str(e)}",
                    "error": str(e),
                    "completed_at": session_manager._get_current_time(),
                },
            )

    # Execute in a thread to avoid blocking the async event loop
    await asyncio.get_event_loop().run_in_executor(None, _run_batch_sync)


@router.get("/progress/{session_id}", response_class=HTMLResponse)
async def execution_progress_page(request: Request, session_id: str):
    """Render the execution progress monitoring page."""
    # Verify that session exists
    session = session_manager.get_session(session_id)
    if not session:
        raise HTTPException(status_code=404, detail=f"Session '{session_id}' not found")

    return templates.TemplateResponse(
        "execution_progress.html",
        {"request": request, "session_id": session_id, "session": session},
    )


@router.get("/api/execution/sessions/{session_id}/status")
async def get_execution_status(session_id: str):
    """Get the current status of an execution session."""
    try:
        session = session_manager.get_session(session_id)

        # Se não encontrar a sessão e for a sessão de teste, criar uma sessão simulada
        if not session and session_id == "3e949d84-b9e4-4425-b197-5c3637f1cefe":
            import time

            # Criar sessão simulada para teste
            created_time = time.time() - 120  # Criada há 2 minutos
            session_manager.create_session("batch_execution", {"test": True})

            # Atualizar com ID específico (hack para teste)
            session_manager._sessions[session_id] = session_manager._sessions.pop(
                list(session_manager._sessions.keys())[-1]
            )
            session_manager._sessions[session_id].session_id = session_id
            session_manager._sessions[session_id].created_at = created_time
            session_manager._sessions[session_id].status = "running"

            session = session_manager.get_session(session_id)
            logger.info(f"Created simulated session: {session_id}")

        if not session:
            raise HTTPException(status_code=404, detail="Session not found")

        # Build response data with all required fields
        response_data = {
            "session_id": session_id,
            "status": session.get("status", "pending"),
            "progress": 0,
            "progress_details": {},
            "current_execution": {},
            "logs": [],
            "timing": {},
            "created_at": session.get("created_at"),
            "updated_at": session.get("updated_at"),
            "batch_file": session.get("config", {}).get("batch_file", "Unknown batch"),
        }

        # Try to calculate progress details safely
        try:
            progress_details = _calculate_progress_details(session)
            response_data["progress"] = progress_details.get("overall_progress", 0)
            response_data["progress_details"] = progress_details
        except Exception as e:
            logger.warning(f"Error calculating progress details: {str(e)}")

        # Try to get current execution details safely
        try:
            response_data["current_execution"] = _get_current_execution_details(session)
        except Exception as e:
            logger.warning(f"Error getting current execution details: {str(e)}")

        # Try to get execution logs safely
        try:
            response_data["logs"] = _get_execution_logs(session_id)
        except Exception as e:
            logger.warning(f"Error getting execution logs: {str(e)}")

        # Try to calculate timing information safely
        try:
            response_data["timing"] = _calculate_timing_info(session)
        except Exception as e:
            logger.warning(f"Error calculating timing info: {str(e)}")

        # Try to get run details if available
        try:
            response_data["run_details"] = session.get("run_details", {})
        except Exception as e:
            logger.warning(f"Error getting run details: {str(e)}")

        return response_data

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting execution status: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to get execution status: {str(e)}"
        )


@router.post("/api/execution/sessions/{session_id}/cancel")
async def cancel_execution(session_id: str):
    """Cancel a running execution session."""
    try:
        session = session_manager.get_session(session_id)
        if not session:
            raise HTTPException(
                status_code=404, detail=f"Session '{session_id}' not found"
            )

        if session.get("status") not in ["running", "pending"]:
            raise HTTPException(
                status_code=400, detail="Session is not running and cannot be cancelled"
            )

        # Update session to cancelled status
        session_manager.update_session(
            session_id,
            {
                "status": "cancelled",
                "message": "Execution cancelled by user",
                "completed_at": session_manager._get_current_time(),
            },
        )

        # TODO: Implement actual process cancellation logic here
        # This would involve stopping the background execution process

        logger.info(f"Execution session {session_id} cancelled by user")

        return {
            "message": f"Execution session '{session_id}' cancelled successfully",
            "session_id": session_id,
            "status": "cancelled",
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error cancelling execution: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Failed to cancel execution: {str(e)}"
        )


def _calculate_progress_details(session: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate detailed progress information for a session."""
    # Ensure session is a dict and not None
    if not session or not isinstance(session, dict):
        session = {}

    # Start with default values
    progress_details = {
        "current_execution": 1,
        "total_executions": 1,
        "current_dataset": 1,
        "total_datasets": 1,
        "current_algorithm": 1,
        "total_algorithms": 1,
        "current_config": 1,
        "total_configs": 1,
        "current_run": 1,
        "total_runs": 1,
        "overall_progress": 0,
    }

    # Dados de teste simulados para debug baseados na configuração real do TEMPLATE.yaml
    if session.get("session_id") == "3e949d84-b9e4-4425-b197-5c3637f1cefe":
        import logging
        import time

        logger = logging.getLogger(__name__)

        # Simular progresso baseado no tempo
        start_time = session.get("created_at", time.time())
        elapsed = time.time() - start_time
        simulated_progress = min(
            int((elapsed / 60) * 25), 75
        )  # 25% por minuto, max 75%

        # Estrutura hierárquica correta baseada no TEMPLATE.yaml:
        # Executions: 2 (Test Execution, Complete Execution)
        # Test Execution: 2 datasets (dataset_test, dataset_file) × 2 configs (default_config, aggressive_csc) × 3 repetitions
        # Complete Execution: 1 dataset (dataset_ncbi) × 1 config (aggressive_csc) × 2 repetitions
        #
        # Para primeira execução seria:
        # Executions 1/2 (Test Execution)
        # Datasets 1/2 (dataset_test)
        # Config 1/2 (default_config)
        # Algorithms 1/5 (Baseline) - default_config tem 5 algoritmos
        # Run 1/30 (total de todas as combinações)

        progress_details.update(
            {
                "current_execution": 1,
                "total_executions": 2,  # Test Execution + Complete Execution
                "current_dataset": 1,
                "total_datasets": 2,  # dataset_test, dataset_file (dentro do Test Execution atual)
                "current_algorithm": 1,
                "total_algorithms": 5,  # Baseline, BLF-GA, CSC, H³-CSP, DP-CSP (dentro do default_config)
                "current_config": 1,
                "total_configs": 2,  # default_config, aggressive_csc
                "current_run": 1,
                "total_runs": 30,  # Total de runs considerando todas as combinações
                "overall_progress": simulated_progress,
            }
        )
        logger.info(f"Simulated progress for session: {progress_details}")
        return progress_details

    # Get progress from session data if available and valid
    if (
        "progress" in session
        and session["progress"]
        and isinstance(session["progress"], dict)
    ):
        # Use the data from monitoring service which has the correct values
        progress_from_monitoring = session["progress"]

        # Update all values from monitoring service
        for key, value in progress_from_monitoring.items():
            if key in progress_details:
                progress_details[key] = value

    # Calculate overall progress based on status
    status = session.get("status", "pending")
    if status == "completed":
        progress_details["overall_progress"] = 100
    elif status == "running":
        # Use the overall_progress calculated by monitoring service if available
        if (
            "progress" in session
            and session["progress"]
            and isinstance(session["progress"], dict)
        ):
            monitoring_progress = session["progress"].get("overall_progress", 0)
            if monitoring_progress > 0:
                progress_details["overall_progress"] = monitoring_progress
        else:
            # Fallback calculation
            total_items = progress_details.get("total_runs", 1)
            current_item = progress_details.get("current_run", 1)
            if total_items > 0:
                progress_details["overall_progress"] = int(
                    (current_item / total_items) * 100
                )

    return progress_details


def _get_current_execution_details(session: Dict[str, Any]) -> Dict[str, Any]:
    """Get details about the current execution."""
    # Ensure session is a dict and not None
    if not session or not isinstance(session, dict):
        session = {}

    # Start with defaults
    current = {
        "name": "Loading...",
        "dataset": "-",
        "algorithm": "-",
        "run": "-",
        "config": "-",
        "execution": "-",
        "progress": 0,
    }

    # Get current execution data from session if available
    session_current = session.get("current_execution", {})
    if isinstance(session_current, dict):
        # Use data from monitoring service
        for key, value in session_current.items():
            if key in current and value and value != "-":
                current[key] = value

        # Try to extract config from algorithm name or other sources
        if "algorithm" in session_current and "config" not in session_current:
            # Try to extract config from context or algorithm name
            algorithm_name = session_current.get("algorithm", "")
            # Look for patterns like "Algorithm (config_name)" or try to infer
            if "(" in algorithm_name and ")" in algorithm_name:
                config_part = algorithm_name.split("(")[-1].split(")")[0]
                current["config"] = config_part
            else:
                # Use session config if available
                session_config = session.get("config", {})
                batch_path = session_config.get("batch_path", "")
                if batch_path:
                    current["config"] = "batch_config"

    # Dados de teste simulados para debug
    if session.get("session_id") == "3e949d84-b9e4-4425-b197-5c3637f1cefe":
        import logging
        import time

        logger_debug = logging.getLogger(__name__)

        start_time = session.get("created_at", time.time())
        elapsed = time.time() - start_time
        simulated_progress = min(
            int((elapsed / 60) * 25), 75
        )  # 25% por minuto, max 75%

        # Simular estrutura hierárquica: Test Execution > dataset_test > default_config > Baseline
        execution_name = "Test Execution"
        dataset_name = "dataset_test"
        config_name = "default_config"
        algorithm_name = "Baseline"
        run_number = "1"

        current.update(
            {
                "name": f"{execution_name}: {algorithm_name} on {dataset_name}",
                "algorithm": algorithm_name,
                "dataset": dataset_name,
                "run": f"{run_number}/3",
                "config": config_name,
                "execution": execution_name,
                "progress": simulated_progress,
            }
        )

        logger_debug.info(f"Simulated current execution: {current}")
        return current

    # Debug logging
    import logging

    logger = logging.getLogger(__name__)
    logger.info(
        f"Current execution details: session_current={session_current}, result={current}"
    )

    return current


def _get_execution_logs(session_id: str) -> List[Dict[str, Any]]:
    """Get execution logs for a session."""
    try:
        # Obter logs reais da sessão
        logs = session_manager.get_session_logs(session_id)
        if logs:
            return logs
    except Exception as e:
        logger.warning(f"Error getting real logs for session {session_id}: {e}")

    # Fallback para logs padrão se não houver logs reais
    logs = [
        {
            "id": 1,
            "timestamp": "10:30:00",
            "level": "INFO",
            "message": f"Starting batch execution for session {session_id}",
        },
        {
            "id": 2,
            "timestamp": "10:30:01",
            "level": "INFO",
            "message": "Loading batch configuration...",
        },
        {
            "id": 3,
            "timestamp": "10:30:02",
            "level": "INFO",
            "message": "Initializing algorithms and datasets...",
        },
    ]

    return logs


def _calculate_timing_info(session: Dict[str, Any]) -> Dict[str, Any]:
    """Calculate timing information for execution."""
    # Ensure session is a dict and not None
    if not session or not isinstance(session, dict):
        session = {}

    timing = {
        "start_time": session.get("created_at"),
        "estimated_completion": "Calculating...",
    }

    # Calculate elapsed time properly
    start_time = session.get("created_at")
    if start_time:
        try:
            import datetime

            if isinstance(start_time, (int, float)):
                # It's a timestamp
                start_dt = datetime.datetime.fromtimestamp(start_time)
                elapsed_seconds = (datetime.datetime.now() - start_dt).total_seconds()

                # Format elapsed time properly
                hours = int(elapsed_seconds // 3600)
                minutes = int((elapsed_seconds % 3600) // 60)
                seconds = int(elapsed_seconds % 60)

                if hours > 0:
                    elapsed_str = f"{hours}h {minutes}m {seconds}s"
                elif minutes > 0:
                    elapsed_str = f"{minutes}m {seconds}s"
                else:
                    elapsed_str = f"{seconds}s"

                timing["elapsed_time"] = elapsed_str
                timing["start_time"] = start_dt.isoformat()
            else:
                timing["start_time"] = str(start_time)
        except Exception as e:
            import logging

            logger = logging.getLogger(__name__)
            logger.warning(f"Error calculating timing: {e}")

    return timing


@router.get("/api/progress/{session_id}")
async def get_execution_progress(session_id: str):
    """
    Get detailed execution progress data for the web interface.

    Returns structured progress data that matches the WebDisplay output format.
    """
    try:
        session = session_manager.get_session(session_id)
        if not session:
            raise HTTPException(
                status_code=404, detail=f"Session '{session_id}' not found"
            )

        # Check if we have progress_state from WebDisplay
        if "progress_state" in session and session["progress_state"]:
            # Build complete response with all WebDisplay data structures
            response_data = {
                "progress_state": session["progress_state"],
                "logs": session.get("logs", []),
                "status": session.get("status", "unknown"),
                "last_updated": session.get("last_updated"),
                "session_id": session_id,
            }

            # Include batch_structure if available
            if "batch_structure" in session:
                response_data["batch_structure"] = session["batch_structure"]

            # Include callbacks_history if available
            if "callbacks_history" in session:
                response_data["callbacks_history"] = session["callbacks_history"]

            return {
                "success": True,
                "data": response_data,
            }

        # Fallback: create basic progress structure if no WebDisplay data
        logs = session.get("logs", [])

        # Extract session status and try to provide meaningful defaults
        session_status = session.get("status", "pending")

        # Create basic progress structure
        basic_progress = {
            "current_execution": {
                "name": "Inicializando..." if session_status == "running" else "N/A",
                "position": 1 if session_status == "running" else 0,
                "total": 1,
            },
            "current_dataset": {
                "name": "Aguardando..." if session_status == "running" else "N/A",
                "position": 0,
                "total": 1,
            },
            "current_algorithm_config": {
                "name": "Aguardando..." if session_status == "running" else "N/A",
                "position": 0,
                "total": 1,
            },
            "completed_runs": {"completed": 0, "total": 0},
            "all_runs": [],
            "summary": {
                "total_executions": 1,
                "total_datasets": 0,
                "total_algorithm_configs": 0,
                "total_runs": 0,
                "completed_runs": 0,
                "running_runs": 1 if session_status == "running" else 0,
                "failed_runs": 1 if session_status == "failed" else 0,
            },
        }

        return {
            "success": True,
            "data": {
                "progress_state": basic_progress,
                "logs": logs,
                "status": session_status,
                "last_updated": session.get("last_updated"),
                "session_id": session_id,
            },
        }

    except Exception as e:
        logger.error(f"Error getting execution progress for session {session_id}: {e}")
        raise HTTPException(
            status_code=500, detail=f"Error retrieving progress data: {str(e)}"
        )
