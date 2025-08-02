"""
Algorithm execution endpoints.
"""

import asyncio
import json
import logging
import tempfile
import uuid
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Dict

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from slowapi import Limiter
from slowapi.util import get_remote_address

from src.domain.algorithms import global_registry

from ..core.models import ExecutionRequest, ExecutionResult
from ..core.security import SecurityValidator

logger = logging.getLogger(__name__)

# Rate limiting
limiter = Limiter(key_func=get_remote_address)

router = APIRouter(prefix="/api/execution", tags=["execution"])


@router.post("/", response_model=ExecutionResult)
@limiter.limit("10/minute")  # Limit execution requests
async def execute_algorithm(request: ExecutionRequest):
    """Execute algorithm with security context."""
    session_id = str(uuid.uuid4())

    try:
        # Validate algorithm exists
        algorithm_cls = global_registry.get(request.algorithm)
        if not algorithm_cls:
            raise HTTPException(
                status_code=400, detail=f"Algorithm '{request.algorithm}' not found"
            )

        # Validate dataset content
        if not SecurityValidator.validate_dataset_content(request.dataset_content):
            raise HTTPException(status_code=400, detail="Invalid dataset content")

        # Get algorithm metadata for parameter validation
        temp_instance = algorithm_cls()
        default_params = getattr(temp_instance, "default_params", {})

        # Validate parameters
        validated_params = SecurityValidator.validate_algorithm_parameters(
            request.parameters, default_params
        )

        # Sanitize filename
        safe_filename = SecurityValidator.sanitize_filename(request.dataset_name)

        # Create secure session directory
        session_dir = Path(tempfile.mkdtemp(prefix=f"exec_{session_id}_"))

        # Execute with security context
        result = await _execute_with_security_context(
            request, validated_params, safe_filename, session_id, session_dir
        )

        return result

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Execution failed for session {session_id}: {e}")
        return ExecutionResult(
            session_id=session_id, status="failed", error=f"Execution failed: {str(e)}"
        )


async def _execute_with_security_context(
    execution_request: ExecutionRequest,
    validated_params: Dict,
    safe_filename: str,
    session_id: str,
    session_dir: Path,
) -> ExecutionResult:
    """Execute algorithm with comprehensive security context."""
    try:
        # 1. Create secure dataset file
        dataset_path = session_dir / safe_filename
        with open(dataset_path, "w", encoding="utf-8") as f:
            f.write(execution_request.dataset_content)

        # 2. Load dataset using infrastructure layer
        from src.infrastructure.persistence.dataset_repository import (
            FileDatasetRepository,
        )

        dataset_repo = FileDatasetRepository(str(session_dir))
        dataset = dataset_repo.load(safe_filename.replace(".fasta", ""))

        # 3. Get algorithm class through registry
        algorithm_cls = global_registry.get(execution_request.algorithm)
        if not algorithm_cls:
            raise ValueError(f"Algorithm not found: {execution_request.algorithm}")

        # 4. Prepare execution parameters
        execution_params = {
            **validated_params,
            "save_history": execution_request.save_history,
            "session_id": session_id,
            "work_dir": str(session_dir),
        }

        # 5. Execute algorithm with timeout and monitoring
        logger.info(f"Starting secure execution for session {session_id}")
        start_time = datetime.now()

        # Create algorithm instance
        algorithm_instance = algorithm_cls()

        # Execute with timeout protection
        def run_algorithm():
            return algorithm_instance.run(dataset, **execution_params)

        # Run in executor with timeout
        best_string, max_distance, metadata = await asyncio.wait_for(
            asyncio.get_event_loop().run_in_executor(None, run_algorithm),
            timeout=execution_request.timeout,
        )

        end_time = datetime.now()
        execution_time = (end_time - start_time).total_seconds()

        # 6. Enhance metadata with security information
        enhanced_metadata = {
            **metadata,
            "execution_context": {
                "session_id": session_id,
                "execution_time": execution_time,
                "algorithm_version": getattr(algorithm_cls, "__version__", "1.0.0"),
                "cspbench_version": "0.1.0",
                "secure_execution": True,
            },
        }

        # 7. Prepare comprehensive results
        results = {
            "algorithm": execution_request.algorithm,
            "best_string": best_string,
            "max_distance": max_distance,
            "metadata": enhanced_metadata,
            "execution_time": execution_time,
            "session_id": session_id,
            "dataset_info": {
                "num_strings": len(dataset.sequences),
                "string_length": len(dataset.sequences[0]) if dataset.sequences else 0,
                "alphabet": dataset.alphabet,
            },
            "security_info": {
                "sanitized_filename": safe_filename,
                "parameters_validated": True,
                "secure_session": True,
            },
        }

        # 8. Create secure results package
        zip_path = await _create_secure_results_zip(session_id, results, session_dir)

        logger.info(f"Execution completed successfully for session {session_id}")

        return ExecutionResult(
            session_id=session_id,
            status="completed",
            result=results,
            download_url=f"/api/execution/download/{session_id}",
        )

    except asyncio.TimeoutError:
        logger.warning(f"Execution timeout for session {session_id}")
        return ExecutionResult(
            session_id=session_id,
            status="timeout",
            error=f"Execution timed out after {execution_request.timeout} seconds",
        )
    except Exception as e:
        logger.error(f"Secure execution failed for session {session_id}: {e}")

        # Clean up on failure
        try:
            import shutil

            if session_dir.exists():
                shutil.rmtree(session_dir)
        except Exception as cleanup_error:
            logger.warning(f"Failed to cleanup session {session_id}: {cleanup_error}")

        return ExecutionResult(
            session_id=session_id, status="failed", error=f"Execution failed: {str(e)}"
        )


async def _create_secure_results_zip(
    session_id: str, results: Dict, session_dir: Path
) -> Path:
    """Create secure ZIP file with execution results."""
    zip_path = session_dir / f"results_{session_id}.zip"

    try:
        with zipfile.ZipFile(
            zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6
        ) as zipf:
            # Add main results as JSON
            results_json = json.dumps(results, indent=2, default=str)
            zipf.writestr("results.json", results_json)

            # Add execution summary
            summary = f"""CSPBench Execution Results
============================

Session ID: {session_id}
Algorithm: {results.get('algorithm', 'Unknown')}
Best String: {results.get('best_string', 'N/A')}
Maximum Distance: {results.get('max_distance', 'N/A')}
Execution Time: {results.get('execution_time', 'N/A')}s

Dataset Information:
- Number of strings: {results.get('dataset_info', {}).get('num_strings', 'N/A')}
- String length: {results.get('dataset_info', {}).get('string_length', 'N/A')}
- Alphabet: {results.get('dataset_info', {}).get('alphabet', 'N/A')}

Security Information:
- Secure execution: {results.get('security_info', {}).get('secure_session', False)}
- Parameters validated: {results.get('security_info', {}).get('parameters_validated', False)}
- Sanitized filename: {results.get('security_info', {}).get('sanitized_filename', 'N/A')}

Metadata:
{json.dumps(results.get('metadata', {}), indent=2, default=str)}
"""
            zipf.writestr("execution_summary.txt", summary)

            # Add algorithm-specific outputs if they exist
            output_dir = session_dir / "algorithm_output"
            if output_dir.exists():
                for file_path in output_dir.rglob("*"):
                    if file_path.is_file():
                        arc_name = (
                            f"algorithm_output/{file_path.relative_to(output_dir)}"
                        )
                        zipf.write(file_path, arc_name)

            # Add execution log if it exists
            log_file = session_dir / "execution.log"
            if log_file.exists():
                zipf.write(log_file, "execution.log")

    except Exception as e:
        logger.error(f"Failed to create results ZIP for session {session_id}: {e}")
        raise

    return zip_path


@router.get("/download/{session_id}")
async def download_results(session_id: str):
    """Download execution results ZIP file."""
    try:
        # Find the session directory (this would need proper session management)
        import tempfile

        temp_dir = Path(tempfile.gettempdir())

        # Look for session directory
        session_dirs = list(temp_dir.glob(f"exec_{session_id}_*"))
        if not session_dirs:
            raise HTTPException(status_code=404, detail="Session not found")

        session_dir = session_dirs[0]
        zip_file = session_dir / f"results_{session_id}.zip"

        if not zip_file.exists():
            raise HTTPException(status_code=404, detail="Results file not found")

        return FileResponse(
            path=zip_file,
            filename=f"cspbench_results_{session_id}.zip",
            media_type="application/zip",
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error downloading results for session {session_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to download results")
