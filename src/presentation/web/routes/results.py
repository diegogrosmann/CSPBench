"""
Results management endpoints.
"""

import json
import logging
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse, StreamingResponse
from pydantic import BaseModel

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/results", tags=["results"])


class ResultSession(BaseModel):
    """Session result model."""
    session_id: str
    timestamp: datetime
    algorithm: Optional[str] = None
    dataset_name: Optional[str] = None
    status: str
    execution_time: Optional[float] = None
    best_string: Optional[str] = None
    max_distance: Optional[float] = None
    download_url: str
    files: List[str] = []
    size_mb: Optional[float] = None


class ResultsListResponse(BaseModel):
    """Response model for results listing."""
    sessions: List[ResultSession]
    total_count: int
    total_size_mb: float


@router.get("/", response_model=ResultsListResponse)
async def list_results(
    limit: int = Query(50, ge=1, le=100),
    offset: int = Query(0, ge=0),
    algorithm: Optional[str] = None,
    status: Optional[str] = None,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
):
    """List all execution results with filtering and pagination."""
    try:
        # Get results from outputs directory
        outputs_dir = Path("outputs")
        sessions = []
        total_size = 0

        if outputs_dir.exists():
            # Get all session directories
            session_dirs = [d for d in outputs_dir.iterdir() if d.is_dir()]
            
            # Sort by timestamp (newest first)
            session_dirs.sort(reverse=True)

            for session_dir in session_dirs:
                try:
                    session_info = await _extract_session_info(session_dir)
                    
                    # Apply filters
                    if algorithm and session_info.algorithm != algorithm:
                        continue
                    if status and session_info.status != status:
                        continue
                    if date_from:
                        from_date = datetime.fromisoformat(date_from.replace('Z', '+00:00'))
                        if session_info.timestamp < from_date:
                            continue
                    if date_to:
                        to_date = datetime.fromisoformat(date_to.replace('Z', '+00:00'))
                        if session_info.timestamp > to_date:
                            continue

                    sessions.append(session_info)
                    total_size += session_info.size_mb or 0

                except Exception as e:
                    logger.warning(f"Failed to extract info from {session_dir}: {e}")
                    continue

        # Apply pagination
        paginated_sessions = sessions[offset:offset + limit]

        return ResultsListResponse(
            sessions=paginated_sessions,
            total_count=len(sessions),
            total_size_mb=round(total_size, 2)
        )

    except Exception as e:
        logger.error(f"Error listing results: {e}")
        raise HTTPException(status_code=500, detail="Failed to list results")


@router.get("/{session_id}")
async def get_result_details(session_id: str):
    """Get detailed information about a specific result session."""
    try:
        session_dir = await _find_session_directory(session_id)
        if not session_dir:
            raise HTTPException(status_code=404, detail="Session not found")

        session_info = await _extract_session_info(session_dir)
        
        # Load detailed results if available
        results_files = list(session_dir.glob("results_*.json"))
        if results_files:
            with open(results_files[0], 'r') as f:
                detailed_results = json.load(f)
            session_info_dict = session_info.dict()
            session_info_dict['detailed_results'] = detailed_results
            return session_info_dict

        return session_info

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting result details for {session_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to get result details")


@router.get("/{session_id}/download")
async def download_session_results(session_id: str):
    """Download results for a specific session as ZIP."""
    try:
        session_dir = await _find_session_directory(session_id)
        if not session_dir:
            raise HTTPException(status_code=404, detail="Session not found")

        # Create ZIP file with all session results
        zip_path = await _create_session_zip(session_id, session_dir)
        
        return FileResponse(
            path=zip_path,
            filename=f"cspbench_results_{session_id}.zip",
            media_type="application/zip",
            headers={"Content-Disposition": f"attachment; filename=cspbench_results_{session_id}.zip"}
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error downloading session {session_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to download session results")


@router.get("/{session_id}/files/{file_type}")
async def download_session_file(session_id: str, file_type: str):
    """Download a specific file type from a session (json, csv, log, etc.)."""
    try:
        session_dir = await _find_session_directory(session_id)
        if not session_dir:
            raise HTTPException(status_code=404, detail="Session not found")

        # Map file types to patterns
        file_patterns = {
            "json": "results_*.json",
            "csv": "results_*.csv", 
            "log": "cspbench.log",
            "report": "report/*",
            "plots": "plots/*"
        }

        if file_type not in file_patterns:
            raise HTTPException(status_code=400, detail=f"Unsupported file type: {file_type}")

        # Find files matching the pattern
        files = list(session_dir.glob(file_patterns[file_type]))
        if not files:
            raise HTTPException(status_code=404, detail=f"No {file_type} files found")

        # If single file, return it directly
        if len(files) == 1 and files[0].is_file():
            return FileResponse(
                path=files[0],
                filename=files[0].name,
                headers={"Content-Disposition": f"attachment; filename={files[0].name}"}
            )

        # If multiple files or directory, create ZIP
        zip_path = await _create_files_zip(session_id, files, file_type)
        return FileResponse(
            path=zip_path,
            filename=f"cspbench_{session_id}_{file_type}.zip",
            media_type="application/zip",
            headers={"Content-Disposition": f"attachment; filename=cspbench_{session_id}_{file_type}.zip"}
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error downloading {file_type} for session {session_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to download {file_type}")


@router.post("/download/batch")
async def download_batch_results(session_ids: List[str]):
    """Download results for multiple sessions as a single ZIP."""
    try:
        if len(session_ids) > 20:  # Limit batch size
            raise HTTPException(status_code=400, detail="Batch size limited to 20 sessions")

        # Find all session directories
        valid_sessions = []
        for session_id in session_ids:
            session_dir = await _find_session_directory(session_id)
            if session_dir:
                valid_sessions.append((session_id, session_dir))

        if not valid_sessions:
            raise HTTPException(status_code=404, detail="No valid sessions found")

        # Create batch ZIP
        batch_zip_path = await _create_batch_zip(valid_sessions)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"cspbench_batch_results_{timestamp}.zip"
        
        return FileResponse(
            path=batch_zip_path,
            filename=filename,
            media_type="application/zip",
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error creating batch download: {e}")
        raise HTTPException(status_code=500, detail="Failed to create batch download")


@router.delete("/{session_id}")
async def delete_session_results(session_id: str):
    """Delete results for a specific session."""
    try:
        session_dir = await _find_session_directory(session_id)
        if not session_dir:
            raise HTTPException(status_code=404, detail="Session not found")

        # Remove session directory
        import shutil
        shutil.rmtree(session_dir)
        
        return {"message": f"Session {session_id} deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error deleting session {session_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to delete session")


# Helper functions

async def _find_session_directory(session_id: str) -> Optional[Path]:
    """Find the directory for a given session ID."""
    outputs_dir = Path("outputs")
    if not outputs_dir.exists():
        return None

    # Look for exact match first
    exact_match = outputs_dir / session_id
    if exact_match.exists():
        return exact_match

    # Look for session ID in directory names or files
    for session_dir in outputs_dir.iterdir():
        if session_dir.is_dir():
            # Check if any result file contains this session ID
            for file_path in session_dir.glob("results_*.json"):
                try:
                    with open(file_path, 'r') as f:
                        data = json.load(f)
                        if data.get('session_id') == session_id:
                            return session_dir
                except Exception:
                    continue

    return None


async def _extract_session_info(session_dir: Path) -> ResultSession:
    """Extract session information from directory."""
    session_id = session_dir.name
    
    # Try to parse timestamp from directory name
    try:
        timestamp = datetime.strptime(session_id, "%Y%m%d_%H%M%S")
    except ValueError:
        # Fallback to directory modification time
        timestamp = datetime.fromtimestamp(session_dir.stat().st_mtime)

    # Initialize session info
    session_info = {
        "session_id": session_id,
        "timestamp": timestamp,
        "status": "completed",
        "download_url": f"/api/results/{session_id}/download",
        "files": []
    }

    # Get file list and calculate size
    total_size = 0
    for file_path in session_dir.rglob("*"):
        if file_path.is_file():
            session_info["files"].append(str(file_path.relative_to(session_dir)))
            total_size += file_path.stat().st_size

    session_info["size_mb"] = round(total_size / (1024 * 1024), 2)

    # Try to extract additional info from results files
    results_files = list(session_dir.glob("results_*.json"))
    if results_files:
        try:
            with open(results_files[0], 'r') as f:
                results_data = json.load(f)
                
            session_info.update({
                "algorithm": results_data.get("algorithm"),
                "dataset_name": results_data.get("dataset_name"),
                "execution_time": results_data.get("execution_time"),
                "best_string": results_data.get("best_string"),
                "max_distance": results_data.get("max_distance")
            })
        except Exception as e:
            logger.warning(f"Failed to parse results file: {e}")

    return ResultSession(**session_info)


async def _create_session_zip(session_id: str, session_dir: Path) -> Path:
    """Create ZIP file containing all files from a session directory."""
    temp_dir = Path(tempfile.gettempdir())
    zip_path = temp_dir / f"session_{session_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"

    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file_path in session_dir.rglob("*"):
            if file_path.is_file():
                arc_name = str(file_path.relative_to(session_dir))
                zipf.write(file_path, arc_name)

    return zip_path


async def _create_files_zip(session_id: str, files: List[Path], file_type: str) -> Path:
    """Create ZIP file containing specific files."""
    temp_dir = Path(tempfile.gettempdir())
    zip_path = temp_dir / f"files_{session_id}_{file_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"

    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file_path in files:
            if file_path.is_file():
                zipf.write(file_path, file_path.name)
            elif file_path.is_dir():
                for sub_file in file_path.rglob("*"):
                    if sub_file.is_file():
                        arc_name = str(sub_file.relative_to(file_path.parent))
                        zipf.write(sub_file, arc_name)

    return zip_path


async def _create_batch_zip(sessions: List[tuple]) -> Path:
    """Create ZIP file containing multiple sessions."""
    temp_dir = Path(tempfile.gettempdir())
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    zip_path = temp_dir / f"batch_results_{timestamp}.zip"

    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for session_id, session_dir in sessions:
            for file_path in session_dir.rglob("*"):
                if file_path.is_file():
                    arc_name = f"{session_id}/{file_path.relative_to(session_dir)}"
                    zipf.write(file_path, arc_name)

    return zip_path
