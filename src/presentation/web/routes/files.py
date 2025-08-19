"""
API routes for file operations and downloads
"""

from fastapi import APIRouter, HTTPException, Response
from fastapi.responses import FileResponse, StreamingResponse
from typing import Optional
import os
import zipfile
import io
from pathlib import Path

from src.infrastructure.logging_config import get_logger

logger = get_logger(__name__)
router = APIRouter()


@router.get("/api/files/download")
async def download_file(path: str):
    """Download a specific file"""
    try:
        # Sanitize path to prevent directory traversal
        file_path = Path(path).resolve()
        
        # Ensure file is within allowed directories (outputs, logs)
        allowed_dirs = [Path("outputs").resolve(), Path("logs").resolve()]
        if not any(str(file_path).startswith(str(allowed_dir)) for allowed_dir in allowed_dirs):
            raise HTTPException(status_code=403, detail="Access denied to this path")
        
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="File not found")
        
        if not file_path.is_file():
            raise HTTPException(status_code=400, detail="Path is not a file")
        
        return FileResponse(
            path=str(file_path),
            filename=file_path.name,
            media_type='application/octet-stream'
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error downloading file {path}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/files/download-zip/{work_id}")
async def download_execution_zip(work_id: str):
    """Download all files from an execution as a ZIP archive"""
    try:
        # Check if work directory exists
        work_dir = Path("outputs") / work_id
        if not work_dir.exists():
            raise HTTPException(status_code=404, detail=f"Execution {work_id} not found")
        
        # Create zip file in memory
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            # Add all files from the work directory
            for file_path in work_dir.rglob("*"):
                if file_path.is_file():
                    # Create relative path for the zip
                    relative_path = file_path.relative_to(work_dir)
                    zip_file.write(file_path, relative_path)
        
        zip_buffer.seek(0)
        
        # Return zip file
        return StreamingResponse(
            io.BytesIO(zip_buffer.read()),
            media_type="application/zip",
            headers={"Content-Disposition": f"attachment; filename=execution_{work_id}_results.zip"}
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error creating zip for execution {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/files/list/{work_id}")
async def list_execution_files(work_id: str):
    """List all files in an execution directory"""
    try:
        # Check if work directory exists
        work_dir = Path("outputs") / work_id
        if not work_dir.exists():
            raise HTTPException(status_code=404, detail=f"Execution {work_id} not found")
        
        files = {
            "logs": [],
            "outputs": [],
            "all": []
        }
        
        for file_path in work_dir.rglob("*"):
            if file_path.is_file():
                file_info = {
                    "name": file_path.name,
                    "path": str(file_path),
                    "size": file_path.stat().st_size,
                    "modified": file_path.stat().st_mtime
                }
                
                files["all"].append(file_info)
                
                # Categorize files
                if file_path.suffix in ['.log', '.txt']:
                    files["logs"].append(file_info)
                elif file_path.suffix in ['.json', '.csv', '.xlsx', '.pdf']:
                    files["outputs"].append(file_info)
        
        return files
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error listing files for execution {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))
