"""
File management routes for accessing work results and downloading files.

This module provides secure file access endpoints for downloading execution
results, creating ZIP archives, and managing file access with proper security
validation and path sanitization.
"""

import io
import logging
import zipfile
from pathlib import Path
from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, StreamingResponse

from src.infrastructure.utils.path_utils import get_output_base_directory

logger = logging.getLogger(__name__)

router = APIRouter(tags=["files"])


def format_file_size(size_bytes: int) -> str:
    """Format file size in human readable format.

    Converts byte values to human-readable format with appropriate
    units (B, KB, MB, GB) for display purposes.

    Args:
        size_bytes (int): File size in bytes.

    Returns:
        str: Human-readable file size with units.

    Example:
        >>> format_file_size(1024)
        "1.0 KB"
        >>> format_file_size(1048576)
        "1.0 MB"
    """
    if size_bytes == 0:
        return "0 B"

    size_names = ["B", "KB", "MB", "GB"]
    i = 0
    while size_bytes >= 1024 and i < len(size_names) - 1:
        size_bytes /= 1024.0
        i += 1

    return f"{size_bytes:.1f} {size_names[i]}"


from src.infrastructure.logging_config import get_logger

logger = get_logger(__name__)
router = APIRouter()


@router.get("/api/files/download")
async def download_file(path: str):
    """Download a specific file with security validation.

    Provides secure file download with path validation to prevent
    directory traversal attacks and unauthorized file access.

    Args:
        path (str): File path to download.

    Returns:
        FileResponse: File download response with proper headers.

    Raises:
        HTTPException: If path is invalid, file not found, or access denied.

    Security:
        - Validates path against allowed directories
        - Prevents directory traversal attacks
        - Ensures file exists and is accessible
    """
    try:
        # Sanitize path to prevent directory traversal
        file_path = Path(path).resolve()

        # Ensure file is within allowed directories (outputs, logs)
        output_base = get_output_base_directory()
        allowed_dirs = [output_base.resolve(), Path("logs").resolve()]
        if not any(
            str(file_path).startswith(str(allowed_dir)) for allowed_dir in allowed_dirs
        ):
            raise HTTPException(status_code=403, detail="Access denied to this path")

        if not file_path.exists():
            raise HTTPException(status_code=404, detail="File not found")

        if not file_path.is_file():
            raise HTTPException(status_code=400, detail="Path is not a file")

        return FileResponse(
            path=str(file_path),
            filename=file_path.name,
            media_type="application/octet-stream",
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error downloading file {path}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/files/download-zip/{work_id}")
async def download_execution_zip(work_id: str):
    """Download all files from an execution as a ZIP archive.

    Creates a compressed ZIP archive containing all files from
    a work execution directory for convenient bulk download.

    Args:
        work_id (str): Unique identifier of the work execution.

    Returns:
        StreamingResponse: ZIP file download stream with appropriate headers.

    Raises:
        HTTPException: If work directory not found or ZIP creation fails.

    Note:
        Files are organized in the ZIP maintaining their relative
        directory structure from the work directory.
    """
    try:
        # Check if work directory exists
        work_dir = get_output_base_directory() / work_id
        if not work_dir.exists():
            raise HTTPException(
                status_code=404, detail=f"Execution {work_id} not found"
            )

        # Create zip file in memory
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
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
            headers={
                "Content-Disposition": f"attachment; filename=execution_{work_id}_results.zip"
            },
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error creating zip for execution {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/files/list/{work_id}")
async def list_execution_files(work_id: str):
    """List all files in an execution directory with categorization.

    Provides a comprehensive file listing with categorization
    (logs, outputs, all) and metadata for display and selection.

    Args:
        work_id (str): Unique identifier of the work execution.

    Returns:
        dict: Categorized file listings with metadata.

    Raises:
        HTTPException: If work directory not found or listing fails.

    Note:
        Files are categorized based on extensions:
        - Logs: .log, .txt files
        - Outputs: .json, .csv, .xlsx, .pdf files
    """
    try:
        # Check if work directory exists
        work_dir = get_output_base_directory() / work_id
        if not work_dir.exists():
            raise HTTPException(
                status_code=404, detail=f"Execution {work_id} not found"
            )

        files = {"logs": [], "outputs": [], "all": []}

        for file_path in work_dir.rglob("*"):
            if file_path.is_file():
                file_info = {
                    "name": file_path.name,
                    "path": str(file_path),
                    "size": file_path.stat().st_size,
                    "modified": file_path.stat().st_mtime,
                }

                files["all"].append(file_info)

                # Categorize files
                if file_path.suffix in [".log", ".txt"]:
                    files["logs"].append(file_info)
                elif file_path.suffix in [".json", ".csv", ".xlsx", ".pdf"]:
                    files["outputs"].append(file_info)

        return files

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error listing files for execution {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/files/works")
async def list_works() -> List[Dict[str, Any]]:
    """List available work result directories under OUTPUT_BASE_DIRECTORY.

    Returns basic metadata for all work directories to enable
    UI selection and navigation of available results.

    Returns:
        List[Dict[str, Any]]: List of work directories with metadata.

    Note:
        Includes manifest information if available for enhanced
        metadata display and organization.
    """
    outputs_dir = get_output_base_directory()
    if not outputs_dir.exists():
        return []
    works: List[Dict[str, Any]] = []
    for child in sorted(outputs_dir.iterdir()):
        if child.is_dir():
            try:
                # Try to read manifest for quick metadata
                manifest_path = child / "manifest.json"
                manifest: Dict[str, Any] | None = None
                if manifest_path.exists():
                    import json

                    try:
                        manifest = json.loads(manifest_path.read_text())
                    except Exception:  # noqa: BLE001
                        manifest = None
                works.append(
                    {
                        "work_id": child.name,
                        "path": str(child),
                        "modified": child.stat().st_mtime,
                        "manifest": manifest,
                    }
                )
            except Exception as e:  # noqa: BLE001
                logger.warning(f"Failed to read work dir {child}: {e}")
    return works


@router.post("/api/files/download-zip-selected/{work_id}")
async def download_selected_files_zip(work_id: str, paths: List[str]):
    """Download a subset of files from a work directory as a ZIP.

    Creates a ZIP archive containing only selected files from
    a work directory with security validation and path checking.

    Args:
        work_id (str): Work directory identifier.
        paths (List[str]): List of file paths to include in ZIP.

    Returns:
        StreamingResponse: ZIP file containing selected files.

    Raises:
        HTTPException: If work not found, paths invalid, or ZIP creation fails.

    Security:
        - Validates all paths are within work directory
        - Prevents directory traversal attacks
        - Ensures all files exist before creating ZIP
    """
    try:
        base_dir = (get_output_base_directory() / work_id).resolve()
        if not base_dir.exists():
            raise HTTPException(status_code=404, detail="Work not found")

        # Normalize and validate paths
        normalized: List[Path] = []
        for p in paths:
            p_path = Path(p)
            if not p_path.is_absolute():
                p_path = base_dir / p_path
            p_path = p_path.resolve()
            if not str(p_path).startswith(str(base_dir)):
                raise HTTPException(
                    status_code=400, detail=f"Invalid path outside work dir: {p}"
                )
            if not p_path.exists() or not p_path.is_file():
                raise HTTPException(status_code=404, detail=f"File not found: {p}")
            normalized.append(p_path)

        if not normalized:
            raise HTTPException(status_code=400, detail="No files selected")

        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for file_path in normalized:
                rel = file_path.relative_to(base_dir)
                zip_file.write(file_path, rel)
        zip_buffer.seek(0)

        return StreamingResponse(
            io.BytesIO(zip_buffer.read()),
            media_type="application/zip",
            headers={
                "Content-Disposition": f"attachment; filename=execution_{work_id}_selected.zip"
            },
        )
    except HTTPException:
        raise
    except Exception as e:  # noqa: BLE001
        logger.error(f"Error creating selected zip for {work_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))
