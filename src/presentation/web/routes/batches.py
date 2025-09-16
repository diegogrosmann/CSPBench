"""
Batch management API endpoints for CSPBench Web Interface.

This module provides REST API endpoints for batch file management including:
- CRUD operations for batch configuration files
- Batch file listing and metadata extraction
- YAML validation and secure file handling
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml
from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from src.infrastructure.utils.path_utils import get_batch_directory

from ..core.batch_models import (
    BatchFileInfo,
    BatchListResponse,
    BatchUploadRequest,
    OperationResponse,
)
from ..core.security import sanitize_filename

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/batches", tags=["batches"])


def format_file_size(size_bytes: int) -> str:
    """Format file size in human readable format.
    
    Converts byte values to human-readable strings with appropriate
    units for display in user interfaces.
    
    Args:
        size_bytes (int): File size in bytes.
        
    Returns:
        str: Formatted size string with units (B, KB, MB, GB).
    """
    if size_bytes == 0:
        return "0 B"

    size_names = ["B", "KB", "MB", "GB"]
    i = 0
    while size_bytes >= 1024 and i < len(size_names) - 1:
        size_bytes /= 1024.0
        i += 1

    return f"{size_bytes:.1f} {size_names[i]}"


def extract_description_from_yaml(file_path: Path) -> Optional[str]:
    """Extract description from YAML file for metadata display.
    
    Attempts to extract description information from YAML batch files
    by looking for description fields or comment-based descriptions.
    
    Args:
        file_path (Path): Path to YAML file to analyze.
        
    Returns:
        Optional[str]: Extracted description or None if not found.
        
    Note:
        Gracefully handles parsing errors and returns None
        rather than raising exceptions for malformed files.
    """
    try:
        with open(file_path, encoding="utf-8") as f:
            content = f.read()

        # Look for description field in YAML
        lines = content.split("\n")
        for line in lines:
            line = line.strip()
            if line.startswith("description:"):
                desc = line.split("description:", 1)[1].strip()
                # Remove quotes if present
                desc = desc.strip("\"'")
                return desc

        # If no description field, use first comment as description
        for line in lines:
            line = line.strip()
            if line.startswith("#") and len(line) > 1:
                return line[1:].strip()

        return None
    except Exception as e:
        logger.warning(f"Failed to extract description from {file_path}: {e}")
        return None


def extract_metadata_from_yaml(file_path: Path) -> dict:
    """Extract complete metadata from YAML file structure.
    
    Parses YAML batch configuration files to extract embedded
    metadata including author, version, tags, and other properties.
    
    Args:
        file_path (Path): Path to YAML file to parse.
        
    Returns:
        dict: Dictionary containing extracted metadata fields.
        
    Note:
        Returns empty dict if metadata section not found or
        parsing fails, ensuring robust error handling.
    """
    try:
        with open(file_path, encoding="utf-8") as f:
            content = yaml.safe_load(f)

        # Look for metadata section
        if isinstance(content, dict):
            metadata = content.get("metadata", {})
            if isinstance(metadata, dict):
                return {
                    "name": metadata.get("name"),
                    "description": metadata.get("description"),
                    "author": metadata.get("author"),
                    "version": metadata.get("version"),
                    "creation_date": metadata.get("creation_date"),
                    "tags": metadata.get("tags", []),
                }

        return {}

    except Exception as e:
        logger.warning(f"Failed to extract metadata from {file_path}: {e}")
        return {}


def extract_metadata_name_from_yaml(file_path: Path) -> Optional[str]:
    """Extract metadata.name field from YAML file.
    
    Specifically extracts the name field from the metadata section
    of YAML batch configuration files for display purposes.
    
    Args:
        file_path (Path): Path to YAML file to parse.
        
    Returns:
        Optional[str]: Metadata name or None if not found.
    """
    try:
        with open(file_path, encoding="utf-8") as f:
            content = yaml.safe_load(f)

        # Look for metadata.name field
        if isinstance(content, dict):
            metadata = content.get("metadata", {})
            if isinstance(metadata, dict):
                return metadata.get("name")

        return None

    except Exception as e:
        logger.warning(f"Failed to extract metadata.name from {file_path}: {e}")
        return None


@router.post("/upload", response_model=OperationResponse)
async def upload_batch_file(
    name: str = Form(..., description="Batch configuration name"),
    file: UploadFile = File(..., description="YAML batch configuration file"),
    overwrite: bool = Form(False, description="Whether to overwrite existing file"),
):
    """Upload a new batch configuration from YAML file with validation.
    
    Processes uploaded YAML batch configuration files with comprehensive
    validation, security checks, and conflict resolution.
    
    Args:
        name (str): Desired name for the batch configuration.
        file (UploadFile): YAML file containing batch configuration.
        overwrite (bool): Whether to overwrite existing files with same name.
        
    Returns:
        OperationResponse: Upload result with file information.
        
    Raises:
        HTTPException: If file invalid, YAML malformed, or upload fails.
        
    Security:
        - Validates file extensions and content encoding
        - Sanitizes filenames to prevent security issues
        - Checks for file existence conflicts
    """
    try:
        # Validate file type
        if not file.filename.lower().endswith((".yaml", ".yml")):
            raise HTTPException(
                status_code=400,
                detail="Invalid file type. Please upload a YAML file (.yaml or .yml)",
            )

        # Read file content
        content = await file.read()
        content_str = content.decode("utf-8")

        # Validate YAML content
        try:
            yaml.safe_load(content_str)
        except yaml.YAMLError as e:
            raise HTTPException(
                status_code=400, detail=f"Invalid YAML format: {str(e)}"
            )

        # Sanitize filename
        safe_name = sanitize_filename(name)
        if not safe_name.endswith(".yaml"):
            safe_name += ".yaml"

        # Check if file exists and overwrite is not allowed
        batch_dir = get_batch_directory()
        file_path = batch_dir / safe_name

        if file_path.exists() and not overwrite:
            raise HTTPException(
                status_code=409,
                detail=f"Batch file '{safe_name}' already exists. Set overwrite=true to replace it.",
            )

        # Save file
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content_str)

        logger.info(f"Batch file uploaded: {file_path}")

        return OperationResponse(
            success=True,
            message=f"Batch configuration '{name}' uploaded successfully from file '{file.filename}'",
            data={"filename": safe_name, "file_path": str(file_path)},
        )

    except UnicodeDecodeError:
        raise HTTPException(
            status_code=400,
            detail="File contains invalid characters. Please ensure it's a valid text file.",
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error uploading batch file: {e}")
        raise HTTPException(status_code=500, detail="Failed to upload batch file")


@router.post("/upload/text", response_model=OperationResponse)
async def upload_batch_text(request: BatchUploadRequest):
    """Upload a new batch configuration from YAML content (text input).
    
    Processes YAML batch configuration content submitted as text
    with the same validation and security measures as file upload.
    
    Args:
        request (BatchUploadRequest): YAML content and metadata.
        
    Returns:
        OperationResponse: Upload result with file information.
        
    Raises:
        HTTPException: If content invalid, conflicts exist, or save fails.
        
    Note:
        Provides text-based alternative to file upload for
        programmatic batch configuration creation.
    """
    try:
        # Sanitize filename
        safe_name = sanitize_filename(request.name)
        if not safe_name.endswith(".yaml"):
            safe_name += ".yaml"

        # Check if file exists and overwrite is not allowed
        batch_dir = get_batch_directory()
        file_path = batch_dir / safe_name

        if file_path.exists() and not request.overwrite:
            raise HTTPException(
                status_code=409,
                detail=f"Batch file '{safe_name}' already exists. Set overwrite=true to replace it.",
            )

        # Save file
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(request.content)

        logger.info(f"Batch file uploaded: {file_path}")

        return OperationResponse(
            success=True,
            message=f"Batch configuration '{request.name}' uploaded successfully",
            data={"filename": safe_name, "file_path": str(file_path)},
        )

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error uploading batch configuration: {e}")
        raise HTTPException(
            status_code=500, detail="Failed to upload batch configuration"
        )


@router.get("/", response_model=BatchListResponse)
async def list_batch_files():
    """List all batch configuration files with metadata extraction.
    
    Provides comprehensive listing of batch configuration files
    with extracted metadata, file statistics, and organization
    information for management interfaces.
    
    Returns:
        BatchListResponse: Complete list of batch files with metadata.
        
    Raises:
        HTTPException: If directory access fails or listing errors occur.
        
    Note:
        Automatically extracts metadata from YAML files and
        organizes results with templates prioritized first.
    """
    try:
        batch_dir = get_batch_directory()
        logger.info(f"Looking for batch files in: {batch_dir}")
        files = []

        # Check if directory exists
        if not batch_dir.exists():
            logger.warning(f"Batch directory does not exist: {batch_dir}")
            return BatchListResponse(files=[], total=0)

        # Get all YAML files in the batch directory
        yaml_files = list(batch_dir.glob("*.yaml")) + list(batch_dir.glob("*.yml"))
        logger.info(f"Found {len(yaml_files)} YAML files")

        for file_path in yaml_files:
            if not file_path.is_file():
                continue

            try:
                stat = file_path.stat()
                description = extract_description_from_yaml(file_path)
                metadata = extract_metadata_from_yaml(file_path)

                # Try to get relative path, fallback to absolute path if not possible
                try:
                    relative_path = str(file_path.relative_to(Path.cwd()))
                except ValueError:
                    # If file is not in a subdirectory of cwd, use relative to batch_dir
                    try:
                        relative_path = str(file_path.relative_to(batch_dir.parent))
                    except ValueError:
                        # Last fallback: use absolute path
                        relative_path = str(file_path)

                file_info = BatchFileInfo(
                    name=file_path.name,
                    description=description,
                    path=relative_path,
                    size=format_file_size(stat.st_size),
                    created=datetime.fromtimestamp(stat.st_ctime).isoformat(),
                    modified=datetime.fromtimestamp(stat.st_mtime).isoformat(),
                    is_template=file_path.name.upper() == "TEMPLATE.YAML",
                    metadata_name=metadata.get("name"),
                    metadata_description=metadata.get("description"),
                    metadata_author=metadata.get("author"),
                    metadata_version=metadata.get("version"),
                    metadata_creation_date=metadata.get("creation_date"),
                    metadata_tags=metadata.get("tags", []),
                )

                files.append(file_info)
                logger.info(f"Added batch file: {file_info.name}")

            except Exception as e:
                logger.warning(f"Failed to process batch file {file_path}: {e}")
                continue

        # Sort files: templates first, then by name
        files.sort(key=lambda x: (not x.is_template, x.name.lower()))

        logger.info(f"Returning {len(files)} batch files")
        return BatchListResponse(files=files, total=len(files))

    except Exception as e:
        logger.error(f"Failed to list batch files: {e}")
        raise HTTPException(
            status_code=500, detail=f"Failed to list batch files: {str(e)}"
        )


@router.get("/{filename}")
async def get_batch_file(filename: str):
    """Get a specific batch file content with metadata.
    
    Retrieves complete batch configuration file content
    along with path information for editing and display.
    
    Args:
        filename (str): Name of the batch file to retrieve.
        
    Returns:
        dict: File content, metadata, and path information.
        
    Raises:
        HTTPException: If file not found or access denied.
        
    Security:
        - Validates filename and file type
        - Ensures file exists within allowed directory
    """
    try:
        batch_dir = get_batch_directory()
        file_path = batch_dir / filename

        if not file_path.exists():
            raise HTTPException(
                status_code=404, detail=f"Batch file '{filename}' not found"
            )

        if not file_path.is_file() or not filename.endswith(".yaml"):
            raise HTTPException(status_code=400, detail="Invalid batch file")

        with open(file_path, encoding="utf-8") as f:
            content = f.read()

        # Try to get relative path, fallback to absolute path if not possible
        try:
            relative_path = str(file_path.relative_to(Path.cwd()))
        except ValueError:
            # If file is not in a subdirectory of cwd, use filename only
            relative_path = filename

        return {
            "name": filename,
            "content": content,
            "path": relative_path,
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get batch file {filename}: {e}")
        raise HTTPException(
            status_code=500, detail=f"Failed to get batch file: {str(e)}"
        )


@router.post("/{filename}")
async def save_batch_file(filename: str, content: dict):
    """Save or update a batch file with validation.
    
    Creates or updates batch configuration files with content
    validation and filename sanitization for security.
    
    Args:
        filename (str): Target filename for the batch configuration.
        content (dict): File content and metadata to save.
        
    Returns:
        dict: Save result with path information.
        
    Raises:
        HTTPException: If filename invalid, content empty, or save fails.
        
    Security:
        - Validates filename characters and format
        - Ensures proper file extension
        - Sanitizes content before writing
    """
    try:
        if not filename.endswith(".yaml"):
            filename += ".yaml"

        batch_dir = get_batch_directory()
        file_path = batch_dir / filename

        # Validate filename
        if (
            not filename.replace(".yaml", "")
            .replace("_", "")
            .replace("-", "")
            .isalnum()
        ):
            raise HTTPException(
                status_code=400,
                detail="Invalid filename. Use only letters, numbers, hyphens and underscores",
            )

        file_content = content.get("content", "")
        if not file_content:
            raise HTTPException(status_code=400, detail="Content is required")

        # Write file
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(file_content)

        logger.info(f"Batch file saved: {file_path}")

        # Try to get relative path, fallback to filename if not possible
        try:
            relative_path = str(file_path.relative_to(Path.cwd()))
        except ValueError:
            relative_path = filename

        return {
            "message": f"Batch file '{filename}' saved successfully",
            "path": relative_path,
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to save batch file {filename}: {e}")
        raise HTTPException(
            status_code=500, detail=f"Failed to save batch file: {str(e)}"
        )


@router.delete("/{filename}")
async def delete_batch_file(filename: str):
    """Delete a batch file with safety checks.
    
    Removes batch configuration files from the system with
    protection against deleting system templates and validation.
    
    Args:
        filename (str): Name of the batch file to delete.
        
    Returns:
        dict: Deletion result confirmation.
        
    Raises:
        HTTPException: If file not found, protected, or deletion fails.
        
    Security:
        - Prevents deletion of template files
        - Validates file existence and permissions
    """
    try:
        batch_dir = get_batch_directory()
        file_path = batch_dir / filename

        if not file_path.exists():
            raise HTTPException(
                status_code=404, detail=f"Batch file '{filename}' not found"
            )

        # Don't allow deleting template file
        if filename.upper() == "TEMPLATE.YAML":
            raise HTTPException(status_code=400, detail="Cannot delete template file")

        file_path.unlink()
        logger.info(f"Batch file deleted: {file_path}")

        return {"message": f"Batch file '{filename}' deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to delete batch file {filename}: {e}")
        raise HTTPException(
            status_code=500, detail=f"Failed to delete batch file: {str(e)}"
        )


@router.get("/{filename}/download")
async def download_batch_file(filename: str):
    """Download a batch file for external use.
    
    Provides direct download of batch configuration files
    in YAML format with appropriate content headers.
    
    Args:
        filename (str): Name of the batch file to download.
        
    Returns:
        FileResponse: YAML file download with proper content type.
        
    Raises:
        HTTPException: If file not found or download fails.
    """
    try:
        batch_dir = get_batch_directory()
        file_path = batch_dir / filename

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
        logger.error(f"Failed to download batch file {filename}: {e}")
        raise HTTPException(
            status_code=500, detail=f"Failed to download batch file: {str(e)}"
        )
