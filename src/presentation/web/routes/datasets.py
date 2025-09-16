"""
Dataset CRUD API endpoints for CSPBench Web Interface.

This module provides REST API endpoints for dataset management including:
- CRUD operations (Create, Read, Update, Delete)
- Synthetic dataset generation
- NCBI dataset download
- Dataset preview and statistics
"""

import logging
import os
import uuid
from datetime import datetime

from fastapi import APIRouter, BackgroundTasks, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from src.application.services.dataset_generator import SyntheticDatasetGenerator
from src.domain.config import (
    EntrezDatasetConfig,
    SyntheticDatasetConfig,
)
from src.domain.dataset import Dataset
from src.domain.errors import DatasetNotFoundError
from src.domain.status import BaseStatus
from src.infrastructure.external.dataset_entrez import EntrezDatasetDownloader
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository

from ..core.dataset_models import (
    DatasetGenerationStatus,
    DatasetInfo,
    DatasetListResponse,
    DatasetPreview,
    DatasetType,
    DatasetUpdateRequest,
    DatasetUploadRequest,
    NCBIDatasetRequest,
    OperationResponse,
    SyntheticDatasetRequest,
)

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/datasets", tags=["datasets"])

# In-memory storage for generation status (in production, use Redis or database)
_generation_status = {}


def _dataset_to_info(
    dataset: Dataset,
    dataset_id: str,
    file_path: str = None,
    generation_params: dict = None,
) -> DatasetInfo:
    """Convert Dataset domain object to DatasetInfo response model.

    Transforms internal Dataset objects into API response models with
    comprehensive metadata, statistics, and type information.

    Args:
        dataset (Dataset): Domain dataset object to convert.
        dataset_id (str): Unique identifier for the dataset.
        file_path (str, optional): Filesystem path if dataset is persisted.
        generation_params (dict, optional): Parameters used for dataset generation.

    Returns:
        DatasetInfo: Complete dataset information for API responses.

    Note:
        Automatically determines dataset type based on generation parameters
        and calculates human-readable file sizes when applicable.
    """
    stats = dataset.get_statistics()

    # Convert file_path to string if it's a Path object
    file_path_str = str(file_path) if file_path else None

    # Calculate file size if file exists
    file_size = None
    if file_path_str and os.path.exists(file_path_str):
        size_bytes = os.path.getsize(file_path_str)
        if size_bytes < 1024:
            file_size = f"{size_bytes} B"
        elif size_bytes < 1024 * 1024:
            file_size = f"{size_bytes / 1024:.1f} KB"
        else:
            file_size = f"{size_bytes / (1024 * 1024):.1f} MB"

    # Determine dataset type
    dataset_type = DatasetType.FILE
    if generation_params:
        if "mode" in generation_params or "method" in generation_params:
            dataset_type = DatasetType.SYNTHETIC
        elif "query" in generation_params:
            dataset_type = DatasetType.NCBI

    return DatasetInfo(
        id=dataset_id,
        name=dataset.name or dataset_id,
        type=dataset_type,
        size=stats["size"],
        min_length=stats["min_length"],
        max_length=stats["max_length"],
        average_length=stats["average_length"],
        alphabet=stats["alphabet"],
        alphabet_size=stats["alphabet_size"],
        diversity=stats["diversity"],
        uniform_lengths=stats["uniform_lengths"],
        file_path=file_path_str,
        file_size=file_size,
        created_at=datetime.now(),
        generation_params=generation_params,
    )


@router.get("/", response_model=DatasetListResponse)
async def list_datasets():
    """List all available datasets with metadata and statistics.

    Provides a comprehensive listing of all datasets in the system
    with essential metadata for browsing and selection purposes.

    Returns:
        DatasetListResponse: Paginated list of datasets with total count.

    Raises:
        HTTPException: If dataset listing fails due to persistence errors.

    Note:
        Continues processing even if individual dataset loading fails,
        logging warnings for problematic datasets.
    """
    try:
        dataset_names = FileDatasetRepository.list_available()
        datasets = []

        for name in dataset_names:
            try:
                dataset, params = FileDatasetRepository.load(name)
                file_path = params.get("file_path")
                dataset_info = _dataset_to_info(dataset, name, file_path, params)
                datasets.append(dataset_info)
            except Exception as e:
                logger.warning(f"Error loading dataset {name}: {e}")
                continue

        return DatasetListResponse(datasets=datasets, total=len(datasets))

    except Exception as e:
        logger.error(f"Error listing datasets: {e}")
        raise HTTPException(status_code=500, detail="Failed to list datasets")


@router.get("/{dataset_id}", response_model=DatasetInfo)
async def get_dataset(dataset_id: str):
    """Get detailed information about a specific dataset.

    Retrieves comprehensive metadata, statistics, and configuration
    information for a single dataset by its unique identifier.

    Args:
        dataset_id (str): Unique identifier of the dataset to retrieve.

    Returns:
        DatasetInfo: Complete dataset information including statistics.

    Raises:
        HTTPException: If dataset not found or retrieval fails.
    """
    try :
        if not FileDatasetRepository.exists(dataset_id):
            raise HTTPException(
                status_code=404, detail=f"Dataset '{dataset_id}' not found"
            )

        dataset, params = FileDatasetRepository.load(dataset_id)
        file_path = params.get("file_path")

        return _dataset_to_info(dataset, dataset_id, file_path, params)

    except DatasetNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found")
    except Exception as e:
        logger.error(f"Error getting dataset {dataset_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to retrieve dataset")

@router.get("/{dataset_id}/preview", response_model=DatasetPreview)
async def get_dataset_preview(dataset_id: str, max_sequences: int = 5):
    """Get dataset preview with sample sequences for inspection.

    Provides a preview of dataset contents including metadata
    and a configurable number of sample sequences for user inspection
    before processing or analysis.

    Args:
        dataset_id (str): Unique identifier of the dataset to preview.
        max_sequences (int): Maximum number of sequences to include in preview.

    Returns:
        DatasetPreview: Dataset information with sample sequence data.

    Raises:
        HTTPException: If dataset not found or preview generation fails.
    """
    try :
        if not FileDatasetRepository.exists(dataset_id):
            raise HTTPException(
                status_code=404, detail=f"Dataset '{dataset_id}' not found"
            )

        dataset, params = FileDatasetRepository.load(dataset_id)
        file_path = params.get("file_path")

        dataset_info = _dataset_to_info(dataset, dataset_id, file_path, params)

        # Get sample sequences
        sample_sequences = dataset.sequences[:max_sequences]

        return DatasetPreview(
            info=dataset_info,
            sample_sequences=sample_sequences,
            total_sequences=len(dataset.sequences),
        )

    except DatasetNotFoundError:
        raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found")
    except Exception as e:
        logger.error(f"Error getting dataset preview {dataset_id}: {e}")
        raise HTTPException(
            status_code=500, detail="Failed to retrieve dataset preview"
        )


@router.post("/upload", response_model=OperationResponse)
async def upload_dataset_file(
    name: str = Form(..., description="Dataset name"),
    file: UploadFile = File(..., description="FASTA file"),
):
    """Upload a new dataset from FASTA file with validation.

    Processes uploaded FASTA files, validates format, parses sequences,
    and creates a new dataset entry in the system with proper metadata.

    Args:
        name (str): Desired name for the new dataset.
        file (UploadFile): FASTA file upload containing sequence data.

    Returns:
        OperationResponse: Upload result with dataset ID and file path.

    Raises:
        HTTPException: If file format invalid, parsing fails, or upload errors.

    Security:
        - Validates file extensions for allowed types
        - Checks file content encoding and format
        - Sanitizes input data before processing
    """
    try:
        # Validate file type
        if not file.filename.lower().endswith((".fasta", ".fa", ".txt")):
            raise HTTPException(
                status_code=400,
                detail="Invalid file type. Please upload a FASTA file (.fasta, .fa, or .txt)",
            )

        # Read file content
        content = await file.read()
        content_str = content.decode("utf-8")

        # Parse FASTA content
        sequences = []
        current_seq = []

        for line in content_str.strip().split("\n"):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)

        if current_seq:
            sequences.append("".join(current_seq))

        if not sequences:
            raise HTTPException(
                status_code=400, detail="No sequences found in FASTA file"
            )

        # Create dataset
        dataset = Dataset(name=name, sequences=sequences)

        # Generate unique ID
        dataset_id = f"{name}_{uuid.uuid4().hex[:8]}"

        # Save dataset
        file_path = FileDatasetRepository.save(dataset, dataset_id)

        return OperationResponse(
            success=True,
            message=f"Dataset '{name}' uploaded successfully from file '{file.filename}'",
            data={"dataset_id": dataset_id, "file_path": file_path},
        )

    except UnicodeDecodeError:
        raise HTTPException(
            status_code=400,
            detail="File contains invalid characters. Please ensure it's a valid text file.",
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error uploading dataset file: {e}")
        raise HTTPException(status_code=500, detail="Failed to upload dataset file")

@router.post("/upload/text", response_model=OperationResponse)
async def upload_dataset_text(request: DatasetUploadRequest):
    """Upload a new dataset from FASTA content (text input).

    Processes FASTA content submitted as text, validates format,
    and creates a new dataset entry with proper parsing and validation.

    Args:
        request (DatasetUploadRequest): FASTA content and metadata.

    Returns:
        OperationResponse: Upload result with dataset ID and storage info.

    Raises:
        HTTPException: If content format invalid or processing fails.

    Note:
        Provides the same functionality as file upload but accepts
        text input directly from web forms or API clients.
    """
    try :
        # Parse FASTA content
        sequences = []
        current_seq = []

        for line in request.content.strip().split("\n"):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)

        if current_seq:
            sequences.append("".join(current_seq))

        if not sequences:
            raise HTTPException(
                status_code=400, detail="No sequences found in FASTA content"
            )

        # Create dataset
        dataset = Dataset(name=request.name, sequences=sequences)

        # Generate unique ID
        dataset_id = f"{request.name}_{uuid.uuid4().hex[:8]}"

        # Save dataset
        file_path = FileDatasetRepository.save(dataset, dataset_id)

        return OperationResponse(
            success=True,
            message=f"Dataset '{request.name}' uploaded successfully",
            data={"dataset_id": dataset_id, "file_path": file_path},
        )

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error uploading dataset: {e}")
        raise HTTPException(status_code=500, detail="Failed to upload dataset")


@router.post("/generate/synthetic", response_model=OperationResponse)
async def generate_synthetic_dataset(
    request: SyntheticDatasetRequest, background_tasks: BackgroundTasks
):
    """Generate a synthetic dataset using algorithmic methods.

    Initiates background generation of synthetic biological sequences
    using various algorithms and parameters with progress tracking.

    Args:
        request (SyntheticDatasetRequest): Generation parameters and configuration.
        background_tasks (BackgroundTasks): FastAPI background task manager.

    Returns:
        OperationResponse: Generation initiation result with task ID.

    Raises:
        HTTPException: If generation setup fails or parameters invalid.

    Note:
        Generation runs in background with status tracking via task ID.
        Use generation status endpoint to monitor progress.
    """
    try :
        # Generate unique task ID
        task_id = str(uuid.uuid4())

        # Initialize status
        _generation_status[task_id] = DatasetGenerationStatus(
            status=BaseStatus.RUNNING.value,
            progress=0,
            message="Starting synthetic dataset generation...",
        )

        # Start background generation
        background_tasks.add_task(_generate_synthetic_background, task_id, request)

        return OperationResponse(
            success=True,
            message="Synthetic dataset generation started",
            data={"task_id": task_id},
        )

    except Exception as e :
        logger.error(f"Error starting synthetic generation: {e}")
        raise HTTPException(
            status_code=500, detail="Failed to start dataset generation"
        )

@router.post("/generate/ncbi", response_model=OperationResponse)
async def generate_ncbi_dataset(
    request: NCBIDatasetRequest, background_tasks: BackgroundTasks
):
    """Download dataset from NCBI databases with filtering options.

    Initiates background download of biological sequences from NCBI
    databases with query filtering and processing parameters.

    Args:
        request (NCBIDatasetRequest): NCBI query parameters and filters.
        background_tasks (BackgroundTasks): FastAPI background task manager.

    Returns:
        OperationResponse: Download initiation result with task ID.

    Raises:
        HTTPException: If download setup fails or parameters invalid.

    Note:
        Download runs in background with progress tracking.
        Large queries may take significant time to complete.
    """
    try :
        # Generate unique task ID
        task_id = str(uuid.uuid4())

        # Initialize status
        _generation_status[task_id] = DatasetGenerationStatus(
            status=BaseStatus.RUNNING.value,
            progress=0,
            message="Starting NCBI dataset download...",
        )

        # Start background download
        background_tasks.add_task(_generate_ncbi_background, task_id, request)

        return OperationResponse(
            success=True,
            message="NCBI dataset download started",
            data={"task_id": task_id},
        )

    except Exception as e :
        logger.error(f"Error starting NCBI download: {e}")
        raise HTTPException(status_code=500, detail="Failed to start dataset download")

@router.get("/generation/{task_id}/status", response_model=DatasetGenerationStatus)
async def get_generation_status(task_id: str):
    """Get dataset generation status and progress information.

    Provides real-time status updates for background dataset generation
    or download tasks with progress indication and error reporting.

    Args:
        task_id (str): Unique identifier of the generation task.

    Returns:
        DatasetGenerationStatus: Current status, progress, and result information.

    Raises:
        HTTPException: If task ID not found or status retrieval fails.

    Note:
        Completed or failed tasks are removed after being accessed
        to prevent memory leaks in the status storage.
    """
    if task_id not in _generation_status:
        raise HTTPException(status_code=404, detail="Generation task not found")

    status = _generation_status[task_id]

    # Remove completed or failed tasks after being accessed to prevent loop
    if status.status in [BaseStatus.COMPLETED.value, BaseStatus.FAILED.value]:
        # Store the status before removing
        final_status = _generation_status.pop(task_id)
        return final_status

    return status


@router.put("/{dataset_id}", response_model=OperationResponse)
async def update_dataset(dataset_id: str, request: DatasetUpdateRequest):
    """Update dataset metadata and properties.

    Modifies mutable dataset properties such as name and description
    while preserving sequence data and statistics.

    Args:
        dataset_id (str): Unique identifier of dataset to update.
        request (DatasetUpdateRequest): Updated metadata values.

    Returns:
        OperationResponse: Update result with success confirmation.

    Raises:
        HTTPException: If dataset not found or update fails.
    """
    try :
        if not FileDatasetRepository.exists(dataset_id):
            raise HTTPException(
                status_code=404, detail=f"Dataset '{dataset_id}' not found"
            )

        # Load current dataset
        dataset, params = FileDatasetRepository.load(dataset_id)

        # Update name if provided
        if request.name:
            dataset.name = request.name

            # Save with new name (keeping same ID for consistency)
            FileDatasetRepository.save(dataset, dataset_id)

        return OperationResponse(
            success=True, message=f"Dataset '{dataset_id}' updated successfully"
        )

    except DatasetNotFoundError :
        raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found")
    except Exception as e :
        logger.error(f"Error updating dataset {dataset_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to update dataset")
    
@router.delete("/{dataset_id}", response_model=OperationResponse)
async def delete_dataset(dataset_id: str):
    """Delete a dataset and its associated files.

    Removes dataset from the system including metadata and
    associated files with proper cleanup and validation.

    Args:
        dataset_id (str): Unique identifier of dataset to delete.

    Returns:
        OperationResponse: Deletion result with success confirmation.

    Raises:
        HTTPException: If dataset not found or deletion fails.

    Warning:
        This operation is irreversible. Ensure proper backup
        procedures if dataset preservation is required.
    """
    try :
        if not FileDatasetRepository.exists(dataset_id):
            raise HTTPException(
                status_code=404, detail=f"Dataset '{dataset_id}' not found"
            )

        success = FileDatasetRepository.delete(dataset_id)

        if success:
            return OperationResponse(
                success=True, message=f"Dataset '{dataset_id}' deleted successfully"
            )
        else:
            raise HTTPException(status_code=500, detail="Failed to delete dataset")

    except Exception as e :
        logger.error(f"Error deleting dataset {dataset_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to delete dataset")
    
@router.get("/{dataset_id}/download", response_class=FileResponse)
async def download_dataset(dataset_id: str):
    """Download dataset file in FASTA format.

    Provides direct file download of dataset in standard FASTA
    format for use in external tools and applications.

    Args:
        dataset_id (str): Unique identifier of dataset to download.

    Returns:
        FileResponse: FASTA file download with appropriate headers.

    Raises:
        HTTPException: If dataset or file not found.
    """
    try :
        if not FileDatasetRepository.exists(dataset_id):
            raise HTTPException(
                status_code=404, detail=f"Dataset '{dataset_id}' not found"
            )

        # Get file path
        _, params = FileDatasetRepository.load(dataset_id)
        file_path = params.get("file_path")

        if not file_path or not os.path.exists(file_path):
            raise HTTPException(status_code=404, detail="Dataset file not found")

        return FileResponse(
            path=file_path, filename=f"{dataset_id}.fasta", media_type="text/plain"
        )

    except Exception as e :
        logger.error(f"Error downloading dataset {dataset_id}: {e}")
        raise HTTPException(status_code=500, detail="Failed to download dataset")
    
async def _generate_synthetic_background(
    task_id: str, request: SyntheticDatasetRequest
):
    """Background task for synthetic dataset generation with progress tracking.

    Executes synthetic dataset generation in background with progress
    updates and error handling for various generation methods.

    Args:
        task_id (str): Unique identifier for progress tracking.
        request (SyntheticDatasetRequest): Generation configuration and parameters.

    Note:
        Updates global status dictionary with progress information
        and handles all generation methods with appropriate parameter mapping.
    """
    try :
        # Update status
        _generation_status[task_id].progress = 25
        _generation_status[task_id].message = "Configuring generation parameters..."

        # Create configuration using the correct field names
        config = SyntheticDatasetConfig(
            id=f"synthetic_{uuid.uuid4().hex[:8]}",
            name=request.name,
            mode=request.method.value,
            n=request.n,
            alphabet=request.alphabet,
            seed=request.seed,
        )

        # Add method-specific parameters
        if request.method.value == "random" or request.method.value == "clustered":
            config.L = request.length or 50

        if request.method.value == "noise" and request.base_sequence:
            config.parameters_mode["base_sequence"] = request.base_sequence

        if request.method.value == "noise" and request.noise_rate is not None:
            config.parameters_mode["noise_rate"] = request.noise_rate

        if request.method.value == "clustered":
            if request.num_clusters is not None:
                config.parameters_mode["num_clusters"] = request.num_clusters
            if request.cluster_distance is not None:
                config.parameters_mode["cluster_distance"] = request.cluster_distance

        if request.method.value == "mutations":
            if request.base_sequence:
                config.parameters_mode["base_sequence"] = request.base_sequence
            if request.mutation_rate is not None:
                config.parameters_mode["mutation_rate"] = request.mutation_rate
            if request.mutation_types:
                config.parameters_mode["mutation_types"] = request.mutation_types

        # Update status
        _generation_status[task_id].progress = 50
        _generation_status[task_id].message = "Generating sequences..."

        # Generate dataset
        dataset, params = SyntheticDatasetGenerator.generate_from_config(config)

        # Update status
        _generation_status[task_id].progress = 75
        _generation_status[task_id].message = "Saving dataset..."

        # Generate unique ID and save
        dataset_id = f"{request.name}_{uuid.uuid4().hex[:8]}"
        file_path = FileDatasetRepository.save(dataset, dataset_id)

        # Update status - completed
        _generation_status[task_id] = DatasetGenerationStatus(
            status=BaseStatus.COMPLETED.value,
            progress=100,
            message="Synthetic dataset generated successfully",
            dataset_id=dataset_id,
        )

    except Exception as e :
        logger.error(f"Error in synthetic generation task {task_id}: {e}")
        _generation_status[task_id] = DatasetGenerationStatus(
            status=BaseStatus.FAILED.value,
            progress=0,
            message="Generation failed",
            error=str(e),
        )
    
async def _generate_ncbi_background(task_id: str, request: NCBIDatasetRequest):
    """Background task for NCBI dataset download with progress tracking.

    Executes NCBI database queries and sequence downloads in background
    with progress updates and comprehensive error handling.

    Args:
        task_id (str): Unique identifier for progress tracking.
        request (NCBIDatasetRequest): NCBI query configuration and filters.

    Note:
        Handles NCBI API communication, sequence filtering, and
        dataset creation with appropriate error recovery.
    """
    try :
        # Update status
        _generation_status[task_id].progress = 25
        _generation_status[task_id].message = "Connecting to NCBI..."

        # Create configuration
        config = EntrezDatasetConfig(
            id=f"ncbi_{uuid.uuid4().hex[:8]}",
            name=request.name,
            query=request.query,
            db=request.db,
            retmax=request.max_sequences or 100,
            min_length=request.min_length,
            max_length=request.max_length,
            uniform_policy=request.uniform_policy,
        )

        # Update status
        _generation_status[task_id].progress = 50
        _generation_status[task_id].message = "Downloading sequences from NCBI..."

        # Download dataset
        dataset, params = EntrezDatasetDownloader.download(config)

        # Update status
        _generation_status[task_id].progress = 75
        _generation_status[task_id].message = "Saving dataset..."

        # Generate unique ID and save
        dataset_id = f"{request.name}_{uuid.uuid4().hex[:8]}"
        file_path = FileDatasetRepository.save(dataset, dataset_id)

        # Update status - completed
        _generation_status[task_id] = DatasetGenerationStatus(
            status=BaseStatus.COMPLETED.value,
            progress=100,
            message="NCBI dataset downloaded successfully",
            dataset_id=dataset_id,
        )

    except Exception as e :
        logger.error(f"Error in NCBI download task {task_id}: {e}")
        _generation_status[task_id] = DatasetGenerationStatus(
            status=BaseStatus.FAILED.value,
            progress=0,
            message="Download failed",
            error=str(e),
        )
    
