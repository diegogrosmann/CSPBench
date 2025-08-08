"""
Dataset-related API endpoints.
"""

import asyncio
import logging
import tempfile
import uuid
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from ..core.config import web_config
from ..core.models import (
    DatasetGenerationRequest,
    DatasetGenerationResult,
    NCBIDatasetRequest,
)

logger = logging.getLogger(__name__)


def dataset_to_fasta(dataset) -> str:
    """Convert Dataset object to FASTA format string."""
    fasta_lines = []
    for i, sequence in enumerate(dataset.sequences):
        header = f">seq_{i+1}"
        fasta_lines.append(header)
        fasta_lines.append(sequence)
    return "\n".join(fasta_lines)


router = APIRouter(prefix="/api/datasets", tags=["datasets"])


@router.get("/samples")
async def get_sample_datasets():
    """Get list of sample datasets."""
    try:
        datasets_dir = Path("datasets")
        if not datasets_dir.exists():
            return {"samples": []}

        sample_files = []
        for file_path in datasets_dir.glob("*.fasta"):
            try:
                # Get basic file info
                file_size = file_path.stat().st_size

                # Read first few lines to get basic info
                with open(file_path) as f:
                    first_lines = f.read(1000)
                    num_sequences = first_lines.count(">")

                sample_files.append(
                    {
                        "name": file_path.name,
                        "size": file_size,
                        "estimated_sequences": num_sequences,
                        "path": str(file_path.relative_to(datasets_dir)),
                    }
                )
            except Exception as e:
                logger.warning(f"Error reading sample file {file_path}: {e}")

        return {"samples": sample_files}

    except Exception as e:
        logger.error(f"Error getting sample datasets: {e}")
        raise HTTPException(
            status_code=500, detail="Failed to retrieve sample datasets"
        )


@router.post("/generate/synthetic", response_model=DatasetGenerationResult)
async def generate_synthetic_dataset(request: DatasetGenerationRequest):
    """Generate synthetic dataset."""
    try:
        logger.info(f"Starting synthetic dataset generation with params: {request}")
        session_id = str(uuid.uuid4())

        # Create session manager if not available
        session_manager = web_config.get_session_manager()
        if not session_manager:
            raise HTTPException(status_code=503, detail="Session manager not available")

        # Create session directory
        session_dir = Path(tempfile.mkdtemp(prefix=f"dataset_gen_{session_id}_"))
        logger.info(f"Created session directory: {session_dir}")

        # Generate dataset using the domain service
        from src.application.services.dataset_generation_service import (
            DatasetGenerationService,
        )
        from src.infrastructure.persistence.dataset_repository import (
            FileDatasetRepository,
        )

        # Create temporary repository for this session
        temp_repo = FileDatasetRepository(str(session_dir))
        dataset_service = DatasetGenerationService(temp_repo)
        logger.info("Created dataset service")

        # Prepare parameters
        params = {
            "method": request.generation_method,
            "n": request.num_strings,
            "length": request.string_length,
            "alphabet": request.alphabet,
            "seed": request.seed,
        }
        if request.max_distance:
            params["max_distance"] = request.max_distance

        # Generate dataset
        result = await asyncio.get_event_loop().run_in_executor(
            None, dataset_service.generate_synthetic_dataset, params
        )

        # Debug: Check result structure
        logger.info(f"Result type: {type(result)}")
        logger.info(f"Result attributes: {dir(result)}")
        if hasattr(result, "sequences"):
            logger.info(f"Number of sequences: {len(result.sequences)}")
            logger.info(
                f"First sequence: {result.sequences[0] if result.sequences else 'None'}"
            )
        else:
            logger.info("Result does not have 'sequences' attribute")

        # Save to file
        filename = f"synthetic_n{request.num_strings}_L{request.string_length}_{request.alphabet}_{session_id}.fasta"
        file_path = session_dir / filename

        with open(file_path, "w") as f:
            f.write(dataset_to_fasta(result))

        # Also save to datasets folder for persistence
        datasets_dir = Path("datasets")
        datasets_dir.mkdir(exist_ok=True)
        persistent_file_path = datasets_dir / filename

        with open(persistent_file_path, "w") as f:
            f.write(dataset_to_fasta(result))

        logger.info(f"Dataset saved to: {persistent_file_path}")

        # Prepare dataset info for response
        file_size = persistent_file_path.stat().st_size
        dataset_info = {
            "filename": filename,
            "num_sequences": len(result.sequences),
            "sequence_length": request.string_length,
            "alphabet": request.alphabet,
            "type": "synthetic",
            "file_size": f"{file_size} bytes",
            "created_at": persistent_file_path.stat().st_mtime,
        }

        # Get first 5 sequences for preview
        preview_sequences = result.sequences[:5] if result.sequences else []

        logger.info(f"Returning dataset info: {dataset_info}")
        logger.info(f"Preview sequences: {preview_sequences}")

        return DatasetGenerationResult(
            session_id=session_id,
            status="completed",
            filename=filename,
            download_url=f"/api/datasets/download/{session_id}/{filename}",
            metadata={
                "num_strings": request.num_strings,
                "string_length": request.string_length,
                "alphabet": request.alphabet,
                "generation_method": request.generation_method,
                "actual_distance": getattr(result, "actual_distance", None),
            },
            dataset_info=dataset_info,
            sequences=preview_sequences,
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error generating synthetic dataset: {e}")
        return DatasetGenerationResult(
            session_id=session_id if "session_id" in locals() else "unknown",
            status="failed",
            error=str(e),
        )


@router.post("/generate/ncbi", response_model=DatasetGenerationResult)
async def generate_ncbi_dataset(request: NCBIDatasetRequest):
    """Generate dataset from NCBI."""
    try:
        logger.info(f"Starting NCBI dataset generation with params: {request}")
        session_id = str(uuid.uuid4())

        # Create session manager if not available
        session_manager = web_config.get_session_manager()
        if not session_manager:
            raise HTTPException(status_code=503, detail="Session manager not available")

        # Create session directory
        session_dir = Path(tempfile.mkdtemp(prefix=f"ncbi_dataset_{session_id}_"))
        logger.info(f"Created session directory: {session_dir}")

        # Import NCBI repository
        from src.infrastructure.persistence.entrez_dataset_repository import (
            NCBIEntrezDatasetRepository,
        )

        # Create NCBI repository
        try:
            ncbi_repo = NCBIEntrezDatasetRepository(email=request.email)
        except Exception as e:
            logger.error(f"Failed to initialize NCBI repository: {e}")
            raise HTTPException(
                status_code=400, detail=f"NCBI initialization failed: {str(e)}"
            )

        # Fetch dataset from NCBI
        try:
            sequences, metadata = await asyncio.get_event_loop().run_in_executor(
                None,
                ncbi_repo.fetch_dataset,
                request.query,
                request.sequence_type,  # Changed from database to sequence_type
                request.max_sequences,
            )

            if not sequences:
                raise HTTPException(
                    status_code=404, detail="No sequences found for the given query"
                )

            # Note: Length filtering removed as it's not in the current model
            # This can be added back if min_length/max_length are added to NCBIDatasetRequest

            # Create Dataset object
            from src.domain.dataset import Dataset

            dataset = Dataset(sequences=sequences)

            # Generate filename
            filename = f"ncbi_{request.sequence_type}_{len(sequences)}seq_{session_id[:8]}.fasta"

            # Convert to FASTA and save
            fasta_content = dataset_to_fasta(dataset)
            file_path = session_dir / filename

            with open(file_path, "w") as f:
                f.write(fasta_content)

            # Also save to datasets directory for persistence
            datasets_dir = Path("datasets")
            datasets_dir.mkdir(exist_ok=True)
            persistent_path = datasets_dir / filename

            with open(persistent_path, "w") as f:
                f.write(fasta_content)

            logger.info(f"NCBI dataset generated successfully: {filename}")

            # Prepare dataset info
            dataset_info = {
                "num_sequences": len(sequences),
                "sequence_length": len(sequences[0]) if sequences else 0,
                "alphabet": "".join(sorted(set("".join(sequences)))),
                "source": "NCBI",
                "query": request.query,
                "sequence_type": request.sequence_type,  # Changed from database to sequence_type
            }

            return DatasetGenerationResult(
                session_id=session_id,
                status="success",
                filename=filename,
                download_url=f"/datasets/{filename}",
                metadata=metadata,
                dataset_info=dataset_info,
                sequences=sequences,  # Changed from len(sequences) to sequences
            )

        except Exception as e:
            logger.error(f"NCBI fetch failed: {e}")
            error_msg = str(e)
            if "No sequences found" in error_msg:
                raise HTTPException(status_code=404, detail=error_msg)
            elif "Biopython" in error_msg:
                raise HTTPException(
                    status_code=500, detail="NCBI service dependencies not available"
                )
            else:
                raise HTTPException(
                    status_code=500, detail=f"NCBI fetch failed: {error_msg}"
                )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error generating NCBI dataset: {e}")
        return DatasetGenerationResult(
            session_id=session_id if "session_id" in locals() else "unknown",
            status="failed",
            error=str(e),
        )


@router.get("/download/{session_id}/{filename}")
async def download_dataset(session_id: str, filename: str):
    """Download generated dataset."""
    try:
        # Find the session directory (this would need proper session management)
        import tempfile

        temp_dir = Path(tempfile.gettempdir())

        # Look for session directory
        session_dirs = list(temp_dir.glob(f"*{session_id}*"))
        if not session_dirs:
            raise HTTPException(status_code=404, detail="Session not found")

        session_dir = session_dirs[0]
        file_path = session_dir / filename

        if not file_path.exists():
            raise HTTPException(status_code=404, detail="File not found")

        return FileResponse(
            path=file_path, filename=filename, media_type="application/octet-stream"
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error downloading dataset: {e}")
        raise HTTPException(status_code=500, detail="Failed to download dataset")


@router.get("/saved")
async def get_saved_datasets():
    """Get list of saved datasets."""
    try:
        datasets_dir = Path("datasets")
        if not datasets_dir.exists():
            return {"saved_datasets": []}

        saved_datasets = []
        for file_path in datasets_dir.glob("*.fasta"):
            try:
                # Get file stats
                stat = file_path.stat()
                file_size = stat.st_size
                created_at = stat.st_mtime

                # Read first few lines to get basic info
                with open(file_path) as f:
                    lines = f.readlines()

                # Count sequences
                num_sequences = sum(1 for line in lines if line.startswith(">"))

                # Get first sequence length (if any)
                sequence_length = 0
                for line in lines:
                    if not line.startswith(">") and line.strip():
                        sequence_length = len(line.strip())
                        break

                # Determine alphabet from first sequence
                alphabet = "Unknown"
                first_seq = ""
                for line in lines:
                    if not line.startswith(">") and line.strip():
                        first_seq = line.strip()
                        break

                if first_seq:
                    unique_chars = set(first_seq.upper())
                    if unique_chars.issubset(set("ACGT")):
                        alphabet = "ACGT"
                    elif unique_chars.issubset(set("ACDEFGHIKLMNPQRSTVWY")):
                        alphabet = "Protein"
                    else:
                        alphabet = "".join(sorted(unique_chars))

                dataset_info = {
                    "filename": file_path.name,
                    "num_sequences": num_sequences,
                    "sequence_length": sequence_length,
                    "alphabet": alphabet,
                    "type": "fasta",
                    "file_size": f"{file_size} bytes",
                    "created_at": created_at,  # Keep in seconds, JavaScript will convert
                }

                saved_datasets.append(dataset_info)

            except Exception as e:
                logger.warning(f"Error processing file {file_path}: {e}")
                continue

        # Sort by creation time (newest first)
        saved_datasets.sort(key=lambda x: x["created_at"], reverse=True)

        return {"saved_datasets": saved_datasets}

    except Exception as e:
        logger.error(f"Error getting saved datasets: {e}")
        raise HTTPException(status_code=500, detail="Failed to retrieve saved datasets")


@router.delete("/saved/{filename}")
async def delete_saved_dataset(filename: str):
    """Delete a saved dataset file."""
    try:
        # Validate filename for security
        if ".." in filename or "/" in filename or "\\" in filename:
            raise HTTPException(status_code=400, detail="Invalid filename")

        datasets_dir = Path("datasets")
        file_path = datasets_dir / filename

        if not file_path.exists():
            raise HTTPException(status_code=404, detail="Dataset not found")

        # Delete the file
        file_path.unlink()

        logger.info(f"Deleted dataset: {filename}")
        return {"message": f"Dataset {filename} deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error deleting dataset {filename}: {e}")
        raise HTTPException(status_code=500, detail="Failed to delete dataset")


@router.post("/upload")
async def upload_dataset(file: UploadFile = File(...), name: str = Form("")):
    """Upload a FASTA dataset file."""
    try:
        # Validate file type
        if not file.filename or not file.filename.lower().endswith((".fasta", ".fa")):
            raise HTTPException(
                status_code=400, detail="File must be a FASTA file (.fasta or .fa)"
            )

        # Use provided name or original filename
        filename = name.strip() if name else file.filename

        # Ensure .fasta extension
        if not filename.lower().endswith(".fasta"):
            filename = filename.rsplit(".", 1)[0] + ".fasta"

        datasets_dir = Path("datasets")
        datasets_dir.mkdir(exist_ok=True)

        file_path = datasets_dir / filename

        # Check if file already exists
        if file_path.exists():
            raise HTTPException(
                status_code=409, detail=f"Dataset with name '{filename}' already exists"
            )

        # Read and validate FASTA content
        content = await file.read()
        content_str = content.decode("utf-8")

        # Basic FASTA validation
        if not content_str.strip().startswith(">"):
            raise HTTPException(
                status_code=400, detail="Invalid FASTA format: file must start with '>'"
            )

        # Count sequences
        sequence_count = content_str.count(">")
        if sequence_count == 0:
            raise HTTPException(
                status_code=400, detail="No sequences found in FASTA file"
            )

        # Write file
        with open(file_path, "w") as f:
            f.write(content_str)

        logger.info(f"Uploaded dataset: {filename} with {sequence_count} sequences")

        return {
            "message": f"Dataset '{filename}' uploaded successfully",
            "filename": filename,
            "sequences": sequence_count,
            "size": len(content),
        }

    except HTTPException:
        raise
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="File must be a valid text file")
    except Exception as e:
        logger.error(f"Error uploading dataset: {e}")
        raise HTTPException(status_code=500, detail="Failed to upload dataset")
