"""
Simple dataset generation endpoint for testing.
"""

import logging
import random
import uuid
from pathlib import Path

from fastapi import APIRouter

from ..core.models import DatasetGenerationRequest, DatasetGenerationResult

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/datasets-simple", tags=["datasets-simple"])


@router.post("/generate/synthetic", response_model=DatasetGenerationResult)
async def generate_synthetic_dataset_simple(request: DatasetGenerationRequest):
    """Generate synthetic dataset with simple implementation."""
    try:
        logger.info(
            f"Starting simple synthetic dataset generation with params: {request}"
        )
        session_id = str(uuid.uuid4())

        # Set seed if provided
        if request.seed:
            random.seed(request.seed)

        # Generate random sequences
        sequences = []
        for i in range(request.num_strings):
            sequence = "".join(
                random.choice(request.alphabet) for _ in range(request.string_length)
            )
            sequences.append(sequence)

        logger.info(f"Generated {len(sequences)} sequences")

        # Save to datasets folder
        datasets_dir = Path("datasets")
        datasets_dir.mkdir(exist_ok=True)

        filename = f"synthetic_n{request.num_strings}_L{request.string_length}_{request.alphabet}_{session_id}.fasta"
        file_path = datasets_dir / filename

        # Write FASTA file
        with open(file_path, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq_{i+1}\n{seq}\n")

        logger.info(f"Dataset saved to: {file_path}")

        # Prepare dataset info for response
        file_size = file_path.stat().st_size
        created_at = file_path.stat().st_mtime
        dataset_info = {
            "filename": filename,
            "num_sequences": len(sequences),
            "sequence_length": request.string_length,
            "alphabet": request.alphabet,
            "type": "synthetic",
            "file_size": f"{file_size} bytes",
            "created_at": created_at,  # Keep in seconds, JavaScript will convert
        }

        # Get first 5 sequences for preview
        preview_sequences = sequences[:5]

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
                "actual_distance": None,
            },
            dataset_info=dataset_info,
            sequences=preview_sequences,
        )

    except Exception as e:
        logger.error(f"Error generating simple synthetic dataset: {e}")
        return DatasetGenerationResult(
            session_id=session_id if "session_id" in locals() else "unknown",
            status="failed",
            error=str(e),
        )
