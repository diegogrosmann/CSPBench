"""
Test endpoint for debugging saved datasets.
"""

from fastapi import APIRouter

router = APIRouter(prefix="/api/test", tags=["test"])


@router.get("/saved-datasets-debug")
async def debug_saved_datasets():
    """Debug endpoint for saved datasets."""
    return {
        "status": "ok",
        "saved_datasets": [
            {
                "filename": "test1.fasta",
                "num_sequences": 10,
                "sequence_length": 20,
                "alphabet": "ACGT",
                "type": "test",
                "file_size": "100 bytes",
                "created_at": 1754053000000,
            },
            {
                "filename": "test2.fasta",
                "num_sequences": 5,
                "sequence_length": 15,
                "alphabet": "ACGT",
                "type": "test",
                "file_size": "50 bytes",
                "created_at": 1754053100000,
            },
        ],
    }
