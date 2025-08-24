"""
File-based Dataset Repository

Implements DatasetRepository using file system.
Respects DATASET_DIRECTORY from environment for base path.
"""

import os
from pathlib import Path
from typing import List, Any

from src.domain import Dataset
from src.domain.errors import DatasetNotFoundError, DatasetValidationError


class FileDatasetRepository:
    """FASTA file-based dataset repository."""

    @staticmethod
    def _get_base_path(base_path: str | None = None) -> Path:
        """Get base path, using environment variable as default."""
        if base_path is None:
            base_path = os.getenv("DATASET_DIRECTORY", "./datasets")
        path = Path(base_path)
        path.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def save(dataset: Dataset, name: str, base_path: str | None = None) -> str:
        """Save dataset to FASTA file."""
        base = FileDatasetRepository._get_base_path(base_path)
        file_path = base / f"{name}.fasta"

        with open(file_path, "w") as f:
            for i, sequence in enumerate(dataset.sequences):
                f.write(f">seq_{i}\n{sequence}\n")

        return str(file_path)

    @staticmethod
    def load(
        filename: str, base_path: str | None = None
    ) -> tuple[Dataset, dict[str, Any]]:
        """Load dataset from file and also return a dict with used parameters."""
        base = FileDatasetRepository._get_base_path(base_path)
        file_path = FileDatasetRepository._resolve_path(filename, base)

        if not file_path.exists():
            raise DatasetNotFoundError(f"Dataset not found: {filename}")

        sequences = FileDatasetRepository._parse_fasta(file_path)

        params: dict[str, Any] = {"file_path": file_path}

        # Use filename without extension as default name and id
        dataset_name = file_path.stem
        dataset_id = f"file_{dataset_name}"

        return Dataset(id=dataset_id, name=dataset_name, sequences=sequences), params

    @staticmethod
    def list_available(base_path: str | None = None) -> List[str]:
        """List available datasets."""
        base = FileDatasetRepository._get_base_path(base_path)
        files = list(base.glob("*.fasta"))
        return [f.stem for f in files]

    @staticmethod
    def exists(filename: str, base_path: str | None = None) -> bool:
        """Check if dataset exists."""
        base = FileDatasetRepository._get_base_path(base_path)
        file_path = FileDatasetRepository._resolve_path(filename, base)
        return file_path.exists()

    @staticmethod
    def delete(filename: str, base_path: str | None = None) -> bool:
        """Remove dataset."""
        base = FileDatasetRepository._get_base_path(base_path)
        file_path = FileDatasetRepository._resolve_path(filename, base)
        if file_path.exists():
            file_path.unlink()
            return True
        return False

    @staticmethod
    def _resolve_path(filename: str, base_path: Path) -> Path:
        """Resolve filename to file path.

        If 'filename' is an absolute path, return it directly.
        Otherwise, treat it as a name relative to base path, with optional .fasta.
        """
        p = Path(filename)
        if p.is_absolute():
            return p
        if filename.endswith(".fasta"):
            return base_path / filename
        return base_path / f"{filename}.fasta"

    @staticmethod
    def _parse_fasta(file_path: Path) -> List[str]:
        """Parse FASTA file (multi-line sequences supported)."""
        sequences: List[str] = []
        current: list[str] = []

        with file_path.open("r", encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current:
                        sequences.append("".join(current))
                        current = []
                else:
                    current.append(line)
            if current:
                sequences.append("".join(current))

        if not sequences:
            raise DatasetValidationError(f"No sequences found in {file_path}")

        return sequences


# Backward compatibility alias for legacy imports
# Some modules/tests may still import FastaDatasetRepository
FastaDatasetRepository = FileDatasetRepository
