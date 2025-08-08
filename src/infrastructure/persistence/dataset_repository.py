"""
File-based Dataset Repository

Implements DatasetRepository using file system.
"""

from pathlib import Path
from typing import List

from src.domain import Dataset
from src.domain.errors import DatasetNotFoundError, DatasetValidationError


class FileDatasetRepository:
    """FASTA file-based dataset repository."""

    def __init__(self, base_path: str = "saved_datasets"):
        self.base_path = Path(base_path)
        self.base_path.mkdir(exist_ok=True)

    def save(self, dataset: Dataset, name: str) -> str:
        """Save dataset to FASTA file."""
        file_path = self.base_path / f"{name}.fasta"

        with open(file_path, "w") as f:
            for i, sequence in enumerate(dataset.sequences):
                f.write(f">seq_{i}\n{sequence}\n")

        return str(file_path)

    def load(self, identifier: str) -> Dataset:
        """Load dataset from file."""
        file_path = self._resolve_path(identifier)

        if not file_path.exists():
            raise DatasetNotFoundError(f"Dataset not found: {identifier}")

        sequences = self._parse_fasta(file_path)
        metadata = {"source": str(file_path), "format": "fasta"}

        return Dataset(sequences, metadata)

    def list_available(self) -> List[str]:
        """List available datasets."""
        files = list(self.base_path.glob("*.fasta"))
        return [f.stem for f in files]

    def exists(self, identifier: str) -> bool:
        """Check if dataset exists."""
        file_path = self._resolve_path(identifier)
        return file_path.exists()

    def delete(self, identifier: str) -> bool:
        """Remove dataset."""
        file_path = self._resolve_path(identifier)
        if file_path.exists():
            file_path.unlink()
            return True
        return False

    def _resolve_path(self, identifier: str) -> Path:
        """Resolve identifier to file path."""
        if identifier.endswith(".fasta"):
            return self.base_path / identifier
        return self.base_path / f"{identifier}.fasta"

    def _parse_fasta(self, file_path: Path) -> List[str]:
        """Parse FASTA file."""
        sequences = []
        current_sequence = ""

        with open(file_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_sequence:
                        sequences.append(current_sequence)
                        current_sequence = ""
                else:
                    current_sequence += line

            if current_sequence:
                sequences.append(current_sequence)

        if not sequences:
            raise DatasetValidationError(f"No sequences found in {file_path}")

        return sequences


class FastaDatasetRepository(FileDatasetRepository):
    """Alias for compatibility."""

    pass
