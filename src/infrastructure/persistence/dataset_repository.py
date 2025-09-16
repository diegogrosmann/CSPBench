"""
File-based Dataset Repository

Implements DatasetRepository using file system.
Respects DATASET_DIRECTORY from environment for base path.
"""

from pathlib import Path
from typing import Any, List

from src.domain import Dataset
from src.domain.errors import DatasetNotFoundError, DatasetValidationError
from src.infrastructure.utils.path_utils import get_dataset_directory


class FileDatasetRepository:
    """
    FASTA file-based dataset repository.

    This class provides methods to save, load, and manage datasets stored as FASTA files
    in the file system. It automatically handles path resolution and supports both
    relative and absolute file paths.
    """

    @staticmethod
    def _get_base_path(base_path: str | None = None) -> Path:
        """
        Get base path, using environment variable as default.

        Args:
            base_path: Optional base path. If None, uses environment variable

        Returns:
            Path: Absolute path to the base directory
        """
        if base_path is None:
            return get_dataset_directory()

        path = Path(base_path)

        # Convert to absolute path if relative
        if not path.is_absolute():
            path = path.resolve()

        path.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def save(dataset: Dataset, name: str, base_path: str | None = None) -> str:
        """
        Save dataset to FASTA file.

        Args:
            dataset: Dataset object containing sequences to save
            name: Name for the FASTA file (without extension)
            base_path: Optional base directory path

        Returns:
            str: Full path to the saved file
        """
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
        """
        Load dataset from file and return dataset with metadata.

        Args:
            filename: Name of the FASTA file to load
            base_path: Optional base directory path

        Returns:
            tuple: Dataset object and dictionary with used parameters

        Raises:
            DatasetNotFoundError: If the specified file doesn't exist
        """
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
        """
        List available datasets in the repository.

        Args:
            base_path: Optional base directory path

        Returns:
            List[str]: List of dataset names (without .fasta extension)
        """
        base = FileDatasetRepository._get_base_path(base_path)
        files = list(base.glob("*.fasta"))
        return [f.stem for f in files]

    @staticmethod
    def exists(filename: str, base_path: str | None = None) -> bool:
        """
        Check if dataset exists in the repository.

        Args:
            filename: Name of the FASTA file to check
            base_path: Optional base directory path

        Returns:
            bool: True if dataset exists, False otherwise
        """
        base = FileDatasetRepository._get_base_path(base_path)
        file_path = FileDatasetRepository._resolve_path(filename, base)
        return file_path.exists()

    @staticmethod
    def delete(filename: str, base_path: str | None = None) -> bool:
        """
        Remove dataset from the repository.

        Args:
            filename: Name of the FASTA file to delete
            base_path: Optional base directory path

        Returns:
            bool: True if dataset was deleted, False if it didn't exist
        """
        base = FileDatasetRepository._get_base_path(base_path)
        file_path = FileDatasetRepository._resolve_path(filename, base)
        if file_path.exists():
            file_path.unlink()
            return True
        return False

    @staticmethod
    def _resolve_path(filename: str, base_path: Path) -> Path:
        """
        Resolve filename to full file path.

        If 'filename' is an absolute path, return it directly.
        Otherwise, treat it as a name relative to base path, with optional .fasta extension.

        Args:
            filename: File name or path to resolve
            base_path: Base directory for relative paths

        Returns:
            Path: Resolved file path
        """
        p = Path(filename)
        if p.is_absolute():
            return p
        if filename.endswith(".fasta"):
            return base_path / filename
        return base_path / f"{filename}.fasta"

    @staticmethod
    def _parse_fasta(file_path: Path) -> List[str]:
        """
        Parse FASTA file supporting multi-line sequences.

        Args:
            file_path: Path to the FASTA file to parse

        Returns:
            List[str]: List of sequences found in the file

        Raises:
            DatasetValidationError: If no sequences are found in the file
        """
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
