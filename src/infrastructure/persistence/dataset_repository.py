"""
Repositório de Datasets baseado em arquivos

Implementa DatasetRepository usando sistema de arquivos.
"""

from pathlib import Path
from typing import List

from src.domain import Dataset
from src.domain.errors import DatasetNotFoundError, DatasetValidationError


class FileDatasetRepository:
    """Repositório de datasets baseado em arquivos FASTA."""

    def __init__(self, base_path: str = "saved_datasets"):
        self.base_path = Path(base_path)
        self.base_path.mkdir(exist_ok=True)

    def save(self, dataset: Dataset, name: str) -> str:
        """Salva dataset em arquivo FASTA."""
        file_path = self.base_path / f"{name}.fasta"

        with open(file_path, "w") as f:
            for i, sequence in enumerate(dataset.sequences):
                f.write(f">seq_{i}\n{sequence}\n")

        return str(file_path)

    def load(self, identifier: str) -> Dataset:
        """Carrega dataset de arquivo."""
        file_path = self._resolve_path(identifier)

        if not file_path.exists():
            raise DatasetNotFoundError(f"Dataset não encontrado: {identifier}")

        sequences = self._parse_fasta(file_path)
        metadata = {"source": str(file_path), "format": "fasta"}

        return Dataset(sequences, metadata)

    def list_available(self) -> List[str]:
        """Lista datasets disponíveis."""
        files = list(self.base_path.glob("*.fasta"))
        return [f.stem for f in files]

    def exists(self, identifier: str) -> bool:
        """Verifica se dataset existe."""
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
        """Resolve identificador para caminho do arquivo."""
        if identifier.endswith(".fasta"):
            return self.base_path / identifier
        return self.base_path / f"{identifier}.fasta"

    def _parse_fasta(self, file_path: Path) -> List[str]:
        """Parse de arquivo FASTA."""
        sequences = []
        current_sequence = ""

        with open(file_path, "r") as f:
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
            raise DatasetValidationError(f"Nenhuma sequência encontrada em {file_path}")

        return sequences


class FastaDatasetRepository(FileDatasetRepository):
    """Alias para compatibilidade."""

    pass
