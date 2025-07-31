"""
Dataset Generation Service - Application Layer

Orchestrates generation of synthetic and real datasets, coordinating
between different generators and repositories.
"""

from pathlib import Path
from typing import Any, Dict, Optional

from src.domain import Dataset, SyntheticDatasetGenerator
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository


class DatasetGenerationService:
    """Dataset generation service."""

    def __init__(self, dataset_repository: FileDatasetRepository):
        """
        Initialize the service.

        Args:
            dataset_repository: Repository for dataset persistence
        """
        self.dataset_repository = dataset_repository
        self.synthetic_generator = SyntheticDatasetGenerator()

    def generate_synthetic_dataset(self, params: Dict[str, Any]) -> Dataset:
        """
        Generate synthetic dataset based on parameters.

        Args:
            params: Generation parameters (n, length, alphabet, noise, seed)

        Returns:
            Generated dataset
        """
        return self.synthetic_generator.generate_random(
            n=params["n"],
            length=params["length"],
            alphabet=params["alphabet"],
            seed=params.get("seed"),
        )

    def download_real_dataset(self, params: Dict[str, Any]) -> Dataset:
        """
        Download real dataset based on parameters.

        Args:
            params: Download parameters (source, query, etc.)

        Returns:
            Downloaded dataset

        Raises:
            NotImplementedError: For unimplemented sources
        """
        source = params.get("source")

        if source == "ncbi":
            return self._download_from_ncbi(params)
        elif source == "file":
            return self._import_from_file(params)
        else:
            raise NotImplementedError(f"Source '{source}' not implemented")

    def save_dataset(
        self, dataset: Dataset, filename: str, base_path: str = "datasets"
    ) -> str:
        """
        Save dataset to file.

        Args:
            dataset: Dataset to save
            filename: File name
            base_path: Base directory

        Returns:
            Complete path of saved file
        """
        output_path = Path(base_path) / filename

        # Create directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save in FASTA format
        with open(output_path, "w", encoding="utf-8") as f:
            for i, seq in enumerate(dataset.sequences):
                f.write(f">seq_{i+1}\n{seq}\n")

        return str(output_path.absolute())

    def _download_from_ncbi(self, params: Dict[str, Any]) -> Dataset:
        """
        Download dataset from NCBI.

        Args:
            params: NCBI parameters (query, max_sequences, min_length, max_length)

        Returns:
            Downloaded dataset
        """
        # For now, simulate download with synthetic data
        # Real implementation would use Bio.Entrez or similar

        print(f"🌐 Simulating NCBI download...")
        print(f"   Query: {params['query']}")
        print(f"   Max sequences: {params['max_sequences']}")
        print(f"   Length range: {params['min_length']}-{params['max_length']}")

        # Simulate with synthetic dataset based on parameters
        synthetic_params = {
            "n": params["max_sequences"],
            "length": (params["min_length"] + params["max_length"]) // 2,
            "alphabet": "ACTG",  # Standard DNA
            "noise": 0.05,  # Low noise for "real" data
            "seed": 42,  # Deterministic for simulation
        }

        return self.generate_synthetic_dataset(synthetic_params)

    def _import_from_file(self, params: Dict[str, Any]) -> Dataset:
        """
        Import dataset from file.

        Args:
            params: File parameters (file_path)

        Returns:
            Imported dataset
        """
        file_path = params["file_path"]

        print(f"📁 Importing file: {file_path}")

        # Use repository to load file
        dataset_name = Path(file_path).stem
        return self.dataset_repository.load(dataset_name)

    def get_generation_summary(
        self, dataset: Dataset, params: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Generate dataset generation summary.

        Args:
            dataset: Generated dataset
            params: Parameters used

        Returns:
            Generation summary
        """
        return {
            "total_sequences": len(dataset.sequences),
            "sequence_length": len(dataset.sequences[0]) if dataset.sequences else 0,
            "alphabet": dataset.alphabet,
            "parameters_used": params,
            "estimated_size_kb": sum(len(seq) for seq in dataset.sequences) / 1024,
        }
