"""
Dataset Generation Service - Application Layer

Orchestrates generation of synthetic and real datasets, coordinating
between different generators and repositories.
"""

from pathlib import Path
from typing import Any, Dict, List

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
            params: Generation parameters including method and method-specific params

        Returns:
            Generated dataset
        """
        method = params.get("method", "random")

        if method == "random":
            return self.synthetic_generator.generate_random(
                n=params["n"],
                length=params["length"],
                alphabet=params["alphabet"],
                seed=params.get("seed"),
            )
        elif method == "noise":
            center_sequence = params.get("center_sequence")
            if not center_sequence:
                # Generate random center sequence if not provided
                import random

                rng = random.Random(params.get("seed"))
                center_sequence = "".join(
                    rng.choice(params["alphabet"]) for _ in range(params["length"])
                )

            return self.synthetic_generator.generate_from_center(
                center=center_sequence,
                n=params["n"],
                noise_rate=params.get("noise_rate", 0.1),
                alphabet=params["alphabet"],
                seed=params.get("seed"),
            )
        elif method == "clustered":
            # Calculate sequences per cluster
            num_clusters = params.get("num_clusters", 3)
            sequences_per_cluster = max(1, params["n"] // num_clusters)

            return self.synthetic_generator.generate_clustered(
                n_clusters=num_clusters,
                sequences_per_cluster=sequences_per_cluster,
                length=params["length"],
                alphabet=params["alphabet"],
                noise_rate=params.get("cluster_noise", 0.1),
                seed=params.get("seed"),
            )
        elif method == "mutations":
            base_sequence = params.get("base_sequence")
            if not base_sequence:
                # Generate random base sequence if not provided
                import random

                rng = random.Random(params.get("seed"))
                base_sequence = "".join(
                    rng.choice(params["alphabet"]) for _ in range(params["length"])
                )

            # For now, use noise-based generation as mutation approximation
            # TODO: Implement proper mutation-based generation
            return self.synthetic_generator.generate_from_center(
                center=base_sequence,
                n=params["n"],
                noise_rate=params.get("mutation_rate", 0.1),
                alphabet=params["alphabet"],
                seed=params.get("seed"),
            )
        else:
            # Fallback to random generation
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
        Save dataset to file with unique filename.

        Args:
            dataset: Dataset to save
            filename: File name
            base_path: Base directory

        Returns:
            Complete path of saved file
        """
        base_dir = Path(base_path)
        base_dir.mkdir(parents=True, exist_ok=True)

        # Generate unique filename if file already exists
        unique_filename = self._generate_unique_filename(filename, base_dir)
        output_path = base_dir / unique_filename

        # Save in FASTA format
        with open(output_path, "w", encoding="utf-8") as f:
            for i, seq in enumerate(dataset.sequences):
                f.write(f">seq_{i+1}\n{seq}\n")

        return str(output_path.absolute())

    def _generate_unique_filename(self, base_filename: str, directory: Path) -> str:
        """Generate a unique filename by adding a counter if the file already exists."""
        filename = base_filename
        counter = 1

        # Extract name and extension
        if "." in filename:
            name_part, ext_part = filename.rsplit(".", 1)
        else:
            name_part, ext_part = filename, ""

        # Keep trying until we find a unique name
        while (directory / filename).exists():
            if ext_part:
                filename = f"{name_part}_{counter}.{ext_part}"
            else:
                filename = f"{name_part}_{counter}"
            counter += 1

        return filename

    def save_sequences_as_fasta(
        self, sequences: List[str], filename: str, base_path: str = "datasets"
    ) -> str:
        """
        Save sequences directly as FASTA file with unique filename.

        Args:
            sequences: List of sequences to save
            filename: File name
            base_path: Base directory

        Returns:
            Complete path of saved file
        """
        base_dir = Path(base_path)
        base_dir.mkdir(parents=True, exist_ok=True)

        # Generate unique filename if file already exists
        unique_filename = self._generate_unique_filename(filename, base_dir)
        output_path = base_dir / unique_filename

        # Save in FASTA format
        with open(output_path, "w", encoding="utf-8") as f:
            for i, seq in enumerate(sequences):
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

        print("🌐 Simulating NCBI download...")
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
