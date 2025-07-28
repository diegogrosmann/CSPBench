"""
Domain: Dataset

This module contains the main Dataset entity and its operations,
representing string sequences and associated metadata.
Pure implementation without external dependencies.
"""

import random
from typing import Any, Dict, List, Optional


class Dataset:
    """
    Dataset entity representing a set of strings for CSP.

    Attributes:
        sequences: List of dataset strings
        metadata: Dataset metadata (size, origin, etc.)
    """

    def __init__(self, sequences: List[str], metadata: Optional[Dict[str, Any]] = None):
        """
        Initialize a dataset with sequences and metadata.

        Args:
            sequences: List of strings
            metadata: Dictionary with optional metadata
        """
        if not sequences:
            raise ValueError("Dataset cannot be empty")

        # Validate that all strings have the same length
        length = len(sequences[0])
        if not all(len(seq) == length for seq in sequences):
            raise ValueError("All strings must have the same length")

        self.sequences = sequences
        self.metadata = metadata or {}

        # Infer basic metadata
        self.metadata.update(
            {
                "n": len(sequences),
                "L": length,
                "alphabet": self._infer_alphabet(),
                "diversity": self._calculate_diversity(),
            }
        )

    def _infer_alphabet(self) -> str:
        """Infer alphabet from sequences."""
        alphabet_set = set()
        for seq in self.sequences:
            alphabet_set.update(seq)
        return "".join(sorted(alphabet_set))

    def _calculate_diversity(self) -> float:
        """Calculate average dataset diversity."""
        if len(self.sequences) < 2:
            return 0.0

        total_distance = 0
        total_pairs = 0

        for i in range(len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                distance = sum(
                    c1 != c2 for c1, c2 in zip(self.sequences[i], self.sequences[j])
                )
                total_distance += distance
                total_pairs += 1

        avg_distance = total_distance / total_pairs if total_pairs > 0 else 0
        max_possible = len(self.sequences[0])
        return avg_distance / max_possible if max_possible > 0 else 0

    @property
    def size(self) -> int:
        """Return number of sequences."""
        return len(self.sequences)

    @property
    def length(self) -> int:
        """Return length of sequences."""
        return len(self.sequences[0]) if self.sequences else 0

    @property
    def alphabet(self) -> str:
        """Return dataset alphabet."""
        return self.metadata.get("alphabet", "")

    def validate(self) -> bool:
        """
        Validate dataset consistency.

        Returns:
            bool: True if dataset is valid
        """
        if not self.sequences:
            return False

        # Check uniform lengths
        length = len(self.sequences[0])
        if not all(len(seq) == length for seq in self.sequences):
            return False

        # Check valid characters
        alphabet_set = set(self.alphabet)
        for seq in self.sequences:
            if not all(c in alphabet_set for c in seq):
                return False

        return True

    def get_statistics(self) -> Dict[str, Any]:
        """
        Return detailed dataset statistics.

        Returns:
            dict: Dataset statistics
        """
        return {
            "size": self.size,
            "length": self.length,
            "alphabet": self.alphabet,
            "alphabet_size": len(self.alphabet),
            "diversity": self.metadata.get("diversity", 0),
            "total_characters": self.size * self.length,
            "metadata": self.metadata.copy(),
        }

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert dataset to dictionary.

        Returns:
            dict: Dictionary representation
        """
        return {"sequences": self.sequences.copy(), "metadata": self.metadata.copy()}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Dataset":
        """
        Create dataset from dictionary.

        Args:
            data: Dictionary with sequences and metadata

        Returns:
            Dataset: Created instance
        """
        return cls(data["sequences"], data.get("metadata"))

    def add_sequence(self, sequence: str) -> None:
        """
        Add new sequence to dataset.

        Args:
            sequence: String to be added

        Raises:
            ValueError: If sequence has different length
        """
        if len(sequence) != self.length:
            raise ValueError(f"Sequence must have length {self.length}")

        self.sequences.append(sequence)

        # Update metadata
        self.metadata["n"] = len(self.sequences)
        self.metadata["alphabet"] = self._infer_alphabet()
        self.metadata["diversity"] = self._calculate_diversity()

    def remove_sequence(self, index: int) -> str:
        """
        Remove sequence by index.

        Args:
            index: Index of sequence to be removed

        Returns:
            str: Removed sequence

        Raises:
            IndexError: If invalid index
        """
        if not 0 <= index < len(self.sequences):
            raise IndexError("Index out of range")

        removed = self.sequences.pop(index)

        # Update metadata
        self.metadata["n"] = len(self.sequences)
        if self.sequences:
            self.metadata["diversity"] = self._calculate_diversity()

        return removed

    def filter_by_pattern(self, pattern: str, position: int) -> "Dataset":
        """
        Filter sequences that have specific pattern at position.

        Args:
            pattern: Character or pattern to search
            position: Position to check

        Returns:
            Dataset: New dataset with filtered sequences
        """
        filtered_sequences = [seq for seq in self.sequences if seq[position] == pattern]

        new_metadata = self.metadata.copy()
        new_metadata["filter_applied"] = f"position_{position}={pattern}"

        return Dataset(filtered_sequences, new_metadata)

    def sample(self, n: int, seed: Optional[int] = None) -> "Dataset":
        """
        Return random sample from dataset.

        Args:
            n: Number of sequences in sample
            seed: Seed for reproducibility

        Returns:
            Dataset: New dataset with sample
        """
        if n > len(self.sequences):
            raise ValueError("Sample size larger than dataset")

        rng = random.Random(seed)
        sampled_sequences = rng.sample(self.sequences, n)

        new_metadata = self.metadata.copy()
        new_metadata["sampled_from"] = len(self.sequences)
        new_metadata["sample_seed"] = seed

        return Dataset(sampled_sequences, new_metadata)


class SyntheticDatasetGenerator:
    """Synthetic dataset generator for CSP algorithm testing."""

    @staticmethod
    def generate_from_center(
        center: str,
        n: int,
        noise_rate: float,
        alphabet: str,
        seed: Optional[int] = None,
    ) -> Dataset:
        """
        Generate dataset based on center string with noise.

        Args:
            center: Center string
            n: Number of sequences to generate
            noise_rate: Noise rate (0-1)
            alphabet: Valid alphabet
            seed: Seed for reproducibility

        Returns:
            Dataset: Generated synthetic dataset
        """
        rng = random.Random(seed)
        sequences = []

        for _ in range(n):
            sequence = list(center)

            # Aplicar ruído
            for i in range(len(sequence)):
                if rng.random() < noise_rate:
                    # Trocar por caractere diferente do alfabeto
                    available = [c for c in alphabet if c != sequence[i]]
                    if available:
                        sequence[i] = rng.choice(available)

            sequences.append("".join(sequence))

        metadata = {
            "type": "synthetic",
            "center_string": center,
            "noise_rate": noise_rate,
            "generation_seed": seed,
            "alphabet_used": alphabet,
        }

        return Dataset(sequences, metadata)

    @staticmethod
    def generate_random(
        n: int, length: int, alphabet: str, seed: Optional[int] = None
    ) -> Dataset:
        """
        Generate completely random dataset.

        Args:
            n: Number of sequences
            length: Length of sequences
            alphabet: Valid alphabet
            seed: Seed for reproducibility

        Returns:
            Dataset: Generated random dataset
        """
        rng = random.Random(seed)
        sequences = []

        for _ in range(n):
            sequence = "".join(rng.choice(alphabet) for _ in range(length))
            sequences.append(sequence)

        metadata = {
            "type": "random",
            "generation_seed": seed,
            "alphabet_used": alphabet,
        }

        return Dataset(sequences, metadata)

    @staticmethod
    def generate_clustered(
        n_clusters: int,
        sequences_per_cluster: int,
        length: int,
        alphabet: str,
        noise_rate: float = 0.1,
        seed: Optional[int] = None,
    ) -> Dataset:
        """
        Generate dataset with clusters of similar sequences.

        Args:
            n_clusters: Number of clusters
            sequences_per_cluster: Sequences per cluster
            length: Length of sequences
            alphabet: Valid alphabet
            noise_rate: Noise rate within clusters
            seed: Seed for reproducibility

        Returns:
            Dataset: Dataset with clusters
        """
        rng = random.Random(seed)
        all_sequences = []
        cluster_centers = []

        # Generate cluster centers
        for _ in range(n_clusters):
            center = "".join(rng.choice(alphabet) for _ in range(length))
            cluster_centers.append(center)

            # Generate sequences for this cluster
            cluster_dataset = SyntheticDatasetGenerator.generate_from_center(
                center,
                sequences_per_cluster,
                noise_rate,
                alphabet,
                rng.randint(0, 1000000),
            )
            all_sequences.extend(cluster_dataset.sequences)

        metadata = {
            "type": "clustered",
            "n_clusters": n_clusters,
            "sequences_per_cluster": sequences_per_cluster,
            "noise_rate": noise_rate,
            "cluster_centers": cluster_centers,
            "generation_seed": seed,
            "alphabet_used": alphabet,
        }

        return Dataset(all_sequences, metadata)
