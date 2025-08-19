"""
Domain: Dataset

Conjunto de sequências com metadados. Agora permite comprimentos variados
quando desejado (ex.: dados do Entrez sem uniformização).
"""

import random
from typing import Any, Dict, List, Optional


class Dataset:
    """
    Dataset entity representing a set of strings for CSP.

    Attributes:
        sequences: List of dataset strings
        metadata: Dataset metadata (size, origin, etc.)
        name: Optional dataset name for identification
    """

    def __init__(
        self, name: str, sequences: List[str] = [], alphabet: Optional[str] = None
    ):
        """
        Initialize a dataset with sequences.

        Args:
            sequences: List of strings
            alphabet: Optional alphabet. If None, will be inferred from sequences
            name: Optional dataset name for identification
        """
        self.sequences = sequences
        self._alphabet = alphabet
        self.name = name
        self._statistics_cache = None
        self._diversity_cache = None
        self._inferred_alphabet_cache = None
        if not self.validate():
            raise ValueError(
                "Invalid dataset: sequences contain characters not in alphabet"
            )

    def _infer_alphabet(self) -> str:
        """Infer alphabet from sequences."""
        if self._inferred_alphabet_cache is not None:
            return self._inferred_alphabet_cache

        alphabet_set = set()
        for seq in self.sequences:
            alphabet_set.update(seq)
        self._inferred_alphabet_cache = "".join(sorted(alphabet_set))
        return self._inferred_alphabet_cache

    def _calculate_diversity(self) -> float:
        """Calculate average dataset diversity."""
        if self._diversity_cache is not None:
            return self._diversity_cache

        if len(self.sequences) < 2:
            self._diversity_cache = 0.0
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
        diversity = avg_distance / max_possible if max_possible > 0 else 0
        self._diversity_cache = diversity
        return diversity

    @property
    def size(self) -> int:
        """Return number of sequences."""
        return len(self.sequences)

    @property
    def min_length(self) -> int:
        """Return minimum sequence length."""
        if not self.sequences:
            return 0
        return min(len(s) for s in self.sequences)

    @property
    def max_length(self) -> int:
        """Return maximum sequence length."""
        if not self.sequences:
            return 0
        return max(len(s) for s in self.sequences)

    @property
    def average_length(self) -> float:
        """Return average sequence length."""
        if not self.sequences:
            return 0.0
        total_length = sum(len(s) for s in self.sequences)
        return total_length / len(self.sequences)

    @property
    def alphabet(self) -> str:
        """Return dataset alphabet."""
        if self._alphabet is None:
            return self._infer_alphabet()
        return self._alphabet

    def validate(self) -> bool:
        """
        Validate dataset consistency.

        Returns:
            bool: True if dataset is valid
        """

        # Check valid characters only if alphabet was provided
        if self._alphabet is not None:
            alphabet_set = set(self._alphabet)
            for seq in self.sequences:
                if not all(c in alphabet_set for c in seq):
                    return False

        return True

    def _invalidate_cache(self) -> None:
        """Invalidate statistics and diversity cache."""
        self._statistics_cache = None
        self._diversity_cache = None
        self._inferred_alphabet_cache = None

    def get_statistics(self) -> Dict[str, Any]:
        """
        Return detailed dataset statistics.

        Returns:
            dict: Dataset statistics
        """
        if self._statistics_cache is not None:
            return self._statistics_cache

        lengths = [len(s) for s in self.sequences]
        is_uniform = len(set(lengths)) == 1

        self._statistics_cache = {
            "size": self.size,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "average_length": self.average_length,
            "alphabet": self._infer_alphabet(),
            "alphabet_size": len(self._infer_alphabet()),
            "diversity": self._calculate_diversity(),
            "total_characters": sum(len(s) for s in self.sequences),
            "uniform_lengths": is_uniform,
            "lengths": lengths,
            "n": len(self.sequences),
            "L": self.max_length,
        }
        return self._statistics_cache

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert dataset to dictionary.

        Returns:
            dict: Dictionary representation
        """
        metadata = self.get_statistics()
        if self.name:
            metadata["name"] = self.name
        return {"sequences": self.sequences.copy(), "metadata": metadata}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Dataset":
        """
        Create dataset from dictionary.

        Args:
            data: Dictionary with sequences and metadata

        Returns:
            Dataset: Created instance
        """
        alphabet = data.get("metadata", {}).get("alphabet")
        name = data.get("metadata", {}).get("name", "unnamed_dataset")
        return cls(name=name, sequences=data["sequences"], alphabet=alphabet)

    def add_sequence(self, sequence: str) -> None:
        """
        Add new sequence to dataset.

        Args:
            sequence: String to be added

        Raises:
            ValueError: If sequence contains invalid characters
        """
        # Validate sequence if alphabet is provided
        if self._alphabet is not None:
            alphabet_set = set(self._alphabet)
            if not all(c in alphabet_set for c in sequence):
                raise ValueError("Sequence contains characters not in alphabet")

        self.sequences.append(sequence)
        self._invalidate_cache()

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
        self._invalidate_cache()
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
        filtered_name = f"{self.name}_filtered" if self.name else "filtered_dataset"
        return Dataset(
            name=filtered_name, sequences=filtered_sequences, alphabet=self._alphabet
        )

    def get_sequences(self) -> List[str]:
        """
        Return list of sequences.

        Returns:
            List[str]: Copy of sequences list
        """
        return self.sequences.copy()
