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
    """

    def __init__(self, sequences: List[str] = [], alphabet: Optional[str] = None):
        """
        Initialize a dataset with sequences.

        Args:
            sequences: List of strings
            alphabet: Optional alphabet. If None, will be inferred from sequences
        """
        self.sequences = sequences
        self._alphabet = alphabet
        if not self.validate():
            raise ValueError(
                "Invalid dataset: sequences contain characters not in alphabet"
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
        """Return representative length.

        Para conjuntos não uniformes, retorna o maior comprimento.
        """
        if not self.sequences:
            return 0
        lengths = [len(s) for s in self.sequences]
        is_uniform = len(set(lengths)) == 1
        if is_uniform:
            return len(self.sequences[0])
        return max(lengths)

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

    def get_statistics(self) -> Dict[str, Any]:
        """
        Return detailed dataset statistics.

        Returns:
            dict: Dataset statistics
        """
        lengths = [len(s) for s in self.sequences]
        is_uniform = len(set(lengths)) == 1
        min_length = min(lengths) if lengths else 0
        max_length = max(lengths) if lengths else 0

        return {
            "size": self.size,
            "min_length": min_length,
            "max_length": max_length,
            "alphabet": self._infer_alphabet(),
            "alphabet_size": len(self._infer_alphabet()),
            "diversity": self._calculate_diversity(),
            "total_characters": sum(len(s) for s in self.sequences),
            "uniform_lengths": is_uniform,
            "lengths": lengths,
            "n": len(self.sequences),
            "L": max_length,
        }

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert dataset to dictionary.

        Returns:
            dict: Dictionary representation
        """
        return {"sequences": self.sequences.copy(), "metadata": self.get_statistics()}

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
        return cls(data["sequences"], alphabet)

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
        return Dataset(filtered_sequences)

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
        return Dataset(sampled_sequences)

    def get_sequences(self) -> List[str]:
        """
        Return list of sequences.

        Returns:
            List[str]: Copy of sequences list
        """
        return self.sequences.copy()
