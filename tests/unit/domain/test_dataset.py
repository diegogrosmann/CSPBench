"""
Unit tests for Domain Dataset class.

Tests the core Dataset entity with various scenarios including
creation, validation, manipulation, and edge cases.
"""

import os
import tempfile
from pathlib import Path
from typing import List

import pytest

from src.domain.dataset import Dataset


class TestDatasetCreation:
    """Tests for dataset creation and initialization."""

    def test_create_empty_dataset(self):
        """Test creating an empty dataset."""
        dataset = Dataset()

        assert dataset.sequences == []
        assert dataset.size == 0
        assert dataset.length == 0
        assert dataset.alphabet == ""

    def test_create_dataset_with_sequences(self):
        """Test creating a dataset with sequences."""
        sequences = ["ACGT", "ATGT", "GCGT"]
        dataset = Dataset(sequences=sequences)

        assert dataset.sequences == sequences
        assert dataset.size == 3
        assert dataset.length == 4
        assert set(dataset.alphabet) == {"A", "C", "G", "T"}

    def test_create_dataset_with_explicit_alphabet(self):
        """Test creating a dataset with explicit alphabet."""
        sequences = ["ACG", "ATG", "GCG"]
        alphabet = "ACGT"
        dataset = Dataset(sequences=sequences, alphabet=alphabet)

        assert dataset.alphabet == alphabet
        assert dataset.validate()

    def test_create_dataset_invalid_alphabet(self):
        """Test creating dataset with invalid alphabet raises error."""
        sequences = ["ACGX", "ATGT"]  # X not in alphabet
        alphabet = "ACGT"

        with pytest.raises(
            ValueError,
            match="Invalid dataset: sequences contain characters not in alphabet",
        ):
            Dataset(sequences=sequences, alphabet=alphabet)

    def test_create_dataset_mixed_case(self):
        """Test dataset with mixed case sequences."""
        sequences = ["acgtACGT", "ATGTatgt", "GcGtCgCt"]
        dataset = Dataset(sequences=sequences)

        # Should preserve case
        assert dataset.sequences == sequences
        assert set(dataset.alphabet) == {"A", "C", "G", "T", "a", "c", "g", "t"}


class TestDatasetProperties:
    """Tests for dataset properties and computed values."""

    def test_size_property(self):
        """Test size property."""
        assert Dataset([]).size == 0
        assert Dataset(["A"]).size == 1
        assert Dataset(["A", "T", "G"]).size == 3

    def test_length_property_uniform(self):
        """Test length property with uniform sequences."""
        sequences = ["ACGT", "ATGT", "GCGT"]
        dataset = Dataset(sequences=sequences)

        assert dataset.length == 4

    def test_length_property_non_uniform(self):
        """Test length property with non-uniform sequences."""
        sequences = ["ACG", "ATGT", "GCGTACG"]  # lengths 3, 4, 7
        dataset = Dataset(sequences=sequences)

        # Should return max length for non-uniform
        assert dataset.length == 7

    def test_alphabet_inference(self):
        """Test alphabet inference from sequences."""
        sequences = ["ACGT", "CGTA", "TACG"]
        dataset = Dataset(sequences=sequences)

        # Should be sorted
        assert dataset.alphabet == "ACGT"

    def test_alphabet_inference_empty(self):
        """Test alphabet inference with empty dataset."""
        dataset = Dataset([])
        assert dataset.alphabet == ""

    def test_alphabet_explicit_vs_inferred(self):
        """Test explicit alphabet vs inferred."""
        sequences = ["ACG", "ATG"]

        # Inferred
        dataset1 = Dataset(sequences=sequences)
        assert dataset1.alphabet == "ACGT"

        # Explicit (different order)
        dataset2 = Dataset(sequences=sequences, alphabet="TAGC")
        assert dataset2.alphabet == "TAGC"


class TestDatasetStatistics:
    """Tests for dataset statistics calculation."""

    def test_get_statistics_basic(self):
        """Test basic statistics calculation."""
        sequences = ["ACGT", "ATGT", "GCGT"]
        dataset = Dataset(sequences=sequences)

        stats = dataset.get_statistics()

        assert stats["size"] == 3
        assert stats["min_length"] == 4
        assert stats["max_length"] == 4
        assert stats["alphabet"] == "ACGT"
        assert stats["alphabet_size"] == 4
        assert stats["uniform_lengths"] is True
        assert stats["total_characters"] == 12
        assert stats["n"] == 3
        assert stats["L"] == 4

        # Diversity should be > 0 (sequences differ)
        assert 0 < stats["diversity"] <= 1

    def test_get_statistics_non_uniform(self):
        """Test statistics with non-uniform lengths."""
        sequences = ["AC", "ATGT", "GCGTAC"]
        dataset = Dataset(sequences=sequences)

        stats = dataset.get_statistics()

        assert stats["min_length"] == 2
        assert stats["max_length"] == 6
        assert stats["uniform_lengths"] is False
        assert stats["lengths"] == [2, 4, 6]
        assert stats["total_characters"] == 12

    def test_get_statistics_single_sequence(self):
        """Test statistics with single sequence."""
        dataset = Dataset(sequences=["ACGTACGT"])

        stats = dataset.get_statistics()

        assert stats["size"] == 1
        assert stats["diversity"] == 0.0  # Single sequence has no diversity

    def test_get_statistics_identical_sequences(self):
        """Test statistics with identical sequences."""
        dataset = Dataset(sequences=["ACGT", "ACGT", "ACGT"])

        stats = dataset.get_statistics()

        assert stats["diversity"] == 0.0  # Identical sequences


class TestDatasetManipulation:
    """Tests for dataset manipulation methods."""

    def test_add_sequence(self):
        """Test adding sequences to dataset."""
        dataset = Dataset(["ACGT"])

        dataset.add_sequence("ATGT")

        assert len(dataset.sequences) == 2
        assert "ATGT" in dataset.sequences

    def test_add_sequence_invalid_alphabet(self):
        """Test adding sequence with invalid character."""
        dataset = Dataset(["ACGT"], alphabet="ACGT")

        with pytest.raises(
            ValueError, match="Sequence contains characters not in alphabet"
        ):
            dataset.add_sequence("ACGX")

    def test_remove_sequence(self):
        """Test removing sequences by index."""
        sequences = ["ACGT", "ATGT", "GCGT"]
        dataset = Dataset(sequences=sequences)

        removed = dataset.remove_sequence(1)

        assert removed == "ATGT"
        assert len(dataset.sequences) == 2
        assert "ATGT" not in dataset.sequences

    def test_remove_sequence_invalid_index(self):
        """Test removing with invalid index."""
        dataset = Dataset(["ACGT"])

        with pytest.raises(IndexError, match="Index out of range"):
            dataset.remove_sequence(5)

        with pytest.raises(IndexError, match="Index out of range"):
            dataset.remove_sequence(-2)

    def test_filter_by_pattern(self):
        """Test filtering sequences by pattern."""
        sequences = ["ACGT", "ATGT", "GCGT", "ATAT"]
        dataset = Dataset(sequences=sequences)

        filtered = dataset.filter_by_pattern("A", 0)  # First position is A

        assert len(filtered.sequences) == 3  # ACGT, ATGT, ATAT
        assert all(seq[0] == "A" for seq in filtered.sequences)

    def test_sample(self):
        """Test sampling from dataset."""
        sequences = ["ACGT", "ATGT", "GCGT", "ATAT", "CGTA"]
        dataset = Dataset(sequences=sequences)

        sample = dataset.sample(3, seed=42)

        assert len(sample.sequences) == 3
        assert all(seq in sequences for seq in sample.sequences)

    def test_sample_too_large(self):
        """Test sampling more than dataset size."""
        dataset = Dataset(["ACGT", "ATGT"])

        with pytest.raises(ValueError, match="Sample size larger than dataset"):
            dataset.sample(5)

    def test_sample_reproducible(self):
        """Test that sampling is reproducible with same seed."""
        sequences = ["ACGT", "ATGT", "GCGT", "ATAT", "CGTA"]
        dataset = Dataset(sequences=sequences)

        sample1 = dataset.sample(3, seed=42)
        sample2 = dataset.sample(3, seed=42)

        assert sample1.sequences == sample2.sequences


class TestDatasetSerialization:
    """Tests for dataset serialization/deserialization."""

    def test_to_dict(self):
        """Test converting dataset to dictionary."""
        sequences = ["ACGT", "ATGT"]
        dataset = Dataset(sequences=sequences)

        data = dataset.to_dict()

        assert "sequences" in data
        assert "metadata" in data
        assert data["sequences"] == sequences
        assert isinstance(data["metadata"], dict)

    def test_from_dict(self):
        """Test creating dataset from dictionary."""
        data = {"sequences": ["ACGT", "ATGT"], "metadata": {"alphabet": "ACGT"}}

        dataset = Dataset.from_dict(data)

        assert dataset.sequences == data["sequences"]
        assert dataset.alphabet == "ACGT"

    def test_from_dict_without_metadata(self):
        """Test creating dataset from dict without metadata."""
        data = {"sequences": ["ACGT", "ATGT"]}

        dataset = Dataset.from_dict(data)

        assert dataset.sequences == data["sequences"]

    def test_round_trip_serialization(self):
        """Test full serialization round trip."""
        original = Dataset(["ACGT", "ATGT", "GCGT"])

        data = original.to_dict()
        reconstructed = Dataset.from_dict(data)

        assert original.sequences == reconstructed.sequences


class TestDatasetValidation:
    """Tests for dataset validation."""

    def test_validate_valid_dataset(self):
        """Test validation of valid dataset."""
        dataset = Dataset(["ACGT", "ATGT"], alphabet="ACGT")

        assert dataset.validate() is True

    def test_validate_no_alphabet_specified(self):
        """Test validation when no alphabet is specified."""
        dataset = Dataset(["ACGT", "ATGX"])  # X not in inferred alphabet

        # Should still be valid as no explicit alphabet constraint
        assert dataset.validate() is True

    def test_validate_with_explicit_alphabet(self):
        """Test validation with explicit alphabet constraint."""
        sequences = ["ACGT", "ATGT"]
        dataset = Dataset(sequences, alphabet="ACGT")

        assert dataset.validate() is True


class TestDatasetUtilityMethods:
    """Tests for utility methods."""

    def test_get_sequences(self):
        """Test getting copy of sequences."""
        sequences = ["ACGT", "ATGT"]
        dataset = Dataset(sequences=sequences)

        copied_sequences = dataset.get_sequences()

        assert copied_sequences == sequences
        assert copied_sequences is not sequences  # Should be a copy

    def test_get_sequences_modification(self):
        """Test that modifying returned sequences doesn't affect original."""
        sequences = ["ACGT", "ATGT"]
        dataset = Dataset(sequences=sequences)

        copied_sequences = dataset.get_sequences()
        copied_sequences.append("GCGT")

        assert len(dataset.sequences) == 2  # Original unchanged


class TestDatasetEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_empty_sequences_list(self):
        """Test creating dataset with empty sequences list."""
        dataset = Dataset(sequences=[])

        assert dataset.size == 0
        assert dataset.length == 0
        assert dataset.alphabet == ""

    def test_sequences_with_empty_strings(self):
        """Test dataset with empty string sequences."""
        dataset = Dataset(sequences=["", "", ""])

        assert dataset.size == 3
        assert dataset.length == 0
        assert dataset.alphabet == ""

    def test_mixed_empty_and_non_empty_sequences(self):
        """Test dataset with mix of empty and non-empty sequences."""
        sequences = ["", "ACGT", "", "AT"]
        dataset = Dataset(sequences=sequences)

        assert dataset.size == 4
        assert dataset.length == 4  # Max length
        assert set(dataset.alphabet) == {"A", "C", "G", "T"}

    def test_unicode_sequences(self):
        """Test dataset with unicode characters."""
        sequences = ["αβγδ", "αγβδ"]
        dataset = Dataset(sequences=sequences)

        assert dataset.size == 2
        assert dataset.length == 4
        assert set(dataset.alphabet) == {"α", "β", "γ", "δ"}

    def test_numeric_sequences(self):
        """Test dataset with numeric characters."""
        sequences = ["0123", "0213", "1023"]
        dataset = Dataset(sequences=sequences)

        assert dataset.size == 3
        assert dataset.length == 4
        assert set(dataset.alphabet) == {"0", "1", "2", "3"}
