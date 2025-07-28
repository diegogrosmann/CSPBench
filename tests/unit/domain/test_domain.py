"""
Testes unitários para o domínio da aplicação.

Testa entidades de domínio como Dataset, métricas e geradores.
"""

import pytest
import tempfile
import os
from pathlib import Path

from src.domain import (
    Dataset, 
    SyntheticDatasetGenerator,
    hamming_distance
)
from src.domain.errors import DatasetValidationError


def calculate_coverage(sequences, center):
    """Calculate coverage metrics for a center string."""
    distances = [hamming_distance(center, seq) for seq in sequences]
    return {
        "max_distance": max(distances),
        "avg_distance": sum(distances) / len(distances),
        "distances": distances
    }


def calculate_entropy(sequences):
    """Calculate entropy for sequences."""
    import math
    
    if not sequences:
        raise ValueError("Cannot calculate entropy for empty sequences")
    
    if len(sequences) == 1:
        return {
            "total_entropy": 0.0,
            "position_entropies": [0.0] * len(sequences[0])
        }
    
    position_entropies = []
    seq_length = len(sequences[0])
    
    for pos in range(seq_length):
        # Count character frequencies at this position
        char_counts = {}
        for seq in sequences:
            char = seq[pos]
            char_counts[char] = char_counts.get(char, 0) + 1
        
        # Calculate entropy for this position
        total_chars = len(sequences)
        entropy = 0.0
        
        for count in char_counts.values():
            if count > 0:
                prob = count / total_chars
                entropy -= prob * math.log2(prob)
        
        position_entropies.append(entropy)
    
    return {
        "total_entropy": sum(position_entropies),
        "position_entropies": position_entropies
    }


class TestDataset:
    """Tests for Dataset class."""

    def test_dataset_creation_from_sequences(self):
        """Test dataset creation from sequences."""
        sequences = ["ACGT", "ATGT", "ACCT"]
        metadata = {"source": "test", "alphabet": "ACGT"}
        
        dataset = Dataset(sequences=sequences, metadata=metadata)
        
        assert dataset.sequences == sequences
        assert dataset.metadata == metadata
        assert dataset.size == 3
        assert dataset.length == 4
        assert dataset.alphabet == "ACGT"

    def test_dataset_creation_from_file(self):
        """Test dataset creation from FASTA file."""
        fasta_content = """>seq1
ACGTACGT
>seq2
ATGTACGT
>seq3
ACGTCCGT
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_content)
            f.flush()
            
            try:
                dataset = Dataset.from_file(f.name)
                
                assert len(dataset.sequences) == 3
                assert dataset.sequences[0] == "ACGTACGT"
                assert dataset.sequences[1] == "ATGTACGT" 
                assert dataset.sequences[2] == "ACGTCCGT"
                assert dataset.size == 3
                assert dataset.length == 8
                assert "filename" in dataset.metadata
            finally:
                os.unlink(f.name)

    def test_dataset_with_empty_sequences(self):
        """Test dataset with empty sequences list."""
        with pytest.raises(DatasetValidationError):
            Dataset(sequences=[], metadata={})

    def test_dataset_with_different_length_sequences(self):
        """Test dataset with sequences of different lengths."""
        sequences = ["ACGT", "ATGTCC", "AC"]
        
        with pytest.raises(DatasetValidationError):
            Dataset(sequences=sequences, metadata={})

    def test_dataset_alphabet_detection(self):
        """Test automatic alphabet detection."""
        sequences = ["ACGT", "ATGT", "GCTA"]
        dataset = Dataset(sequences=sequences, metadata={})
        
        assert set(dataset.alphabet) == {"A", "C", "G", "T"}

    def test_dataset_metadata_defaults(self):
        """Test dataset with default metadata."""
        sequences = ["ACGT", "ATGT"]
        dataset = Dataset(sequences=sequences)
        
        assert isinstance(dataset.metadata, dict)
        assert dataset.metadata.get("size") == 2
        assert dataset.metadata.get("length") == 4

    def test_dataset_str_representation(self):
        """Test string representation of dataset."""
        sequences = ["ACGT", "ATGT"]
        dataset = Dataset(sequences=sequences, metadata={"name": "test"})
        
        str_repr = str(dataset)
        assert "2 sequences" in str_repr
        assert "length 4" in str_repr

    def test_dataset_equality(self):
        """Test dataset equality comparison."""
        sequences1 = ["ACGT", "ATGT"]
        sequences2 = ["ACGT", "ATGT"]
        sequences3 = ["ACGT", "ATCC"]
        
        dataset1 = Dataset(sequences=sequences1, metadata={"name": "test"})
        dataset2 = Dataset(sequences=sequences2, metadata={"name": "test"})
        dataset3 = Dataset(sequences=sequences3, metadata={"name": "test"})
        
        assert dataset1 == dataset2
        assert dataset1 != dataset3

    def test_dataset_get_sequence_by_index(self):
        """Test getting sequence by index."""
        sequences = ["ACGT", "ATGT", "GCTA"]
        dataset = Dataset(sequences=sequences, metadata={})
        
        assert dataset[0] == "ACGT"
        assert dataset[1] == "ATGT"
        assert dataset[2] == "GCTA"
        
        with pytest.raises(IndexError):
            _ = dataset[3]

    def test_dataset_iteration(self):
        """Test dataset iteration."""
        sequences = ["ACGT", "ATGT", "GCTA"]
        dataset = Dataset(sequences=sequences, metadata={})
        
        iterated_sequences = list(dataset)
        assert iterated_sequences == sequences

    def test_dataset_contains(self):
        """Test membership testing."""
        sequences = ["ACGT", "ATGT", "GCTA"]
        dataset = Dataset(sequences=sequences, metadata={})
        
        assert "ACGT" in dataset
        assert "ATGT" in dataset
        assert "TTTT" not in dataset

    def test_dataset_copy(self):
        """Test dataset copying."""
        sequences = ["ACGT", "ATGT"]
        metadata = {"name": "original"}
        original = Dataset(sequences=sequences, metadata=metadata)
        
        copied = original.copy()
        
        assert copied == original
        assert copied is not original
        assert copied.metadata is not original.metadata
        
        # Modifying copy shouldn't affect original
        copied.metadata["name"] = "copy"
        assert original.metadata["name"] == "original"


class TestSyntheticDatasetGenerator:
    """Tests for SyntheticDatasetGenerator."""

    def test_generate_random_basic(self):
        """Test basic random dataset generation."""
        dataset = SyntheticDatasetGenerator.generate_random(
            n=5, length=10, alphabet="ACGT", seed=42
        )
        
        assert dataset.size == 5
        assert dataset.length == 10
        assert all(len(seq) == 10 for seq in dataset.sequences)
        assert all(all(c in "ACGT" for c in seq) for seq in dataset.sequences)
        assert "generation_method" in dataset.metadata
        assert dataset.metadata["generation_method"] == "random"

    def test_generate_random_reproducible(self):
        """Test random generation reproducibility with seed."""
        dataset1 = SyntheticDatasetGenerator.generate_random(
            n=3, length=8, alphabet="ACGT", seed=123
        )
        dataset2 = SyntheticDatasetGenerator.generate_random(
            n=3, length=8, alphabet="ACGT", seed=123
        )
        
        assert dataset1.sequences == dataset2.sequences

    def test_generate_random_different_seeds(self):
        """Test random generation with different seeds."""
        dataset1 = SyntheticDatasetGenerator.generate_random(
            n=5, length=10, alphabet="ACGT", seed=42
        )
        dataset2 = SyntheticDatasetGenerator.generate_random(
            n=5, length=10, alphabet="ACGT", seed=43
        )
        
        # Should be different with different seeds
        assert dataset1.sequences != dataset2.sequences

    def test_generate_with_noise_basic(self):
        """Test generation with noise."""
        base_sequence = "ACGTACGTACGT"
        dataset = SyntheticDatasetGenerator.generate_with_noise(
            base_sequence=base_sequence,
            n=4,
            noise_rate=0.2,
            alphabet="ACGT",
            seed=42
        )
        
        assert dataset.size == 4
        assert dataset.length == len(base_sequence)
        assert all(len(seq) == len(base_sequence) for seq in dataset.sequences)
        assert "generation_method" in dataset.metadata
        assert dataset.metadata["generation_method"] == "noise_based"
        assert dataset.metadata["base_sequence"] == base_sequence
        assert dataset.metadata["noise_rate"] == 0.2

    def test_generate_with_noise_zero_noise(self):
        """Test generation with zero noise."""
        base_sequence = "ACGTACGT"
        dataset = SyntheticDatasetGenerator.generate_with_noise(
            base_sequence=base_sequence,
            n=3,
            noise_rate=0.0,
            alphabet="ACGT",
            seed=42
        )
        
        # With zero noise, all sequences should be identical to base
        assert all(seq == base_sequence for seq in dataset.sequences)

    def test_generate_with_noise_high_noise(self):
        """Test generation with high noise."""
        base_sequence = "AAAAAAAAAA"
        dataset = SyntheticDatasetGenerator.generate_with_noise(
            base_sequence=base_sequence,
            n=5,
            noise_rate=0.8,
            alphabet="ACGT",
            seed=42
        )
        
        # With high noise, sequences should be very different from base
        different_sequences = [
            seq for seq in dataset.sequences if seq != base_sequence
        ]
        assert len(different_sequences) > 0

    def test_generate_with_mutations(self):
        """Test generation with specific mutation types."""
        base_sequence = "ACGTACGTACGT"
        dataset = SyntheticDatasetGenerator.generate_with_mutations(
            base_sequence=base_sequence,
            n=3,
            mutation_types=["substitution", "insertion"],
            mutation_rate=0.1,
            alphabet="ACGT",
            seed=42
        )
        
        assert dataset.size == 3
        assert "generation_method" in dataset.metadata
        assert dataset.metadata["generation_method"] == "mutation_based"
        assert dataset.metadata["mutation_types"] == ["substitution", "insertion"]

    def test_generate_clustered(self):
        """Test clustered dataset generation."""
        dataset = SyntheticDatasetGenerator.generate_clustered(
            n=6,
            length=12,
            num_clusters=2,
            alphabet="ACGT",
            cluster_distance=0.3,
            seed=42
        )
        
        assert dataset.size == 6
        assert dataset.length == 12
        assert "generation_method" in dataset.metadata
        assert dataset.metadata["generation_method"] == "clustered"
        assert dataset.metadata["num_clusters"] == 2
        assert dataset.metadata["cluster_distance"] == 0.3

    def test_generator_invalid_parameters(self):
        """Test generator with invalid parameters."""
        with pytest.raises(ValueError):
            SyntheticDatasetGenerator.generate_random(n=0, length=10, alphabet="ACGT")
        
        with pytest.raises(ValueError):
            SyntheticDatasetGenerator.generate_random(n=5, length=0, alphabet="ACGT")
        
        with pytest.raises(ValueError):
            SyntheticDatasetGenerator.generate_random(n=5, length=10, alphabet="")
        
        with pytest.raises(ValueError):
            SyntheticDatasetGenerator.generate_with_noise(
                base_sequence="ACGT", n=5, noise_rate=-0.1, alphabet="ACGT"
            )
        
        with pytest.raises(ValueError):
            SyntheticDatasetGenerator.generate_with_noise(
                base_sequence="ACGT", n=5, noise_rate=1.1, alphabet="ACGT"
            )

    def test_generator_custom_alphabet(self):
        """Test generator with custom alphabet."""
        dataset = SyntheticDatasetGenerator.generate_random(
            n=3, length=8, alphabet="01", seed=42
        )
        
        assert dataset.size == 3
        assert dataset.length == 8
        assert all(all(c in "01" for c in seq) for seq in dataset.sequences)
        assert dataset.alphabet == "01"


class TestMetrics:
    """Tests for metrics functions."""

    def test_hamming_distance_basic(self):
        """Test basic Hamming distance calculation."""
        assert hamming_distance("ACGT", "ACGT") == 0
        assert hamming_distance("ACGT", "ATGT") == 1
        assert hamming_distance("ACGT", "TTTT") == 3
        assert hamming_distance("AAAA", "TTTT") == 4

    def test_hamming_distance_empty_strings(self):
        """Test Hamming distance with empty strings."""
        assert hamming_distance("", "") == 0

    def test_hamming_distance_different_lengths(self):
        """Test Hamming distance with different length strings."""
        with pytest.raises(ValueError):
            hamming_distance("ACGT", "ACG")
        
        with pytest.raises(ValueError):
            hamming_distance("ACG", "ACGTC")

    def test_hamming_distance_case_sensitivity(self):
        """Test Hamming distance case sensitivity."""
        assert hamming_distance("acgt", "ACGT") == 4  # Case sensitive
        assert hamming_distance("AcGt", "ACGT") == 2

    def test_calculate_coverage_basic(self):
        """Test basic coverage calculation."""
        sequences = ["ACGT", "ATGT", "GCGT"]
        center = "ACGT"
        
        coverage = calculate_coverage(sequences, center)
        
        assert isinstance(coverage, dict)
        assert "max_distance" in coverage
        assert "avg_distance" in coverage
        assert "distances" in coverage
        assert len(coverage["distances"]) == 3
        assert coverage["max_distance"] == 1  # GCGT differs by 1 from ACGT
        assert coverage["avg_distance"] == pytest.approx(2/3, rel=1e-3)

    def test_calculate_coverage_perfect_center(self):
        """Test coverage with perfect center."""
        sequences = ["ACGT", "ACGT", "ACGT"]
        center = "ACGT"
        
        coverage = calculate_coverage(sequences, center)
        
        assert coverage["max_distance"] == 0
        assert coverage["avg_distance"] == 0.0
        assert all(d == 0 for d in coverage["distances"])

    def test_calculate_coverage_single_sequence(self):
        """Test coverage with single sequence."""
        sequences = ["ACGT"]
        center = "ACGT"
        
        coverage = calculate_coverage(sequences, center)
        
        assert coverage["max_distance"] == 0
        assert coverage["avg_distance"] == 0.0
        assert coverage["distances"] == [0]

    def test_calculate_coverage_worst_case(self):
        """Test coverage in worst case scenario."""
        sequences = ["AAAA", "TTTT"]
        center = "CCCC"
        
        coverage = calculate_coverage(sequences, center)
        
        assert coverage["max_distance"] == 4
        assert coverage["avg_distance"] == 4.0
        assert all(d == 4 for d in coverage["distances"])

    def test_calculate_entropy_basic(self):
        """Test basic entropy calculation."""
        sequences = ["ACGT", "ATGT", "GCGT"]
        
        entropy = calculate_entropy(sequences)
        
        assert isinstance(entropy, dict)
        assert "total_entropy" in entropy
        assert "position_entropies" in entropy
        assert len(entropy["position_entropies"]) == 4  # Length of sequences
        assert entropy["total_entropy"] >= 0
        assert all(e >= 0 for e in entropy["position_entropies"])

    def test_calculate_entropy_uniform(self):
        """Test entropy with uniform distribution."""
        # All positions have all possible characters
        sequences = ["ACGT", "TGCA", "GTAC", "CATG"]
        
        entropy = calculate_entropy(sequences)
        
        # Each position should have maximum entropy for 4 characters
        expected_max = 2.0  # log2(4) = 2
        
        for pos_entropy in entropy["position_entropies"]:
            assert pos_entropy == pytest.approx(expected_max, rel=1e-3)

    def test_calculate_entropy_no_variation(self):
        """Test entropy with no variation."""
        sequences = ["AAAA", "AAAA", "AAAA"]
        
        entropy = calculate_entropy(sequences)
        
        # No variation means zero entropy
        assert entropy["total_entropy"] == 0.0
        assert all(e == 0.0 for e in entropy["position_entropies"])

    def test_calculate_entropy_partial_variation(self):
        """Test entropy with partial variation."""
        sequences = ["ACGT", "ACGT", "ATGT"]  # Only position 1 varies
        
        entropy = calculate_entropy(sequences)
        
        # Positions 0, 2, 3 should have zero entropy
        # Position 1 should have some entropy
        assert entropy["position_entropies"][0] == 0.0  # All A
        assert entropy["position_entropies"][1] > 0.0   # C and T
        assert entropy["position_entropies"][2] == 0.0  # All G
        assert entropy["position_entropies"][3] == 0.0  # All T

    def test_calculate_entropy_single_sequence(self):
        """Test entropy with single sequence."""
        sequences = ["ACGT"]
        
        entropy = calculate_entropy(sequences)
        
        # Single sequence has zero entropy
        assert entropy["total_entropy"] == 0.0
        assert all(e == 0.0 for e in entropy["position_entropies"])

    def test_calculate_entropy_empty_sequences(self):
        """Test entropy with empty sequences."""
        with pytest.raises(ValueError):
            calculate_entropy([])

    def test_metrics_integration(self):
        """Test integration of different metrics."""
        # Create a dataset with known properties
        sequences = ["ACGT", "ATGT", "GCGT", "ACCT"]
        center = "ACGT"
        
        # Calculate all metrics
        coverage = calculate_coverage(sequences, center)
        entropy = calculate_entropy(sequences)
        
        # Verify consistency
        assert len(coverage["distances"]) == len(sequences)
        assert len(entropy["position_entropies"]) == len(sequences[0])
        
        # Verify relationships
        assert coverage["max_distance"] >= coverage["avg_distance"]
        assert entropy["total_entropy"] >= 0
        
        # Check specific values for this dataset
        expected_distances = [0, 1, 1, 2]  # Distances from ACGT
        assert coverage["distances"] == expected_distances
        assert coverage["max_distance"] == 2
        assert coverage["avg_distance"] == 1.0


class TestDomainIntegration:
    """Integration tests for domain components."""

    def test_dataset_with_metrics(self):
        """Test dataset creation and metrics calculation."""
        # Generate synthetic dataset
        dataset = SyntheticDatasetGenerator.generate_random(
            n=5, length=8, alphabet="ACGT", seed=42
        )
        
        # Calculate metrics
        entropy = calculate_entropy(dataset.sequences)
        
        # Verify integration
        assert len(entropy["position_entropies"]) == dataset.length
        assert entropy["total_entropy"] >= 0
        
        # Test with a center string
        center = dataset.sequences[0]  # Use first sequence as center
        coverage = calculate_coverage(dataset.sequences, center)
        
        assert len(coverage["distances"]) == dataset.size
        assert coverage["distances"][0] == 0  # Distance to itself

    def test_noisy_dataset_properties(self):
        """Test properties of datasets generated with noise."""
        base_sequence = "ACGTACGTACGT"
        
        # Generate datasets with different noise levels
        low_noise = SyntheticDatasetGenerator.generate_with_noise(
            base_sequence, n=10, noise_rate=0.1, alphabet="ACGT", seed=42
        )
        high_noise = SyntheticDatasetGenerator.generate_with_noise(
            base_sequence, n=10, noise_rate=0.5, alphabet="ACGT", seed=42
        )
        
        # Calculate metrics
        low_entropy = calculate_entropy(low_noise.sequences)
        high_entropy = calculate_entropy(high_noise.sequences)
        
        low_coverage = calculate_coverage(low_noise.sequences, base_sequence)
        high_coverage = calculate_coverage(high_noise.sequences, base_sequence)
        
        # High noise should result in higher entropy and worse coverage
        assert high_entropy["total_entropy"] >= low_entropy["total_entropy"]
        assert high_coverage["avg_distance"] >= low_coverage["avg_distance"]

    def test_file_io_integration(self):
        """Test file I/O integration with domain objects."""
        # Create a synthetic dataset
        dataset = SyntheticDatasetGenerator.generate_random(
            n=3, length=6, alphabet="ACGT", seed=42
        )
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")
            f.write(dataset.sequences[0] + "\n")
            f.write(">seq2\n")
            f.write(dataset.sequences[1] + "\n")
            f.write(">seq3\n")
            f.write(dataset.sequences[2] + "\n")
            f.flush()
            
            try:
                # Read back from file
                loaded_dataset = Dataset.from_file(f.name)
                
                # Verify equivalence
                assert loaded_dataset.sequences == dataset.sequences
                assert loaded_dataset.size == dataset.size
                assert loaded_dataset.length == dataset.length
                
                # Verify metrics are consistent
                original_entropy = calculate_entropy(dataset.sequences)
                loaded_entropy = calculate_entropy(loaded_dataset.sequences)
                
                assert original_entropy == loaded_entropy
                
            finally:
                os.unlink(f.name)
