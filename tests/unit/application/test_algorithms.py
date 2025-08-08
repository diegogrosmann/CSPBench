"""
Testes unitários para algoritmos CSP.

Testa funcionalidades básicas dos algoritmos implementados.
"""

from unittest.mock import Mock, patch

import pytest

from algorithms.baseline.algorithm import BaselineAlg
from algorithms.blf_ga.algorithm import BLFGAAlgorithm
from algorithms.csc.algorithm import CSCAlgorithm
from algorithms.dp_csp.algorithm import DPCSPAlgorithm
from algorithms.h3_csp.algorithm import H3CSPAlgorithm
from src.domain import Dataset, SyntheticDatasetGenerator


class TestBaselineAlgorithm:
    """Tests for BaselineAlg."""

    @pytest.fixture
    def small_dataset(self):
        """Small dataset for testing."""
        return SyntheticDatasetGenerator.generate_random(
            n=5, length=10, alphabet="ACGT", seed=42
        )

    @pytest.fixture
    def algorithm(self, small_dataset):
        """Baseline algorithm instance."""
        return BaselineAlg(small_dataset.sequences, "ACGT")

    def test_algorithm_initialization(self, algorithm):
        """Test algorithm initialization."""
        assert algorithm.name == "Baseline"
        assert algorithm.alphabet == "ACGT"
        assert len(algorithm.strings) == 5
        assert algorithm.is_deterministic is True

    def test_algorithm_run(self, algorithm):
        """Test algorithm execution."""
        result, distance, metadata = algorithm.run()

        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert isinstance(metadata, dict)
        assert distance >= 0
        assert len(result) == len(algorithm.strings[0])
        assert all(c in algorithm.alphabet for c in result)

    def test_algorithm_with_different_tie_break(self, small_dataset):
        """Test algorithm with different tie breaking strategies."""
        algorithms = {
            "lex": BaselineAlg(small_dataset.sequences, "ACGT", tie_break="lex"),
            "random": BaselineAlg(small_dataset.sequences, "ACGT", tie_break="random"),
            "first": BaselineAlg(small_dataset.sequences, "ACGT", tie_break="first"),
        }

        for name, algo in algorithms.items():
            result, distance, metadata = algo.run()
            assert isinstance(result, str)
            assert isinstance(distance, int)
            assert "tie_break" in metadata
            assert metadata["tie_break"] == name

    def test_algorithm_empty_strings(self):
        """Test algorithm with empty strings list."""
        with pytest.raises(ValueError):
            BaselineAlg([], "ACGT")

    def test_algorithm_invalid_alphabet(self, small_dataset):
        """Test algorithm with invalid alphabet."""
        with pytest.raises(ValueError):
            BaselineAlg(small_dataset.sequences, "")


class TestBLFGAAlgorithm:
    """Tests for BLF-GA Algorithm."""

    @pytest.fixture
    def small_dataset(self):
        """Small dataset for testing."""
        return SyntheticDatasetGenerator.generate_random(
            n=8, length=20, alphabet="ACGT", seed=42
        )

    @pytest.fixture
    def algorithm(self, small_dataset):
        """BLF-GA algorithm instance."""
        return BLFGAAlgorithm(
            small_dataset.sequences,
            "ACGT",
            pop_size=20,
            max_gens=10,  # Small for testing
            cross_prob=0.8,
            mut_prob=0.1,
        )

    def test_algorithm_initialization(self, algorithm):
        """Test algorithm initialization."""
        assert algorithm.name == "BLF-GA"
        assert algorithm.alphabet == "ACGT"
        assert algorithm.pop_size == 20
        assert algorithm.max_gens == 10

    def test_algorithm_run(self, algorithm):
        """Test algorithm execution."""
        result, distance, metadata = algorithm.run()

        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert isinstance(metadata, dict)
        assert distance >= 0
        assert len(result) == len(algorithm.strings[0])
        assert all(c in algorithm.alphabet for c in result)
        assert "generations" in metadata
        assert "final_population_size" in metadata

    def test_algorithm_with_different_params(self, small_dataset):
        """Test algorithm with different parameters."""
        algorithm = BLFGAAlgorithm(
            small_dataset.sequences,
            "ACGT",
            pop_size=10,
            max_gens=5,
            cross_prob=0.9,
            mut_prob=0.2,
            elite_rate=0.1,
            initial_blocks=0.3,
            crossover_type="uniform",
            mutation_type="inversion",
        )

        result, distance, metadata = algorithm.run()
        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert "crossover_type" in metadata
        assert metadata["crossover_type"] == "uniform"

    def test_algorithm_deterministic_with_seed(self, small_dataset):
        """Test algorithm determinism with seed."""
        algorithm1 = BLFGAAlgorithm(
            small_dataset.sequences, "ACGT", pop_size=10, max_gens=5, seed=123
        )
        algorithm2 = BLFGAAlgorithm(
            small_dataset.sequences, "ACGT", pop_size=10, max_gens=5, seed=123
        )

        result1, distance1, _ = algorithm1.run()
        result2, distance2, _ = algorithm2.run()

        # With same seed, should get same results
        assert result1 == result2
        assert distance1 == distance2


class TestCSCAlgorithm:
    """Tests for CSC Algorithm."""

    @pytest.fixture
    def small_dataset(self):
        """Small dataset for testing."""
        return SyntheticDatasetGenerator.generate_random(
            n=6, length=16, alphabet="ACGT", seed=42
        )

    @pytest.fixture
    def algorithm(self, small_dataset):
        """CSC algorithm instance."""
        return CSCAlgorithm(
            small_dataset.sequences,
            "ACGT",
            min_d=1,
            d_factor=0.8,
            min_blocks=2,
            max_blocks=4,
        )

    def test_algorithm_initialization(self, algorithm):
        """Test algorithm initialization."""
        assert algorithm.name == "CSC"
        assert algorithm.alphabet == "ACGT"
        assert algorithm.min_d == 1
        assert algorithm.d_factor == 0.8

    def test_algorithm_run(self, algorithm):
        """Test algorithm execution."""
        result, distance, metadata = algorithm.run()

        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert isinstance(metadata, dict)
        assert distance >= 0
        assert len(result) == len(algorithm.strings[0])
        assert all(c in algorithm.alphabet for c in result)
        assert "blocks_used" in metadata
        assert "iterations" in metadata

    def test_algorithm_with_different_params(self, small_dataset):
        """Test algorithm with different parameters."""
        algorithm = CSCAlgorithm(
            small_dataset.sequences,
            "ACGT",
            min_d=2,
            d_factor=0.75,
            min_blocks=3,
            max_blocks=6,
            l_div=20,
        )

        result, distance, metadata = algorithm.run()
        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert "min_d" in metadata
        assert metadata["min_d"] == 2


class TestH3CSPAlgorithm:
    """Tests for H³-CSP Algorithm."""

    @pytest.fixture
    def small_dataset(self):
        """Small dataset for testing."""
        return SyntheticDatasetGenerator.generate_random(
            n=6, length=12, alphabet="ACGT", seed=42
        )

    @pytest.fixture
    def algorithm(self, small_dataset):
        """H³-CSP algorithm instance."""
        return H3CSPAlgorithm(
            small_dataset.sequences,
            "ACGT",
            auto_blocks=True,
            min_block_size=2,
            k_candidates=3,
            local_iters=2,
        )

    def test_algorithm_initialization(self, algorithm):
        """Test algorithm initialization."""
        assert algorithm.name == "H³-CSP"
        assert algorithm.alphabet == "ACGT"
        assert algorithm.auto_blocks is True
        assert algorithm.k_candidates == 3

    def test_algorithm_run(self, algorithm):
        """Test algorithm execution."""
        result, distance, metadata = algorithm.run()

        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert isinstance(metadata, dict)
        assert distance >= 0
        assert len(result) == len(algorithm.strings[0])
        assert all(c in algorithm.alphabet for c in result)
        assert "blocks_processed" in metadata
        assert "local_iterations" in metadata

    def test_algorithm_manual_blocks(self, small_dataset):
        """Test algorithm with manual block configuration."""
        algorithm = H3CSPAlgorithm(
            small_dataset.sequences,
            "ACGT",
            auto_blocks=False,
            block_size=3,
            k_candidates=4,
            local_iters=3,
            fallback_enabled=True,
        )

        result, distance, metadata = algorithm.run()
        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert "auto_blocks" in metadata
        assert metadata["auto_blocks"] is False


class TestDPCSPAlgorithm:
    """Tests for DP-CSP Algorithm."""

    @pytest.fixture
    def very_small_dataset(self):
        """Very small dataset for DP testing (DP is expensive)."""
        return SyntheticDatasetGenerator.generate_random(
            n=4, length=8, alphabet="ACGT", seed=42
        )

    @pytest.fixture
    def algorithm(self, very_small_dataset):
        """DP-CSP algorithm instance."""
        return DPCSPAlgorithm(
            very_small_dataset.sequences,
            "ACGT",
            max_d=3,
            warn_threshold=10,  # Small threshold for testing
        )

    def test_algorithm_initialization(self, algorithm):
        """Test algorithm initialization."""
        assert algorithm.name == "DP-CSP"
        assert algorithm.alphabet == "ACGT"
        assert algorithm.max_d == 3
        assert algorithm.warn_threshold == 10

    def test_algorithm_run(self, algorithm):
        """Test algorithm execution."""
        result, distance, metadata = algorithm.run()

        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert isinstance(metadata, dict)
        assert distance >= 0
        assert len(result) == len(algorithm.strings[0])
        assert all(c in algorithm.alphabet for c in result)
        assert "max_d_used" in metadata
        assert "optimal" in metadata

    def test_algorithm_with_seed(self, very_small_dataset):
        """Test algorithm with seed for reproducibility."""
        algorithm = DPCSPAlgorithm(
            very_small_dataset.sequences, "ACGT", max_d=2, warn_threshold=20, seed=456
        )

        result, distance, metadata = algorithm.run()
        assert isinstance(result, str)
        assert isinstance(distance, int)
        assert "seed" in metadata
        assert metadata["seed"] == 456

    def test_algorithm_large_dataset_warning(self):
        """Test warning for large datasets."""
        large_dataset = SyntheticDatasetGenerator.generate_random(
            n=15, length=10, alphabet="ACGT", seed=42
        )

        # Should handle large dataset gracefully
        algorithm = DPCSPAlgorithm(
            large_dataset.sequences, "ACGT", max_d=2, warn_threshold=10
        )

        # Should still work but might be slow or issue warnings
        result, distance, metadata = algorithm.run()
        assert isinstance(result, str)
        assert isinstance(distance, int)


class TestAlgorithmCommon:
    """Tests for common algorithm functionality."""

    @pytest.fixture
    def test_dataset(self):
        """Test dataset for common tests."""
        return SyntheticDatasetGenerator.generate_random(
            n=5, length=12, alphabet="ACGT", seed=42
        )

    def test_all_algorithms_basic_functionality(self, test_dataset):
        """Test that all algorithms can run with basic parameters."""
        algorithms = [
            BaselineAlg(test_dataset.sequences, "ACGT"),
            BLFGAAlgorithm(test_dataset.sequences, "ACGT", pop_size=10, max_gens=5),
            CSCAlgorithm(test_dataset.sequences, "ACGT", min_d=1, max_blocks=3),
            H3CSPAlgorithm(
                test_dataset.sequences, "ACGT", k_candidates=2, local_iters=1
            ),
            DPCSPAlgorithm(test_dataset.sequences, "ACGT", max_d=2),
        ]

        for algorithm in algorithms:
            result, distance, metadata = algorithm.run()

            # Common assertions for all algorithms
            assert isinstance(result, str)
            assert isinstance(distance, int)
            assert isinstance(metadata, dict)
            assert distance >= 0
            assert len(result) == len(test_dataset.sequences[0])
            assert all(c in algorithm.alphabet for c in result)
            assert hasattr(algorithm, "name")
            assert hasattr(algorithm, "is_deterministic")

    def test_algorithm_metadata_consistency(self, test_dataset):
        """Test that algorithm metadata is consistent."""
        algorithms = [
            BaselineAlg(test_dataset.sequences, "ACGT"),
            BLFGAAlgorithm(test_dataset.sequences, "ACGT", pop_size=10, max_gens=5),
            CSCAlgorithm(test_dataset.sequences, "ACGT"),
            H3CSPAlgorithm(test_dataset.sequences, "ACGT"),
            DPCSPAlgorithm(test_dataset.sequences, "ACGT", max_d=2),
        ]

        for algorithm in algorithms:
            _, _, metadata = algorithm.run()

            # Common metadata checks
            assert "algorithm" in metadata
            assert metadata["algorithm"] == algorithm.name
            assert "execution_time" in metadata
            assert isinstance(metadata["execution_time"], (int, float))
            assert metadata["execution_time"] >= 0

    def test_algorithm_result_validation(self, test_dataset):
        """Test that algorithm results are valid closest strings."""
        from src.domain.metrics import hamming_distance

        algorithms = [
            BaselineAlg(test_dataset.sequences, "ACGT"),
            BLFGAAlgorithm(test_dataset.sequences, "ACGT", pop_size=10, max_gens=5),
        ]

        for algorithm in algorithms:
            result, reported_distance, _ = algorithm.run()

            # Verify the reported distance is correct
            max_hamming = max(
                hamming_distance(result, seq) for seq in test_dataset.sequences
            )
            assert reported_distance == max_hamming

            # Verify result is in the correct alphabet
            assert all(c in algorithm.alphabet for c in result)

            # Verify result has correct length
            assert len(result) == len(test_dataset.sequences[0])


class TestAlgorithmEdgeCases:
    """Tests for algorithm edge cases."""

    def test_algorithms_with_single_sequence(self):
        """Test algorithms with a single sequence."""
        sequences = ["ACGTACGT"]

        algorithms = [
            BaselineAlg(sequences, "ACGT"),
            BLFGAAlgorithm(sequences, "ACGT", pop_size=5, max_gens=3),
            CSCAlgorithm(sequences, "ACGT"),
            H3CSPAlgorithm(sequences, "ACGT"),
            DPCSPAlgorithm(sequences, "ACGT", max_d=1),
        ]

        for algorithm in algorithms:
            result, distance, metadata = algorithm.run()

            # With single sequence, closest string should be the sequence itself
            assert result == sequences[0]
            assert distance == 0

    def test_algorithms_with_identical_sequences(self):
        """Test algorithms with identical sequences."""
        sequences = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]

        algorithms = [
            BaselineAlg(sequences, "ACGT"),
            BLFGAAlgorithm(sequences, "ACGT", pop_size=5, max_gens=3),
            CSCAlgorithm(sequences, "ACGT"),
            H3CSPAlgorithm(sequences, "ACGT"),
            DPCSPAlgorithm(sequences, "ACGT", max_d=1),
        ]

        for algorithm in algorithms:
            result, distance, metadata = algorithm.run()

            # With identical sequences, closest string should be any of them
            assert result == sequences[0]
            assert distance == 0

    def test_algorithms_with_very_short_sequences(self):
        """Test algorithms with very short sequences."""
        sequences = ["AC", "AT", "AG"]

        algorithms = [
            BaselineAlg(sequences, "ACGT"),
            BLFGAAlgorithm(sequences, "ACGT", pop_size=5, max_gens=3),
            CSCAlgorithm(sequences, "ACGT", min_blocks=1, max_blocks=2),
            H3CSPAlgorithm(sequences, "ACGT", min_block_size=1),
            DPCSPAlgorithm(sequences, "ACGT", max_d=2),
        ]

        for algorithm in algorithms:
            result, distance, metadata = algorithm.run()

            assert len(result) == 2
            assert all(c in "ACGT" for c in result)
            assert distance >= 0
