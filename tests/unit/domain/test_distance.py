"""
Unit tests for src.domain.distance module.

This module tests distance calculation functionality including:
- HammingDistanceCalculator implementation
- LevenshteinDistanceCalculator implementation
- Distance calculator abstract methods and utilities
- Caching functionality
- Quality metrics calculation
- Solution comparison functionality
- Factory method for creating calculators
"""

from unittest.mock import MagicMock, patch

import pytest

from src.domain.distance import (
    DistanceCalculator,
    HammingDistanceCalculator,
    LevenshteinDistanceCalculator,
    create_distance_calculator,
)


class TestHammingDistanceCalculator:
    """Test HammingDistanceCalculator implementation."""

    def test_init_default_parameters(self):
        """Test initialization with default parameters."""
        calc = HammingDistanceCalculator()

        assert calc._strings == []
        assert calc.use_cache is True
        assert calc._cache == {}

    def test_init_with_strings_and_cache_disabled(self):
        """Test initialization with strings and cache disabled."""
        strings = ["ACGT", "AGCT", "ATCG"]
        calc = HammingDistanceCalculator(strings, use_cache=False)

        assert calc._strings == strings
        assert calc.use_cache is False
        assert calc._cache is None

    def test_distance_equal_strings(self):
        """Test distance calculation for equal strings."""
        calc = HammingDistanceCalculator()
        distance = calc.distance("ACGT", "ACGT")
        assert distance == 0

    def test_distance_different_strings(self):
        """Test distance calculation for different strings."""
        calc = HammingDistanceCalculator()
        distance = calc.distance("ACGT", "AGCT")
        assert distance == 2  # Positions 1 and 2 differ

    def test_distance_completely_different_strings(self):
        """Test distance calculation for completely different strings."""
        calc = HammingDistanceCalculator()
        distance = calc.distance("AAAA", "TTTT")
        assert distance == 4  # All positions differ

    def test_distance_single_character_strings(self):
        """Test distance calculation for single character strings."""
        calc = HammingDistanceCalculator()
        assert calc.distance("A", "A") == 0
        assert calc.distance("A", "T") == 1

    def test_distance_empty_strings(self):
        """Test distance calculation for empty strings."""
        calc = HammingDistanceCalculator()
        distance = calc.distance("", "")
        assert distance == 0

    def test_distance_different_lengths_raises_error(self):
        """Test that different length strings raise ValueError."""
        calc = HammingDistanceCalculator()

        with pytest.raises(ValueError, match="Strings must have the same length"):
            calc.distance("ACG", "ACGT")

    def test_max_distance_with_strings(self):
        """Test max_distance calculation with string pool."""
        strings = ["AAAA", "TTTT", "ACGT"]
        calc = HammingDistanceCalculator(strings)

        max_dist = calc.max_distance("AAAA")
        assert max_dist == 4  # Distance to "TTTT"

    def test_max_distance_empty_strings(self):
        """Test max_distance with empty string pool."""
        calc = HammingDistanceCalculator()
        max_dist = calc.max_distance("ACGT")
        assert max_dist == 0

    def test_min_distance_with_strings(self):
        """Test min_distance calculation with string pool."""
        strings = ["AAAA", "AAAG", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        min_dist = calc.min_distance("AAAA")
        assert min_dist == 0  # Distance to itself

    def test_min_distance_empty_strings(self):
        """Test min_distance with empty string pool."""
        calc = HammingDistanceCalculator()
        min_dist = calc.min_distance("ACGT")
        assert min_dist == 0

    def test_average_distance_with_strings(self):
        """Test average_distance calculation."""
        strings = ["AAAA", "TTTT", "ACGT"]
        calc = HammingDistanceCalculator(strings)

        avg_dist = calc.average_distance("AAAA")
        # Distances: 0 (to AAAA), 4 (to TTTT), 3 (to ACGT)
        expected = (0 + 4 + 3) / 3
        assert avg_dist == expected

    def test_average_distance_empty_strings(self):
        """Test average_distance with empty string pool."""
        calc = HammingDistanceCalculator()
        avg_dist = calc.average_distance("ACGT")
        assert avg_dist == 0.0

    def test_total_distance_with_strings(self):
        """Test total_distance calculation."""
        strings = ["AAAA", "TTTT", "ACGT"]
        calc = HammingDistanceCalculator(strings)

        total_dist = calc.total_distance("AAAA")
        # Distances: 0 (to AAAA), 4 (to TTTT), 3 (to ACGT)
        assert total_dist == 7

    def test_median_distance_odd_count(self):
        """Test median_distance with odd number of strings."""
        strings = ["AAAA", "AAAT", "TTTT"]  # Distances: 0, 1, 4
        calc = HammingDistanceCalculator(strings)

        median_dist = calc.median_distance("AAAA")
        assert median_dist == 1.0  # Middle value

    def test_median_distance_even_count(self):
        """Test median_distance with even number of strings."""
        strings = ["AAAA", "AAAT", "AATT", "TTTT"]  # Distances: 0, 1, 2, 4
        calc = HammingDistanceCalculator(strings)

        median_dist = calc.median_distance("AAAA")
        assert median_dist == 1.5  # Average of middle two values

    def test_median_distance_empty_strings(self):
        """Test median_distance with empty string pool."""
        calc = HammingDistanceCalculator()
        median_dist = calc.median_distance("ACGT")
        assert median_dist == 0.0

    def test_distances_to_all_with_cache(self):
        """Test distances_to_all with caching enabled."""
        strings = ["AAAA", "TTTT", "ACGT"]
        calc = HammingDistanceCalculator(strings, use_cache=True)

        distances = calc.distances_to_all("AAAA")
        assert distances == [0, 4, 3]

        # Verify cache is populated
        assert "AAAA" in calc._cache
        assert calc._cache["AAAA"] == [0, 4, 3]

    def test_distances_to_all_without_cache(self):
        """Test distances_to_all with caching disabled."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings, use_cache=False)

        distances = calc.distances_to_all("AAAA")
        assert distances == [0, 4]

        # Verify no cache is used
        assert calc._cache is None

    def test_distances_to_all_uses_cache_on_second_call(self):
        """Test that distances_to_all uses cache on subsequent calls."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings, use_cache=True)

        # First call - populates cache
        distances1 = calc.distances_to_all("AAAA")

        # Mock the distance method to verify cache is used
        with patch.object(calc, "distance") as mock_distance:
            distances2 = calc.distances_to_all("AAAA")

            # Should return cached result without calling distance method
            assert distances2 == distances1
            mock_distance.assert_not_called()


class TestLevenshteinDistanceCalculator:
    """Test LevenshteinDistanceCalculator implementation."""

    def test_init_default_parameters(self):
        """Test initialization with default parameters."""
        calc = LevenshteinDistanceCalculator()

        assert calc._strings == []
        assert calc.use_cache is True
        assert calc._cache == {}

    def test_distance_equal_strings(self):
        """Test distance calculation for equal strings."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("ACGT", "ACGT")
        assert distance == 0

    def test_distance_single_insertion(self):
        """Test distance calculation for single insertion."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("ACT", "ACGT")
        assert distance == 1  # Insert G

    def test_distance_single_deletion(self):
        """Test distance calculation for single deletion."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("ACGT", "ACT")
        assert distance == 1  # Delete G

    def test_distance_single_substitution(self):
        """Test distance calculation for single substitution."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("ACGT", "ACTT")
        assert distance == 1  # Substitute G with T

    def test_distance_multiple_operations(self):
        """Test distance calculation for multiple operations."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("ACGT", "TGCA")
        assert distance == 4  # Multiple substitutions needed

    def test_distance_empty_to_string(self):
        """Test distance from empty string to non-empty."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("", "ACGT")
        assert distance == 4  # 4 insertions

    def test_distance_string_to_empty(self):
        """Test distance from non-empty string to empty."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("ACGT", "")
        assert distance == 4  # 4 deletions

    def test_distance_empty_strings(self):
        """Test distance between empty strings."""
        calc = LevenshteinDistanceCalculator()
        distance = calc.distance("", "")
        assert distance == 0


class TestDistanceCalculatorStringManagement:
    """Test string pool management functionality."""

    def test_set_strings(self):
        """Test setting string pool."""
        calc = HammingDistanceCalculator()
        strings = ["AAAA", "TTTT", "ACGT"]

        calc.set_strings(strings)

        assert calc._strings == strings
        # Verify it's a copy, not reference
        strings.append("GGGG")
        assert len(calc._strings) == 3

    def test_set_strings_clears_cache(self):
        """Test that setting strings clears cache."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        # Populate cache
        calc.distances_to_all("AAAA")
        assert len(calc._cache) > 0

        # Set new strings should clear cache
        calc.set_strings(["CCCC", "GGGG"])
        assert len(calc._cache) == 0

    def test_get_strings(self):
        """Test getting string pool."""
        strings = ["AAAA", "TTTT", "ACGT"]
        calc = HammingDistanceCalculator(strings)

        retrieved = calc.get_strings()

        assert retrieved == strings
        # Verify it's a copy, not reference
        retrieved.append("GGGG")
        assert len(calc._strings) == 3

    def test_cache_size(self):
        """Test cache size reporting."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        assert calc.cache_size() == 0

        # Populate cache
        calc.distances_to_all("AAAA")
        assert calc.cache_size() == 1

        calc.distances_to_all("TTTT")
        assert calc.cache_size() == 2

    def test_cache_size_without_cache(self):
        """Test cache size when cache is disabled."""
        calc = HammingDistanceCalculator(use_cache=False)
        assert calc.cache_size() == 0

    def test_clear_cache(self):
        """Test cache clearing."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        # Populate cache
        calc.distances_to_all("AAAA")
        assert calc.cache_size() > 0

        # Clear cache
        calc.clear_cache()
        assert calc.cache_size() == 0

    def test_clear_cache_without_cache(self):
        """Test clearing cache when cache is disabled."""
        calc = HammingDistanceCalculator(use_cache=False)
        # Should not raise error
        calc.clear_cache()


class TestDistanceCalculatorQualityMetrics:
    """Test quality metrics and solution comparison functionality."""

    def test_diversity_metric_with_strings(self):
        """Test diversity metric calculation."""
        calc = HammingDistanceCalculator()
        strings = ["AAAA", "TTTT", "ACGT"]

        diversity = calc.diversity_metric(strings)

        # Diversity should be between 0 and 1
        assert 0.0 <= diversity <= 1.0
        assert diversity > 0  # Should be positive for different strings

    def test_diversity_metric_with_identical_strings(self):
        """Test diversity metric with identical strings."""
        calc = HammingDistanceCalculator()
        strings = ["AAAA", "AAAA", "AAAA"]

        diversity = calc.diversity_metric(strings)
        assert diversity == 0.0  # No diversity

    def test_diversity_metric_single_string(self):
        """Test diversity metric with single string."""
        calc = HammingDistanceCalculator()
        diversity = calc.diversity_metric(["AAAA"])
        assert diversity == 0.0

    def test_diversity_metric_empty_strings(self):
        """Test diversity metric with empty list."""
        calc = HammingDistanceCalculator()
        diversity = calc.diversity_metric([])
        assert diversity == 0.0

    def test_diversity_metric_uses_internal_strings(self):
        """Test diversity metric uses internal strings when none provided."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        diversity = calc.diversity_metric()
        assert diversity > 0  # Should use internal strings

    def test_consensus_strength_perfect_consensus(self):
        """Test consensus strength with perfect consensus."""
        strings = ["AAAA", "AAAA", "AAAA"]
        calc = HammingDistanceCalculator(strings)

        strength = calc.consensus_strength("AAAA")
        assert strength == 1.0  # Perfect consensus

    def test_consensus_strength_no_consensus(self):
        """Test consensus strength with no consensus."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        strength = calc.consensus_strength("AAAA")
        assert 0.0 <= strength < 1.0  # Some consensus lost

    def test_consensus_strength_empty_strings(self):
        """Test consensus strength with empty string pool."""
        calc = HammingDistanceCalculator()
        strength = calc.consensus_strength("AAAA")
        assert strength == 1.0  # Default to perfect

    def test_consensus_strength_empty_center(self):
        """Test consensus strength with empty center."""
        strings = [""]  # Use empty strings to match empty center
        calc = HammingDistanceCalculator(strings)
        strength = calc.consensus_strength("")
        assert strength == 1.0

    def test_solution_quality_comprehensive(self):
        """Test comprehensive solution quality metrics."""
        strings = ["AAAA", "AAAT", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        quality = calc.solution_quality("AAAA")

        expected_keys = {
            "max_distance",
            "min_distance",
            "average_distance",
            "median_distance",
            "total_distance",
            "consensus_strength",
            "diversity",
            "num_strings",
            "string_length",
        }

        assert set(quality.keys()) == expected_keys
        assert quality["num_strings"] == 3
        assert quality["string_length"] == 4
        assert quality["min_distance"] == 0
        assert isinstance(quality["consensus_strength"], float)

    def test_compare_solutions_better_solution(self):
        """Test solution comparison functionality."""
        strings = ["AAAA", "TTTT"]
        calc = HammingDistanceCalculator(strings)

        comparison = calc.compare_solutions("AAAA", "GGGG")

        expected_keys = {
            "center1",
            "center2",
            "better_solution",
            "max_distance_diff",
            "avg_distance_diff",
            "consensus_diff",
        }

        assert set(comparison.keys()) == expected_keys
        assert comparison["better_solution"] in ["AAAA", "GGGG"]
        assert isinstance(comparison["max_distance_diff"], (int, float))


class TestDistanceCalculatorFactory:
    """Test factory method for creating distance calculators."""

    def test_create_hamming_calculator(self):
        """Test creating Hamming distance calculator."""
        calc = create_distance_calculator("hamming")
        assert isinstance(calc, HammingDistanceCalculator)

    def test_create_hamming_calculator_uppercase(self):
        """Test creating Hamming calculator with uppercase input."""
        calc = create_distance_calculator("HAMMING")
        assert isinstance(calc, HammingDistanceCalculator)

    def test_create_hamming_calculator_with_spaces(self):
        """Test creating Hamming calculator with whitespace."""
        calc = create_distance_calculator("  hamming  ")
        assert isinstance(calc, HammingDistanceCalculator)

    def test_create_levenshtein_calculator(self):
        """Test creating Levenshtein distance calculator."""
        calc = create_distance_calculator("levenshtein")
        assert isinstance(calc, LevenshteinDistanceCalculator)

    def test_create_levenshtein_calculator_uppercase(self):
        """Test creating Levenshtein calculator with uppercase input."""
        calc = create_distance_calculator("LEVENSHTEIN")
        assert isinstance(calc, LevenshteinDistanceCalculator)

    def test_create_calculator_with_strings_and_cache(self):
        """Test creating calculator with strings and cache settings."""
        strings = ["AAAA", "TTTT"]
        calc = create_distance_calculator("hamming", strings, use_cache=False)

        assert isinstance(calc, HammingDistanceCalculator)
        assert calc._strings == strings
        assert calc.use_cache is False

    def test_create_calculator_unsupported_method(self):
        """Test creating calculator with unsupported method."""
        with pytest.raises(ValueError, match="Unsupported distance method"):
            create_distance_calculator("manhattan")

    def test_create_calculator_empty_method(self):
        """Test creating calculator with empty method string."""
        with pytest.raises(ValueError, match="Unsupported distance method"):
            create_distance_calculator("")


class TestDistanceCalculatorEdgeCases:
    """Test edge cases and error conditions."""

    def test_diversity_metric_zero_max_distance(self):
        """Test diversity metric when max possible distance is zero."""
        calc = HammingDistanceCalculator()
        # Empty strings should have zero max distance
        diversity = calc.diversity_metric(["", ""])
        assert diversity == 0.0

    def test_diversity_metric_total_pairs_zero(self):
        """Test diversity metric edge case with total_pairs calculation."""
        calc = HammingDistanceCalculator()
        # Single string should result in 0 total_pairs
        diversity = calc.diversity_metric(["AAAA"])
        assert diversity == 0.0

    def test_hamming_distance_performance_mock(self):
        """Test that caching improves performance for repeated calls."""
        strings = ["AAAA", "TTTT"] * 100  # Large string pool
        calc = HammingDistanceCalculator(strings, use_cache=True)

        # First call should populate cache
        distances1 = calc.distances_to_all("AAAA")

        # Mock distance method to verify cache usage
        with patch.object(HammingDistanceCalculator, "distance") as mock_distance:
            # Second call should use cache
            distances2 = calc.distances_to_all("AAAA")

            assert distances1 == distances2
            # Distance method should not be called due to caching
            mock_distance.assert_not_called()

    def test_levenshtein_matrix_edge_cases(self):
        """Test Levenshtein distance with matrix edge cases."""
        calc = LevenshteinDistanceCalculator()

        # Test with single character strings
        assert calc.distance("A", "B") == 1
        assert calc.distance("A", "A") == 0

        # Test with one empty string
        assert calc.distance("A", "") == 1
        assert calc.distance("", "A") == 1

    def test_string_pool_isolation(self):
        """Test that string pool changes don't affect other instances."""
        strings1 = ["AAAA", "TTTT"]
        strings2 = ["CCCC", "GGGG"]

        calc1 = HammingDistanceCalculator(strings1)
        calc2 = HammingDistanceCalculator(strings2)

        # Verify isolated string pools
        assert calc1._strings != calc2._strings

        # Change one calculator's strings
        calc1.set_strings(["ATAT"])

        # Other calculator should be unaffected
        assert calc2._strings == strings2
