"""
Comprehensive unit tests for BLF-GA genetic operations module.

Tests all genetic operators with parameter validation, reproducibility,
and standardized distance function integration.
"""

import random
import unittest
from unittest.mock import Mock

from algorithms.blf_ga.ops import genetic_ops
from src.domain.distance import HammingDistanceCalculator


class TestGeneticOperations(unittest.TestCase):
    """Test class for all genetic operations in BLF-GA algorithm."""

    def setUp(self):
        """Set up test fixtures."""
        self.alphabet = "ACGT"
        self.test_strings = ["ACGT", "AGCT", "ATCC", "GCGT"]
        self.rng = random.Random(42)
        self.distance_calc = HammingDistanceCalculator(self.test_strings)

    def test_mean_hamming_distance_basic(self):
        """Test basic mean Hamming distance calculation."""
        population = ["ACGT", "AGCT", "ATCT"]

        mean_dist = genetic_ops.mean_hamming_distance(population)
        assert isinstance(mean_dist, float)
        assert mean_dist >= 0.0

    def test_mean_hamming_distance_with_standardized_function(self):
        """Test mean Hamming distance calculation using standardized distance function."""
        population = ["ACGT", "AGCT", "ATCT"]

        # Test with standardized distance function
        mean_dist = genetic_ops.mean_hamming_distance(
            population, distance_func=self.distance_calc.distance
        )
        assert isinstance(mean_dist, float)
        assert mean_dist >= 0.0

        # Test with identical population (should be 0)
        identical_pop = ["ACGT", "ACGT", "ACGT"]
        mean_dist_zero = genetic_ops.mean_hamming_distance(
            identical_pop, distance_func=self.distance_calc.distance
        )
        assert mean_dist_zero == 0.0

    def test_mean_hamming_distance_reproducibility(self):
        """Test that mean Hamming distance is reproducible."""
        population = ["ACGT", "AGCT", "ATCT", "GCGT"]

        # Calculate multiple times - should be identical
        dist1 = genetic_ops.mean_hamming_distance(
            population, distance_func=self.distance_calc.distance
        )
        dist2 = genetic_ops.mean_hamming_distance(
            population, distance_func=self.distance_calc.distance
        )
        assert dist1 == dist2

    def test_mutate_multi_with_rng(self):
        """Test multi-point mutation with controlled RNG."""
        individual = "ACGT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        mutated1 = genetic_ops.mutate_multi(individual, self.alphabet, rng1, 2)

        rng2 = random.Random(42)
        mutated2 = genetic_ops.mutate_multi(individual, self.alphabet, rng2, 2)

        # Results should be identical with same seed
        assert mutated1 == mutated2
        assert len(mutated1) == len(individual)

        # Characters should be from alphabet
        for char in mutated1:
            assert char in self.alphabet

    def test_mutate_multi_different_n_values(self):
        """Test multi-point mutation with different n values."""
        individual = "ACGTACGT"

        for n in [1, 2, 4, 8]:
            mutated = genetic_ops.mutate_multi(individual, self.alphabet, self.rng, n)
            assert len(mutated) == len(individual)

            # Count differences
            differences = sum(1 for a, b in zip(individual, mutated) if a != b)
            # Should not exceed n mutations
            assert differences <= n

    def test_mutate_inversion(self):
        """Test inversion mutation."""
        individual = "ACGTACGT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        mutated1 = genetic_ops.mutate_inversion(individual, rng1)

        rng2 = random.Random(42)
        mutated2 = genetic_ops.mutate_inversion(individual, rng2)

        # Results should be identical with same seed
        assert mutated1 == mutated2
        assert len(mutated1) == len(individual)

        # Characters should be preserved (just reordered)
        assert sorted(mutated1) == sorted(individual)

    def test_mutate_transposition(self):
        """Test transposition mutation."""
        individual = "ACGTACGT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        mutated1 = genetic_ops.mutate_transposition(individual, rng1)

        rng2 = random.Random(42)
        mutated2 = genetic_ops.mutate_transposition(individual, rng2)

        # Results should be identical with same seed
        assert mutated1 == mutated2
        assert len(mutated1) == len(individual)

        # Characters should be preserved (just swapped)
        assert sorted(mutated1) == sorted(individual)

    def test_crossover_one_point(self):
        """Test one-point crossover."""
        parent1 = "AAAA"
        parent2 = "TTTT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        child1_1, child2_1 = genetic_ops.crossover_one_point(parent1, parent2, rng1)

        rng2 = random.Random(42)
        child1_2, child2_2 = genetic_ops.crossover_one_point(parent1, parent2, rng2)

        # Results should be identical with same seed
        assert child1_1 == child1_2
        assert child2_1 == child2_2

        # Length should be preserved
        assert len(child1_1) == len(parent1)
        assert len(child2_1) == len(parent2)

    def test_crossover_uniform(self):
        """Test uniform crossover."""
        parent1 = "AAAA"
        parent2 = "TTTT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        child1_1, child2_1 = genetic_ops.crossover_uniform(parent1, parent2, rng1)

        rng2 = random.Random(42)
        child1_2, child2_2 = genetic_ops.crossover_uniform(parent1, parent2, rng2)

        # Results should be identical with same seed
        assert child1_1 == child1_2
        assert child2_1 == child2_2

        # Length should be preserved
        assert len(child1_1) == len(parent1)
        assert len(child2_1) == len(parent2)

    def test_crossover_blend_blocks(self):
        """Test block blending crossover."""
        parent1 = "AAAACCCC"
        parent2 = "TTTTGGGG"
        blocks = [(0, 4), (4, 8)]  # Two blocks of 4 characters each

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        child1_1, child2_1 = genetic_ops.crossover_blend_blocks(
            parent1, parent2, blocks, rng1
        )

        rng2 = random.Random(42)
        child1_2, child2_2 = genetic_ops.crossover_blend_blocks(
            parent1, parent2, blocks, rng2
        )

        # Results should be identical with same seed
        assert child1_1 == child1_2
        assert child2_1 == child2_2

        # Length should be preserved
        assert len(child1_1) == len(parent1)
        assert len(child2_1) == len(parent2)

    def test_refinement_functions_basic(self):
        """Test all refinement functions basic functionality."""
        individual = "ACGT"
        reference_strings = self.test_strings

        # Test greedy refinement
        refined_greedy = genetic_ops.refine_greedy(individual, reference_strings)
        assert refined_greedy is not None
        assert isinstance(refined_greedy, str)
        assert len(refined_greedy) == len(individual)

        # Test swap refinement
        refined_swap = genetic_ops.refine_swap(individual, reference_strings)
        assert refined_swap is not None
        assert isinstance(refined_swap, str)
        assert len(refined_swap) == len(individual)

        # Test insertion refinement
        refined_insertion = genetic_ops.refine_insertion(individual, reference_strings)
        assert refined_insertion is not None
        assert isinstance(refined_insertion, str)
        assert len(refined_insertion) == len(individual)

        # Test 2-opt refinement
        refined_2opt = genetic_ops.refine_2opt(individual, reference_strings)
        assert refined_2opt is not None
        assert isinstance(refined_2opt, str)
        assert len(refined_2opt) == len(individual)

    def test_refinement_with_max_distance_func(self):
        """Test refinement functions with max distance function."""
        individual = "ACGT"
        reference_strings = self.test_strings
        max_distance_func = self.distance_calc.max_distance

        # Test each refinement function
        refined_greedy = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )
        assert refined_greedy is not None
        assert isinstance(refined_greedy, str)

        refined_swap = genetic_ops.refine_swap(
            individual, reference_strings, max_distance_func
        )
        assert refined_swap is not None
        assert isinstance(refined_swap, str)

        refined_insertion = genetic_ops.refine_insertion(
            individual, reference_strings, max_distance_func
        )
        assert refined_insertion is not None
        assert isinstance(refined_insertion, str)

        refined_2opt = genetic_ops.refine_2opt(
            individual, reference_strings, max_distance_func
        )
        assert refined_2opt is not None
        assert isinstance(refined_2opt, str)

    def test_edge_cases_empty_population(self):
        """Test edge cases with empty or minimal populations."""
        # Test with single individual
        single_pop = ["ACGT"]
        mean_dist = genetic_ops.mean_hamming_distance(
            single_pop, distance_func=self.distance_calc.distance
        )
        assert mean_dist == 0.0

    def test_edge_cases_single_character_strings(self):
        """Test operations with single character strings."""
        individual = "A"
        alphabet = "AT"

        # Test mutation
        mutated = genetic_ops.mutate_multi(individual, alphabet, self.rng, 1)
        assert len(mutated) == 1
        assert mutated in alphabet

        # Test crossover with longer strings
        parent1 = "AA"
        parent2 = "TT"
        child1, child2 = genetic_ops.crossover_one_point(parent1, parent2, self.rng)
        assert len(child1) == 2
        assert len(child2) == 2

    def test_parameter_validation_mutation_n(self):
        """Test parameter validation for mutation operators."""
        individual = "ACGT"

        # Test with n=0 (should not crash)
        mutated = genetic_ops.mutate_multi(individual, self.alphabet, self.rng, 0)
        assert mutated == individual

        # Test with n > length (should handle gracefully)
        mutated = genetic_ops.mutate_multi(individual, self.alphabet, self.rng, 10)
        assert len(mutated) == len(individual)

    def test_refinement_consistency(self):
        """Test that refinement is consistent and deterministic."""
        individual = "ACGT"
        reference_strings = self.test_strings
        max_distance_func = self.distance_calc.max_distance

        # Multiple calls should be deterministic
        refined1 = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )
        refined2 = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )
        assert refined1 == refined2

    def test_fallback_distance_calculation(self):
        """Test that functions work without explicit distance function."""
        population = ["ACGT", "AGCT", "ATCT"]

        # Should work with default implementation
        mean_dist = genetic_ops.mean_hamming_distance(population)
        assert isinstance(mean_dist, float)
        assert mean_dist >= 0.0

    def test_alphabet_consistency(self):
        """Test that mutations respect alphabet constraints."""
        individual = "GGGG"  # Start with alphabet chars
        custom_alphabet = "GC"  # Restricted alphabet

        # Test multiple mutations to ensure alphabet compliance
        for i in range(10):
            mutated = genetic_ops.mutate_multi(individual, custom_alphabet, self.rng, 2)

            # All characters should be from the restricted alphabet
            for char in mutated:
                assert (
                    char in custom_alphabet
                ), f"Character '{char}' not in alphabet '{custom_alphabet}'"

    def test_blocks_parameter_validation(self):
        """Test crossover functions handle blocks parameter correctly."""
        parent1 = "ACGTACGT"
        parent2 = "TTTTTTTT"

        # Test with various block configurations
        block_configs = [
            [(0, 1), (1, 8)],  # One small block, one large
            [(0, 2), (2, 4), (4, 6), (6, 8)],  # Four equal blocks
            [(0, 4), (4, 8)],  # Two equal blocks
            [(0, 8)],  # Single block (whole string)
        ]

        for blocks in block_configs:
            child1, child2 = genetic_ops.crossover_blend_blocks(
                parent1, parent2, blocks, self.rng
            )
            assert len(child1) == len(parent1)
            assert len(child2) == len(parent2)


if __name__ == "__main__":
    unittest.main()
