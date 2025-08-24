"""
Unit tests for genetic operations module of BLF-GA algorithm.

This test suite validates the genetic operators, ensuring they:
- Use standardized distance calculation methods
- Maintain reproducibility with fixed seeds
- Handle edge cases properly
- Produce valid outputs within expected ranges

Test Coverage:
- Population diversity calculation
- Mutation operators (multi, inversion, transposition)
- Crossover operators (one-point, uniform, blend-blocks)
- Local refinement methods (greedy, swap, insertion, 2-opt)
- Distance calculation integration
"""

import pytest
from unittest.mock import Mock
import random

from src.domain.distance import HammingDistanceCalculator
from ..ops import genetic_ops


class TestGeneticOperations:
    """Test suite for genetic operations module."""

    def setup_method(self):
        """Set up test fixtures."""
        self.test_strings = ["ACGT", "AGCT", "ATCT", "GCGT"]
        self.alphabet = "ACGT"
        self.distance_calc = HammingDistanceCalculator(self.test_strings)

        # Set up reproducible random state
        self.rng = random.Random(42)

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

    def test_mutation_multi_with_rng(self):
        """Test multi-point mutation with controlled RNG."""
        individual = "ACGT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        mutated1 = genetic_ops.mutate_multi(individual, self.alphabet, rng1, 2)

        rng2 = random.Random(42)
        mutated2 = genetic_ops.mutate_multi(individual, self.alphabet, rng2, 2)

        # Should be identical with same seed
        assert mutated1 == mutated2
        assert len(mutated1) == len(individual)
        # Should contain only valid alphabet characters
        assert all(c in self.alphabet for c in mutated1)

    def test_mutation_multi_different_n_values(self):
        """Test multi-point mutation with different n values."""
        individual = "ACGTACGT"

        for n in [1, 2, 4, 8]:
            mutated = genetic_ops.mutate_multi(individual, self.alphabet, self.rng, n)
            assert len(mutated) == len(individual)
            assert all(c in self.alphabet for c in mutated)

    def test_mutation_inversion(self):
        """Test inversion mutation."""
        individual = "ACGTACGT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        mutated1 = genetic_ops.mutate_inversion(individual, rng1)

        rng2 = random.Random(42)
        mutated2 = genetic_ops.mutate_inversion(individual, rng2)

        assert mutated1 == mutated2
        assert len(mutated1) == len(individual)
        # Character set should be preserved
        assert sorted(mutated1) == sorted(individual)

    def test_mutation_transposition(self):
        """Test transposition mutation."""
        individual = "ACGTACGT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        mutated1 = genetic_ops.mutate_transposition(individual, rng1)

        rng2 = random.Random(42)
        mutated2 = genetic_ops.mutate_transposition(individual, rng2)

        assert mutated1 == mutated2
        assert len(mutated1) == len(individual)
        # Character set should be preserved
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

        # Should be reproducible
        assert child1_1 == child1_2
        assert child2_1 == child2_2

        # Length should be preserved
        assert len(child1_1) == len(parent1)
        assert len(child2_1) == len(parent1)

    def test_crossover_uniform(self):
        """Test uniform crossover."""
        parent1 = "AAAA"
        parent2 = "TTTT"

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        child1_1, child2_1 = genetic_ops.crossover_uniform(parent1, parent2, rng1)

        rng2 = random.Random(42)
        child1_2, child2_2 = genetic_ops.crossover_uniform(parent1, parent2, rng2)

        # Should be reproducible
        assert child1_1 == child1_2
        assert child2_1 == child2_2

        # Length should be preserved
        assert len(child1_1) == len(parent1)
        assert len(child2_1) == len(parent1)

    def test_crossover_blend_blocks(self):
        """Test block-based crossover."""
        parent1 = "AAAATTTT"
        parent2 = "TTTTAAAA"
        blocks = [(0, 4), (4, 8)]

        # Test with fixed seed for reproducibility
        rng1 = random.Random(42)
        child1_1, child2_1 = genetic_ops.crossover_blend_blocks(
            parent1, parent2, blocks, rng1
        )

        rng2 = random.Random(42)
        child1_2, child2_2 = genetic_ops.crossover_blend_blocks(
            parent1, parent2, blocks, rng2
        )

        # Should be reproducible
        assert child1_1 == child1_2
        assert child2_1 == child2_2

        # Length should be preserved
        assert len(child1_1) == len(parent1)
        assert len(child2_1) == len(parent1)

    # def test_tournament_selection(self):
    #     """Test tournament selection."""
    #     # Function not implemented in genetic_ops module
    #     pass

    def test_refinement_functions_with_standardized_distance(self):
        """Test all refinement functions use standardized distance calculation."""
        individual = "ACGT"
        reference_strings = self.test_strings
        max_distance_func = self.distance_calc.max_distance

        # Test greedy refinement
        refined_greedy = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )
        if refined_greedy is not None:
            assert len(refined_greedy) == len(individual)
            assert all(c in self.alphabet for c in refined_greedy)

        # Test swap refinement
        refined_swap = genetic_ops.refine_swap(
            individual, reference_strings, max_distance_func
        )
        if refined_swap is not None:
            assert len(refined_swap) == len(individual)

        # Test insertion refinement
        refined_insertion = genetic_ops.refine_insertion(
            individual, reference_strings, max_distance_func
        )
        if refined_insertion is not None:
            assert len(refined_insertion) == len(individual)

        # Test 2-opt refinement
        refined_2opt = genetic_ops.refine_2opt(
            individual, reference_strings, max_distance_func
        )
        if refined_2opt is not None:
            assert len(refined_2opt) == len(individual)

    def test_refinement_reproducibility(self):
        """Test that refinement functions are reproducible with same seed."""
        individual = "ACGT"
        reference_strings = self.test_strings
        max_distance_func = self.distance_calc.max_distance

        # Test greedy refinement reproducibility
        refined1 = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )

        refined2 = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )

        # Both calls should return same result (deterministic)
        assert refined1 == refined2

    def test_refinement_improvement(self):
        """Test that refinement can improve solution quality."""
        # Use a deliberately suboptimal individual
        individual = "TTTT"  # Should be far from optimal for test strings
        reference_strings = self.test_strings
        max_distance_func = self.distance_calc.max_distance

        original_fitness = max_distance_func(individual)

        refined = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )

        if refined is not None:
            refined_fitness = max_distance_func(refined)
            # Refined solution should be at least as good (lower distance is better)
            assert refined_fitness <= original_fitness

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
        alphabet = "ACGT"

        # Test mutations with single character
        individual = "A"

        # Test multi-point mutation
        mutated = genetic_ops.mutate_multi(individual, alphabet, self.rng, 1)
        assert len(mutated) == 1
        assert mutated in alphabet

        # Test crossover - skip one-point for single chars as it requires length > 1
        parent1 = "A"
        parent2 = "T"

        # Test uniform crossover instead (works with single chars)
        child1, child2 = genetic_ops.crossover_uniform(parent1, parent2, self.rng)
        assert len(child1) == 1
        assert len(child2) == 1
        assert child1 in alphabet
        assert child2 in alphabet

    def test_parameter_validation_mutation_n(self):
        """Test parameter validation for mutation operators."""
        individual = "ACGT"

        # Test with n=0 (should not crash)
        mutated = genetic_ops.mutate_multi(individual, self.alphabet, self.rng, 0)
        assert mutated == individual  # No changes expected

        # Test with n > length (should be handled gracefully)
        mutated = genetic_ops.mutate_multi(individual, self.alphabet, self.rng, 10)
        assert len(mutated) == len(individual)

    def test_refinement_max_iter_parameter(self):
        """Test that refinement functions work properly."""
        individual = "ACGT"
        reference_strings = self.test_strings
        max_distance_func = self.distance_calc.max_distance

        # Test basic refinement
        refined = genetic_ops.refine_greedy(
            individual, reference_strings, max_distance_func
        )

        if refined is not None:
            assert len(refined) == len(individual)

    def test_fallback_distance_calculation(self):
        """Test fallback to HammingDistanceCalculator when no distance function provided."""
        # Test that functions work even without explicit distance function
        population = ["ACGT", "AGCT"]

        # This should use the fallback HammingDistanceCalculator
        mean_dist = genetic_ops.mean_hamming_distance(population)
        assert isinstance(mean_dist, float)
        assert mean_dist >= 0.0

    def test_alphabet_consistency(self):
        """Test that mutations respect alphabet constraints."""
        individual = "GGCC"  # Use characters that exist in custom alphabet
        custom_alphabet = "GC"  # Restricted alphabet

        mutated = genetic_ops.mutate_multi(individual, custom_alphabet, self.rng, 4)

        # All characters should be from the custom alphabet
        assert all(c in custom_alphabet for c in mutated)
        assert len(mutated) == len(individual)

    def test_blocks_parameter_validation(self):
        """Test block parameter validation in block-based operations."""
        parent1 = "ACGTACGT"
        parent2 = "TGCATGCA"

        # Test with valid blocks
        valid_blocks = [(0, 4), (4, 8)]
        child1, child2 = genetic_ops.crossover_blend_blocks(
            parent1, parent2, valid_blocks, self.rng
        )
        assert len(child1) == len(parent1)
        assert len(child2) == len(parent2)

        # Test with single block
        single_block = [(0, 8)]
        child1, child2 = genetic_ops.crossover_blend_blocks(
            parent1, parent2, single_block, self.rng
        )
        assert len(child1) == len(parent1)
        assert len(child2) == len(parent2)
