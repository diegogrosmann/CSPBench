"""
Simple test to verify BLF-GA basic functionality without full execution.
"""

from unittest.mock import Mock

from src.domain.distance import HammingDistanceCalculator

from ..algorithm import BLFGAAlgorithm
from ..config import BLF_GA_DEFAULTS


class TestBLFGABasicFunctionality:
    """Test basic functionality of BLF-GA without full algorithm execution."""

    def setup_method(self):
        """Set up test fixtures."""
        self.test_strings = ["ACGT", "AGCT", "ATCT"]
        self.alphabet = "ACGT"
        self.distance_calc = HammingDistanceCalculator()
        self.store = Mock()

    def create_algorithm(self, **params) -> BLFGAAlgorithm:
        """Create BLF-GA algorithm instance with test data."""
        return BLFGAAlgorithm(
            strings=self.test_strings,
            alphabet=self.alphabet,
            distance_calculator=self.distance_calc,
            store=self.store,
            **params,
        )

    def test_algorithm_initialization(self):
        """Test basic algorithm initialization."""
        algo = self.create_algorithm(seed=42)

        # Check basic attributes
        assert algo.strings == self.test_strings
        assert algo.alphabet == self.alphabet
        assert algo.name == "BLF-GA"
        assert algo.supports_internal_parallel is True

        # Check that blocks were initialized
        assert hasattr(algo, "blocks")
        assert isinstance(algo.blocks, list)
        assert len(algo.blocks) > 0

        # Check that parameters were normalized
        assert algo.params["pop_size"] >= 20  # Should respect minimum
        assert algo.params["initial_blocks"] >= 1

    def test_parameter_defaults(self):
        """Test that default parameters are applied correctly."""
        algo = self.create_algorithm()

        # Check that all required default parameters are present
        for key, default_value in BLF_GA_DEFAULTS.items():
            assert key in algo.params

    def test_initial_population_generation(self):
        """Test that initial population can be generated."""
        algo = self.create_algorithm(
            seed=42, pop_size=30, min_pop_size=10
        )  # Ensure pop_size > min

        population = algo._build_initial_population()

        assert len(population) == 30
        assert all(len(ind) == len(self.test_strings[0]) for ind in population)
        assert all(isinstance(ind, str) for ind in population)

    def test_population_diversity_calculation(self):
        """Test population diversity calculation."""
        algo = self.create_algorithm(seed=42)

        # Test with diverse population
        diverse_pop = ["AAAA", "TTTT", "GGGG", "CCCC"]
        diversity = algo._population_diversity(diverse_pop)
        assert diversity > 0

        # Test with identical population
        identical_pop = ["AAAA", "AAAA", "AAAA", "AAAA"]
        diversity = algo._population_diversity(identical_pop)
        assert diversity == 0

    def test_blocks_initialization(self):
        """Test block initialization."""
        # Test with specific number of blocks
        algo = self.create_algorithm(initial_blocks=2)
        assert len(algo.blocks) == 2

        # Test that blocks cover entire string
        total_coverage = sum(end - start for start, end in algo.blocks)
        assert total_coverage == len(self.test_strings[0])

        # Test that blocks don't overlap
        blocks = sorted(algo.blocks)
        for i in range(len(blocks) - 1):
            assert blocks[i][1] <= blocks[i + 1][0]

    def test_adaptive_mutation_mechanism(self):
        """Test adaptive mutation mechanism."""
        algo = self.create_algorithm(seed=42)

        # Test that algorithm can run successfully and uses adaptive mutation
        result = algo.run()
        assert result["center_string"] is not None

        # Test that algorithm has adaptive mutation functionality by checking
        # mutation parameters can be accessed and are valid
        assert algo.params["mut_prob"] >= 0.0
        assert algo.params["mut_prob"] <= 1.0

    def test_sort_population(self):
        """Test population sorting by fitness."""
        algo = self.create_algorithm(seed=42)

        population = ["AAAA", "ACGT", "TTTT"]
        sorted_pop = algo._sort_population(population)

        # Check that population is sorted (best first - lower distance)
        assert len(sorted_pop) == len(population)
        for i in range(len(sorted_pop) - 1):
            fitness_i = algo.max_distance(sorted_pop[i])
            fitness_next = algo.max_distance(sorted_pop[i + 1])
            assert fitness_i <= fitness_next

    def test_standardized_distance_usage(self):
        """Test that standardized distance methods are used."""
        algo = self.create_algorithm(seed=42)

        # Test max_distance method
        individual = "ACGT"
        distance = algo.max_distance(individual)
        assert isinstance(distance, (int, float))
        assert distance >= 0

        # Test distance method
        distance = algo.distance("ACGT", "AGCT")
        assert isinstance(distance, (int, float))
        assert distance >= 0

    def test_rng_reproducibility(self):
        """Test random number generator reproducibility."""
        algo1 = self.create_algorithm(seed=42)
        algo2 = self.create_algorithm(seed=42)

        # Generate random values
        pop1 = algo1._build_initial_population()
        pop2 = algo2._build_initial_population()

        # Should be identical with same seed
        assert pop1 == pop2

    def test_parameter_validation(self):
        """Test parameter validation and normalization."""
        # Test various parameter combinations
        test_params = [
            {"pop_size": 50},
            {"pop_size": 1.5, "min_pop_size": 30},
            {"initial_blocks": 3},
            {"initial_blocks": 0.5},
            {"max_gens": 10, "no_improve_patience": 0.2},
            {"crossover_type": "uniform"},
            {"mutation_type": "inversion"},
            {"refinement_type": "swap"},
        ]

        for params in test_params:
            algo = self.create_algorithm(**params)
            # Should not raise any exceptions
            assert algo.params is not None

    def test_crossover_types_availability(self):
        """Test that all crossover types are handled."""
        crossover_types = ["one_point", "uniform", "blend_blocks"]

        for cross_type in crossover_types:
            algo = self.create_algorithm(crossover_type=cross_type, seed=42)
            # Should initialize without errors
            assert algo.params["crossover_type"] == cross_type

    def test_mutation_types_availability(self):
        """Test that all mutation types are handled."""
        mutation_types = ["multi", "inversion", "transposition"]

        for mut_type in mutation_types:
            algo = self.create_algorithm(mutation_type=mut_type, seed=42)
            # Should initialize without errors
            assert algo.params["mutation_type"] == mut_type

    def test_refinement_types_availability(self):
        """Test that all refinement types are handled."""
        refinement_types = ["greedy", "swap", "insertion", "2opt"]

        for ref_type in refinement_types:
            algo = self.create_algorithm(refinement_type=ref_type, seed=42)
            # Should initialize without errors
            assert algo.params["refinement_type"] == ref_type
