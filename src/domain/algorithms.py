"""
Domain: CSP Algorithms

This module contains the core logic for Closest String Problem (CSP) algorithms,
including abstract interfaces, algorithm registry, and genetic operators.
Free from external dependencies following hexagonal architecture.
"""

import random
from abc import ABC, abstractmethod
from typing import Any, Callable, Optional

# =============================================================================
# ALGORITHM REGISTRY
# =============================================================================

global_registry: dict[str, type] = {}


def register_algorithm(cls: type) -> type:
    """
    Decorator to register an algorithm class in the global registry.

    Args:
        cls: Algorithm class to be registered

    Returns:
        type: The class itself, allowing use as decorator
    """
    algorithm_name = getattr(cls, "name", cls.__name__)
    global_registry[algorithm_name] = cls
    return cls


# =============================================================================
# ALGORITHM INTERFACES
# =============================================================================


class CSPAlgorithm(ABC):
    """Abstract base interface for all CSP algorithms."""

    # Required class attributes
    name: str
    default_params: dict
    is_deterministic: bool = False
    supports_internal_parallel: bool = False

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Initialize algorithm with strings and alphabet.

        Args:
            strings: List of dataset strings
            alphabet: Used alphabet
            **params: Algorithm-specific parameters
        """
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}
        self.progress_callback: Optional[Callable[[str, float], None]] = None
        self.warning_callback: Optional[Callable[[str], None]] = None

        # History settings
        self.save_history = params.get("save_history", False)
        self.history_frequency = params.get(
            "history_frequency", 1
        )  # Every N iterations
        self.history = []

    def set_progress_callback(self, callback: Callable[[str, float], None]) -> None:
        """Set callback to report algorithm progress."""
        self.progress_callback = callback

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """Set callback to report algorithm warnings."""
        self.warning_callback = callback

    def _report_progress(self, message: str, progress: float = 0.0) -> None:
        """Report progress if callback is defined."""
        if self.progress_callback:
            self.progress_callback(message, progress)

    def _report_warning(self, message: str) -> None:
        """Report warning if callback is defined."""
        if self.warning_callback:
            self.warning_callback(message)

    def _save_history_entry(self, iteration: int, **data) -> None:
        """
        Save an entry in history if enabled.

        Args:
            iteration: Current iteration number
            **data: Current state data (fitness, best solution, etc.)
        """
        if self.save_history and (iteration % self.history_frequency == 0):
            entry = {"iteration": iteration, "timestamp": self._get_timestamp(), **data}
            self.history.append(entry)

    def _get_timestamp(self) -> float:
        """Return current timestamp for history."""
        import time

        return time.time()

    def get_history(self) -> list[dict]:
        """Return execution history."""
        return self.history.copy()

    def clear_history(self) -> None:
        """Clear current history."""
        self.history.clear()

    @abstractmethod
    def run(self) -> tuple[str, int, dict[str, Any]]:
        """
        Execute algorithm and return structured result.

        Returns:
            tuple: (center_string, max_distance, metadata)
                metadata should contain:
                - history: List of states during execution (if enabled)
                - iterations: Total number of iterations
                - convergence_data: Convergence data (if applicable)
                - other algorithm-specific information
        """

    def set_params(self, **params) -> None:
        """Set new parameters for the algorithm."""
        self.params.update(params)

    def get_metadata(self) -> dict[str, Any]:
        """Return algorithm metadata."""
        return {
            "name": self.name,
            "params": self.params.copy(),
            "is_deterministic": self.is_deterministic,
            "supports_internal_parallel": self.supports_internal_parallel,
            "input_size": len(self.strings),
            "string_length": len(self.strings[0]) if self.strings else 0,
            "alphabet_size": len(self.alphabet),
        }


class Algorithm(ABC):
    """Legacy interface for compatibility with existing code."""

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str):
        """Initialize legacy algorithm."""
        self.strings = strings
        self.alphabet = alphabet

    @abstractmethod
    def solve(self) -> str:
        """Solve CSP problem returning center string."""


# =============================================================================
# GENETIC OPERATORS
# =============================================================================

String = str
Population = list[String]


def mean_hamming_distance(pop: Population) -> float:
    """
    Calculate average Hamming distance between all pairs of individuals.

    Args:
        pop: Population of strings

    Returns:
        float: Average distance between pairs
    """
    if len(pop) < 2:
        return 0.0

    total_distance = 0
    total_pairs = 0

    for i in range(len(pop)):
        for j in range(i + 1, len(pop)):
            distance = sum(c1 != c2 for c1, c2 in zip(pop[i], pop[j]))
            total_distance += distance
            total_pairs += 1

    return total_distance / total_pairs if total_pairs > 0 else 0.0


def mutate_multi(ind: str, alphabet: str, rng: random.Random, n: int = 2) -> str:
    """
    Perform multi-point mutation by changing up to n random positions.

    Args:
        ind: Original string
        alphabet: Set of valid symbols
        rng: Random number generator
        n: Number of mutations to apply

    Returns:
        str: Mutated string
    """
    chars = list(ind)
    L = len(chars)

    for _ in range(n):
        pos = rng.randint(0, L - 1)
        old = chars[pos]

        # Choose symbol different from current
        available = [c for c in alphabet if c != old]
        if available:
            chars[pos] = rng.choice(available)

    return "".join(chars)


def mutate_inversion(ind: str, rng: random.Random) -> str:
    """
    Perform inversion mutation of a random segment.

    Args:
        ind: Original string
        rng: Random number generator

    Returns:
        str: String with inverted segment
    """
    chars = list(ind)
    L = len(chars)

    if L < 2:
        return ind

    # Select two random points
    i, j = sorted(rng.sample(range(L), 2))

    # Invert segment between i and j
    chars[i : j + 1] = chars[i : j + 1][::-1]

    return "".join(chars)


def crossover_one_point(
    parent1: str, parent2: str, rng: random.Random
) -> tuple[str, str]:
    """
    Perform one-point crossover between two parents.

    Args:
        parent1: First parent
        parent2: Second parent
        rng: Random number generator

    Returns:
        tuple: Two resulting children
    """
    L = len(parent1)
    if L <= 1:
        return parent1, parent2

    # Select cut point
    cut_point = rng.randint(1, L - 1)

    # Create children by swapping suffixes
    child1 = parent1[:cut_point] + parent2[cut_point:]
    child2 = parent2[:cut_point] + parent1[cut_point:]

    return child1, child2


def crossover_uniform(
    parent1: str, parent2: str, rng: random.Random, rate: float = 0.5
) -> tuple[str, str]:
    """
    Perform uniform crossover between two parents.

    Args:
        parent1: First parent
        parent2: Second parent
        rng: Random number generator
        rate: Exchange rate per position

    Returns:
        tuple: Two resulting children
    """
    chars1 = list(parent1)
    chars2 = list(parent2)

    for i in range(len(chars1)):
        if rng.random() < rate:
            chars1[i], chars2[i] = chars2[i], chars1[i]

    return "".join(chars1), "".join(chars2)


def refine_greedy(individual: str, strings: list[str]) -> str:
    """
    Greedy local refinement position by position.

    Args:
        individual: String to be refined
        strings: Reference strings

    Returns:
        str: Refined string
    """
    chars = list(individual)
    alphabet = set("".join(strings))

    for pos in range(len(chars)):
        current_char = chars[pos]
        best_char = current_char
        best_distance = _max_distance_to_strings(chars, strings)

        # Test each symbol in alphabet
        for symbol in alphabet:
            if symbol != current_char:
                chars[pos] = symbol
                distance = _max_distance_to_strings(chars, strings)

                if distance < best_distance:
                    best_distance = distance
                    best_char = symbol

        chars[pos] = best_char

    return "".join(chars)


def _max_distance_to_strings(chars: list[str], strings: list[str]) -> int:
    """
    Calculate maximum distance from a candidate to set of strings.

    Args:
        chars: List of characters forming candidate string
        strings: Reference strings

    Returns:
        int: Maximum distance
    """
    candidate = "".join(chars)
    max_dist = 0

    for string in strings:
        dist = sum(c1 != c2 for c1, c2 in zip(candidate, string))
        max_dist = max(max_dist, dist)

    return max_dist


# =============================================================================
# SPECIFIC ALGORITHMS
# =============================================================================

# Specific algorithms have been moved to the algorithms/ folder
# to be dynamically loaded as plug-ins
