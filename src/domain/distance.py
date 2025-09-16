"""
Domain: Distance Calculator.

This module contains the implementation for calculating distances
between strings in the context of the Closest String Problem (CSP).

Features:
- Abstract distance calculator interface
- String pool management with caching
- Multiple distance metrics (Hamming, Levenshtein)
- Quality metrics and solution comparison
- Factory pattern for easy instantiation

Distance Metrics:
- Hamming Distance: For equal-length strings
- Levenshtein Distance: Edit distance with insertions, deletions, substitutions
"""

from abc import ABC, abstractmethod


class DistanceCalculator(ABC):
    """
    Abstract class for distance calculator with string pool and optimized cache.

    Provides a framework for calculating distances between strings with
    support for string pools, caching, and various quality metrics.

    Features:
    - String pool management
    - Distance caching for performance
    - Quality metrics calculation
    - Solution comparison utilities

    Attributes:
        _strings: Pool of strings for calculations
        use_cache: Whether to use cache for repeated calculations
        _cache: Distance calculation cache
    """

    def __init__(self, strings: list[str] | None = None, use_cache: bool = True):
        """
        Initialize the distance calculator.

        Args:
            strings: Pool of strings for calculations
            use_cache: Whether to use cache to optimize repeated calculations
        """
        self._strings = strings or []
        self.use_cache = use_cache
        self._cache = {} if use_cache else None

    @staticmethod
    @abstractmethod
    def distance(str1: str, str2: str) -> int:
        """
        Calculate distance between two strings.

        Args:
            str1: First string
            str2: Second string

        Returns:
            int: Calculated distance
        """
        pass

    def max_distance(self, center: str) -> int:
        """
        Calculate maximum distance from a center to a set of strings.

        Args:
            center: Center string

        Returns:
            int: Maximum distance
        """
        if not self._strings:
            return 0

        return max(self.distances_to_all(center))

    def average_distance(self, center: str) -> float:
        """
        Calculate average distance from a center to a set of strings.

        Args:
            center: Center string

        Returns:
            float: Average distance
        """
        if not self._strings:
            return 0.0

        total_distance = self.total_distance(center)
        return total_distance / len(self._strings)

    def total_distance(self, center: str) -> int:
        """
        Calculate total distance from a center to a set of strings.

        Args:
            center: Center string

        Returns:
            int: Total distance
        """
        return sum(self.distances_to_all(center))

    def distances_to_all(self, center: str) -> list[int]:
        """
        Calculate distances from a center to all strings.

        Args:
            center: Center string

        Returns:
            list[int]: List of distances
        """
        if not self.use_cache:
            return [self.distance(center, s) for s in self._strings]

        if self._cache.get(center) is None:
            self._cache[center] = [self.distance(center, s) for s in self._strings]

        return self._cache[center]

    def median_distance(self, center: str) -> float:
        """
        Calculate median distance from a center to a set of strings.

        Args:
            center: Center string

        Returns:
            float: Median distance
        """
        if not self._strings:
            return 0.0

        distances = self.distances_to_all(center)
        distances.sort()

        n = len(distances)
        if n % 2 == 0:
            return (distances[n // 2 - 1] + distances[n // 2]) / 2.0
        else:
            return float(distances[n // 2])

    def set_strings(self, strings: list[str]) -> None:
        """
        Set the string pool.

        Args:
            strings: New list of strings
        """
        self._strings = strings.copy()
        # Clear cache when changing strings
        if self._cache is not None:
            self._cache.clear()

    def get_strings(self) -> list[str]:
        """
        Return copy of the string pool.

        Returns:
            list[str]: Current string pool
        """
        return self._strings.copy()

    def cache_size(self) -> int:
        """
        Return the current cache size.

        Returns:
            int: Number of cache entries
        """
        return len(self._cache) if self._cache is not None else 0

    def clear_cache(self) -> None:
        """Clear the distance cache."""
        if self._cache is not None:
            self._cache.clear()

    def min_distance(self, center: str) -> int:
        """
        Calculate minimum distance from a center to a set of strings.

        Args:
            center: Center string

        Returns:
            int: Minimum distance
        """
        if not self._strings:
            return 0

        return min(self.distances_to_all(center))

    def diversity_metric(self, strings: list[str] | None = None) -> float:
        """
        Calculate diversity metric of a set of strings.

        Args:
            strings: List of strings (uses self._strings if None)

        Returns:
            float: Diversity value (0-1, where 1 is maximum diversity)
        """
        target_strings = strings if strings is not None else self._strings

        if len(target_strings) < 2:
            return 0.0

        total_distance = 0
        total_pairs = 0

        for i in range(len(target_strings)):
            for j in range(i + 1, len(target_strings)):
                total_distance += self.distance(target_strings[i], target_strings[j])
                total_pairs += 1

        if total_pairs == 0:
            return 0.0

        avg_distance = total_distance / total_pairs
        max_possible_distance = len(target_strings[0]) if target_strings else 0

        return (
            avg_distance / max_possible_distance if max_possible_distance > 0 else 0.0
        )

    def consensus_strength(self, center: str) -> float:
        """
        Calculate consensus strength of a center string.

        Args:
            center: Center string

        Returns:
            float: Consensus strength (0-1, where 1 is perfect consensus)
        """
        if not self._strings:
            return 1.0

        max_dist = self.max_distance(center)
        max_possible_distance = len(center) if center else 0

        if max_possible_distance == 0:
            return 1.0

        return 1.0 - (max_dist / max_possible_distance)

    def solution_quality(self, center: str) -> dict:
        """
        Calculate multiple quality metrics for a solution.

        Args:
            center: Center string

        Returns:
            dict: Dictionary with various metrics
        """
        return {
            "max_distance": self.max_distance(center),
            "min_distance": self.min_distance(center),
            "average_distance": self.average_distance(center),
            "median_distance": self.median_distance(center),
            "total_distance": self.total_distance(center),
            "consensus_strength": self.consensus_strength(center),
            "diversity": self.diversity_metric(self._strings + [center]),
            "num_strings": len(self._strings),
            "string_length": len(center) if center else 0,
        }

    def compare_solutions(self, center1: str, center2: str) -> dict:
        """
        Compare two candidate solutions.

        Args:
            center1: First center string
            center2: Second center string

        Returns:
            dict: Solution comparison
        """
        eval1 = self.solution_quality(center1)
        eval2 = self.solution_quality(center2)

        return {
            "center1": eval1,
            "center2": eval2,
            "better_solution": (
                center1 if eval1["max_distance"] < eval2["max_distance"] else center2
            ),
            "max_distance_diff": eval1["max_distance"] - eval2["max_distance"],
            "avg_distance_diff": eval1["average_distance"] - eval2["average_distance"],
            "consensus_diff": eval1["consensus_strength"] - eval2["consensus_strength"],
        }


# =============================================================================
# CONCRETE IMPLEMENTATIONS
# =============================================================================


class HammingDistanceCalculator(DistanceCalculator):
    """
    Concrete implementation for Hamming distance calculation.

    The Hamming distance is the number of positions at which the corresponding
    characters are different between two strings of equal length.
    """

    @staticmethod
    def distance(str1: str, str2: str) -> int:
        """
        Calculate Hamming distance between two strings.

        Args:
            str1: First string
            str2: Second string

        Returns:
            int: Number of positions where strings differ

        Raises:
            ValueError: If strings have different lengths
        """
        if len(str1) != len(str2):
            raise ValueError(
                f"Strings must have the same length: {len(str1)} != {len(str2)}"
            )

        return sum(c1 != c2 for c1, c2 in zip(str1, str2))


class LevenshteinDistanceCalculator(DistanceCalculator):
    """
    Implementation for Levenshtein distance (edit distance).

    The Levenshtein distance is the minimum number of edit operations
    (insertion, deletion, or substitution) required to transform
    one string into another.

    Note: This implementation is for future extensibility.
    """

    @staticmethod
    def distance(str1: str, str2: str) -> int:
        """
        Calculate Levenshtein distance between two strings.

        Args:
            str1: First string
            str2: Second string

        Returns:
            int: Levenshtein distance

        Note:
            Basic implementation using dynamic programming.
        """
        m, n = len(str1), len(str2)

        # Create DP matrix
        dp = [[0] * (n + 1) for _ in range(m + 1)]

        # Initialize first row and column
        for i in range(m + 1):
            dp[i][0] = i
        for j in range(n + 1):
            dp[0][j] = j

        # Fill DP matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if str1[i - 1] == str2[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1]
                else:
                    dp[i][j] = 1 + min(
                        dp[i - 1][j],  # deletion
                        dp[i][j - 1],  # insertion
                        dp[i - 1][j - 1],  # substitution
                    )

        return dp[m][n]


# =============================================================================
# FACTORY FOR CREATING DISTANCE CALCULATORS
# =============================================================================


def create_distance_calculator(
    distance_method: str, strings: list[str] | None = None, use_cache: bool = True
) -> DistanceCalculator:
    """
    Factory to create DistanceCalculator based on specified method.

    Args:
        distance_method: Distance calculation method ("hamming" or "levenshtein")
        strings: Pool of strings for calculations (optional)
        use_cache: Whether to use cache to optimize repeated calculations

    Returns:
        DistanceCalculator: Concrete calculator instance of appropriate type

    Raises:
        ValueError: If distance method is not supported

    Examples:
        >>> calc = create_distance_calculator("hamming", ["ACGT", "AGCT"], True)
        >>> distance = calc.distance("ACGT", "AGCT")
        >>> print(distance)  # 1
    """
    distance_method = distance_method.lower().strip()

    if distance_method == "hamming":
        return HammingDistanceCalculator(strings, use_cache)
    elif distance_method == "levenshtein":
        return LevenshteinDistanceCalculator(strings, use_cache)
    else:
        raise ValueError(
            f"Unsupported distance method: '{distance_method}'. "
            f"Available methods: 'hamming', 'levenshtein'"
        )
