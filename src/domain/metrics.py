"""
Domain: Metrics for CSP Algorithm Evaluation

This module contains functions and classes for calculating distance metrics
and evaluating solutions in the context of the Closest String Problem.
Pure implementation without external dependencies.
"""

from typing import List


def hamming_distance(str1: str, str2: str) -> int:
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
        raise ValueError("Strings must have the same length")

    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


def max_distance(center: str, strings: List[str]) -> int:
    """
    Calculate maximum distance from a center to a set of strings.

    Args:
        center: Center string
        strings: List of reference strings

    Returns:
        int: Maximum distance
    """
    if not strings:
        return 0

    return max(hamming_distance(center, s) for s in strings)


def average_distance(center: str, strings: List[str]) -> float:
    """
    Calculate average distance from a center to a set of strings.

    Args:
        center: Center string
        strings: List of reference strings

    Returns:
        float: Average distance
    """
    if not strings:
        return 0.0

    total_distance = sum(hamming_distance(center, s) for s in strings)
    return total_distance / len(strings)


def median_distance(center: str, strings: List[str]) -> float:
    """
    Calculate median distance from a center to a set of strings.

    Args:
        center: Center string
        strings: List of reference strings

    Returns:
        float: Median distance
    """
    if not strings:
        return 0.0

    distances = [hamming_distance(center, s) for s in strings]
    distances.sort()

    n = len(distances)
    if n % 2 == 0:
        return (distances[n // 2 - 1] + distances[n // 2]) / 2.0
    else:
        return float(distances[n // 2])


def diversity_metric(strings: List[str]) -> float:
    """
    Calculate diversity metric of a set of strings.

    Args:
        strings: List of strings

    Returns:
        float: Diversity value (0-1, where 1 is maximum diversity)
    """
    if len(strings) < 2:
        return 0.0

    total_distance = 0
    total_pairs = 0

    for i in range(len(strings)):
        for j in range(i + 1, len(strings)):
            total_distance += hamming_distance(strings[i], strings[j])
            total_pairs += 1

    if total_pairs == 0:
        return 0.0

    avg_distance = total_distance / total_pairs
    max_possible_distance = len(strings[0]) if strings else 0

    return avg_distance / max_possible_distance if max_possible_distance > 0 else 0.0


def consensus_strength(center: str, strings: List[str]) -> float:
    """
    Calculate consensus strength of a center string.

    Args:
        center: Center string
        strings: List of reference strings

    Returns:
        float: Consensus strength (0-1, where 1 is perfect consensus)
    """
    if not strings:
        return 1.0

    max_dist = max_distance(center, strings)
    max_possible_distance = len(center) if center else 0

    if max_possible_distance == 0:
        return 1.0

    return 1.0 - (max_dist / max_possible_distance)


def solution_quality(center: str, strings: List[str]) -> dict:
    """
    Calculate multiple quality metrics for a solution.

    Args:
        center: Center string
        strings: List of reference strings

    Returns:
        dict: Dictionary with various metrics
    """
    return {
        "max_distance": max_distance(center, strings),
        "average_distance": average_distance(center, strings),
        "median_distance": median_distance(center, strings),
        "consensus_strength": consensus_strength(center, strings),
        "diversity": diversity_metric(strings + [center]),
    }


class DistanceCalculator:
    """Distance calculator with cache for optimization."""

    def __init__(self):
        self._cache = {}

    def hamming_distance_cached(self, str1: str, str2: str) -> int:
        """Calculate Hamming distance with cache."""
        key = (str1, str2) if str1 <= str2 else (str2, str1)

        if key not in self._cache:
            self._cache[key] = hamming_distance(str1, str2)

        return self._cache[key]

    def clear_cache(self) -> None:
        """Clear distance cache."""
        self._cache.clear()

    def cache_size(self) -> int:
        """Return current cache size."""
        return len(self._cache)


class QualityEvaluator:
    """CSP solution quality evaluator."""

    def __init__(self, strings: List[str]):
        self.strings = strings
        self.calculator = DistanceCalculator()

    def evaluate(self, center: str) -> dict:
        """
        Evaluate quality of a center string.

        Args:
            center: Center string to be evaluated

        Returns:
            dict: Quality metrics
        """
        distances = [
            self.calculator.hamming_distance_cached(center, s) for s in self.strings
        ]

        return {
            "center": center,
            "max_distance": max(distances) if distances else 0,
            "min_distance": min(distances) if distances else 0,
            "average_distance": sum(distances) / len(distances) if distances else 0,
            "total_distance": sum(distances),
            "num_strings": len(self.strings),
            "string_length": len(center),
            "distances": distances,
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
        eval1 = self.evaluate(center1)
        eval2 = self.evaluate(center2)

        return {
            "center1": eval1,
            "center2": eval2,
            "better_solution": (
                center1 if eval1["max_distance"] < eval2["max_distance"] else center2
            ),
            "max_distance_diff": eval1["max_distance"] - eval2["max_distance"],
            "avg_distance_diff": eval1["average_distance"] - eval2["average_distance"],
        }
