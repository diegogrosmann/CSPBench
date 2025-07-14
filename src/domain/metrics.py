"""
Domínio: Métricas para Avaliação de Algoritmos CSP

Este módulo contém funções e classes para cálculo de métricas de distância
e avaliação de soluções no contexto do Closest String Problem.
Implementação pura sem dependências externas.
"""

from typing import List


def hamming_distance(str1: str, str2: str) -> int:
    """
    Calcula distância de Hamming entre duas strings.

    Args:
        str1: Primeira string
        str2: Segunda string

    Returns:
        int: Número de posições onde as strings diferem

    Raises:
        ValueError: Se strings têm comprimentos diferentes
    """
    if len(str1) != len(str2):
        raise ValueError("Strings devem ter mesmo comprimento")

    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


def max_distance(center: str, strings: List[str]) -> int:
    """
    Calcula distância máxima de um centro para conjunto de strings.

    Args:
        center: String central
        strings: Lista de strings de referência

    Returns:
        int: Distância máxima
    """
    if not strings:
        return 0

    return max(hamming_distance(center, s) for s in strings)


def max_hamming(center: str, strings: List[str]) -> int:
    """
    Alias para max_distance para compatibilidade.

    Args:
        center: String central
        strings: Lista de strings de referência

    Returns:
        int: Distância máxima
    """
    return max_distance(center, strings)


def average_distance(center: str, strings: List[str]) -> float:
    """
    Calcula distância média de um centro para conjunto de strings.

    Args:
        center: String central
        strings: Lista de strings de referência

    Returns:
        float: Distância média
    """
    if not strings:
        return 0.0

    total_distance = sum(hamming_distance(center, s) for s in strings)
    return total_distance / len(strings)


def median_distance(center: str, strings: List[str]) -> float:
    """
    Calcula distância mediana de um centro para conjunto de strings.

    Args:
        center: String central
        strings: Lista de strings de referência

    Returns:
        float: Distância mediana
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
    Calcula métrica de diversidade de um conjunto de strings.

    Args:
        strings: Lista de strings

    Returns:
        float: Valor de diversidade (0-1, onde 1 é máxima diversidade)
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
    Calcula força do consenso de uma string central.

    Args:
        center: String central
        strings: Lista de strings de referência

    Returns:
        float: Força do consenso (0-1, onde 1 é consenso perfeito)
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
    Calcula múltiplas métricas de qualidade de uma solução.

    Args:
        center: String central
        strings: Lista de strings de referência

    Returns:
        dict: Dicionário com várias métricas
    """
    return {
        "max_distance": max_distance(center, strings),
        "average_distance": average_distance(center, strings),
        "median_distance": median_distance(center, strings),
        "consensus_strength": consensus_strength(center, strings),
        "diversity": diversity_metric(strings + [center]),
    }


class DistanceCalculator:
    """Calculadora de distâncias com cache para otimização."""

    def __init__(self):
        self._cache = {}

    def hamming_distance_cached(self, str1: str, str2: str) -> int:
        """Calcula distância de Hamming com cache."""
        key = (str1, str2) if str1 <= str2 else (str2, str1)

        if key not in self._cache:
            self._cache[key] = hamming_distance(str1, str2)

        return self._cache[key]

    def clear_cache(self) -> None:
        """Limpa o cache de distâncias."""
        self._cache.clear()

    def cache_size(self) -> int:
        """Retorna tamanho atual do cache."""
        return len(self._cache)


class QualityEvaluator:
    """Avaliador de qualidade de soluções CSP."""

    def __init__(self, strings: List[str]):
        self.strings = strings
        self.calculator = DistanceCalculator()

    def evaluate(self, center: str) -> dict:
        """
        Avalia qualidade de uma string central.

        Args:
            center: String central a ser avaliada

        Returns:
            dict: Métricas de qualidade
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
        Compara duas soluções candidatas.

        Args:
            center1: Primeira string central
            center2: Segunda string central

        Returns:
            dict: Comparação das soluções
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
