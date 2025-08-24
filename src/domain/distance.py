"""
Domain: Distance Calculator

Este módulo contém a implementação para cálculo de distâncias
entre strings no contexto do Problema de String Mais Próxima (CSP).
"""

from abc import ABC, abstractmethod


class DistanceCalculator(ABC):
    """
    Classe abstrata para calculadora de distâncias com pool de strings e cache otimizado.
    """

    def __init__(self, strings: list[str] | None = None, use_cache: bool = True):
        """
        Inicializa o calculador de distâncias.

        Args:
            strings: Pool de strings para cálculos
            use_cache: Se deve usar cache para otimizar cálculos repetidos
        """
        self._strings = strings or []
        self.use_cache = use_cache
        self._cache = {} if use_cache else None

    @staticmethod
    @abstractmethod
    def distance(str1: str, str2: str) -> int:
        """
        Calcula distância entre duas strings.

        Args:
            str1: Primeira string
            str2: Segunda string

        Returns:
            int: Distância calculada
        """
        pass

    def max_distance(self, center: str) -> int:
        """
        Calcula distância máxima de um centro para um conjunto de strings.

        Args:
            center: String central

        Returns:
            int: Distância máxima
        """
        if not self._strings:
            return 0

        return max(self.distances_to_all(center))

    def average_distance(self, center: str) -> float:
        """
        Calcula distância média de um centro para um conjunto de strings.

        Args:
            center: String central

        Returns:
            float: Distância média
        """
        if not self._strings:
            return 0.0

        total_distance = self.total_distance(center)
        return total_distance / len(self._strings)

    def total_distance(self, center: str) -> int:
        """
        Calcula distância total de um centro para um conjunto de strings.

        Args:
            center: String central

        Returns:
            int: Distância total
        """
        return sum(self.distances_to_all(center))

    def distances_to_all(self, center: str) -> list[int]:
        """
        Calcula distâncias de um centro para todas as strings.

        Args:
            center: String central

        Returns:
            list[int]: Lista de distâncias
        """
        if not self.use_cache:
            return [self.distance(center, s) for s in self._strings]

        if self._cache.get(center) is None:
            self._cache[center] = [self.distance(center, s) for s in self._strings]

        return self._cache[center]

    def median_distance(self, center: str) -> float:
        """
        Calcula distância mediana de um centro para um conjunto de strings.

        Args:
            center: String central

        Returns:
            float: Distância mediana
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
        Define o pool de strings.

        Args:
            strings: Nova lista de strings
        """
        self._strings = strings.copy()
        # Limpar cache ao alterar strings
        if self._cache is not None:
            self._cache.clear()

    def get_strings(self) -> list[str]:
        """
        Retorna cópia do pool de strings.

        Returns:
            list[str]: Pool de strings atual
        """
        return self._strings.copy()

    def cache_size(self) -> int:
        """
        Retorna o tamanho atual do cache.

        Returns:
            int: Número de entradas no cache
        """
        return len(self._cache) if self._cache is not None else 0

    def clear_cache(self) -> None:
        """Limpa o cache de distâncias."""
        if self._cache is not None:
            self._cache.clear()

    def min_distance(self, center: str) -> int:
        """
        Calcula distância mínima de um centro para um conjunto de strings.

        Args:
            center: String central

        Returns:
            int: Distância mínima
        """
        if not self._strings:
            return 0

        return min(self.distances_to_all(center))

    def diversity_metric(self, strings: list[str] | None = None) -> float:
        """
        Calcula métrica de diversidade de um conjunto de strings.

        Args:
            strings: Lista de strings (usa self._strings se None)

        Returns:
            float: Valor de diversidade (0-1, onde 1 é diversidade máxima)
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
        Calcula força do consenso de uma string central.

        Args:
            center: String central

        Returns:
            float: Força do consenso (0-1, onde 1 é consenso perfeito)
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
        Calcula múltiplas métricas de qualidade para uma solução.

        Args:
            center: String central

        Returns:
            dict: Dicionário com várias métricas
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
        Compara duas soluções candidatas.

        Args:
            center1: Primeira string central
            center2: Segunda string central

        Returns:
            dict: Comparação das soluções
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
# IMPLEMENTAÇÕES CONCRETAS
# =============================================================================


class HammingDistanceCalculator(DistanceCalculator):
    """
    Implementação concreta para cálculo de distância de Hamming.

    A distância de Hamming é o número de posições onde os caracteres diferem
    entre duas strings de mesmo comprimento.
    """

    @staticmethod
    def distance(str1: str, str2: str) -> int:
        """
        Calcula distância de Hamming entre duas strings.

        Args:
            str1: Primeira string
            str2: Segunda string

        Returns:
            int: Número de posições onde as strings diferem

        Raises:
            ValueError: Se as strings têm comprimentos diferentes
        """
        if len(str1) != len(str2):
            raise ValueError(
                f"Strings devem ter o mesmo comprimento: {len(str1)} != {len(str2)}"
            )

        return sum(c1 != c2 for c1, c2 in zip(str1, str2))


class LevenshteinDistanceCalculator(DistanceCalculator):
    """
    Implementação para distância de Levenshtein (distância de edição).

    A distância de Levenshtein é o número mínimo de operações de edição
    (inserção, deleção ou substituição) necessárias para transformar
    uma string em outra.

    Note: Esta implementação é para futura extensibilidade.
    """

    @staticmethod
    def distance(str1: str, str2: str) -> int:
        """
        Calcula distância de Levenshtein entre duas strings.

        Args:
            str1: Primeira string
            str2: Segunda string

        Returns:
            int: Distância de Levenshtein

        Note:
            Implementação básica usando programação dinâmica.
        """
        m, n = len(str1), len(str2)

        # Criar matriz DP
        dp = [[0] * (n + 1) for _ in range(m + 1)]

        # Inicializar primeira linha e coluna
        for i in range(m + 1):
            dp[i][0] = i
        for j in range(n + 1):
            dp[0][j] = j

        # Preencher matriz DP
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if str1[i - 1] == str2[j - 1]:
                    dp[i][j] = dp[i - 1][j - 1]
                else:
                    dp[i][j] = 1 + min(
                        dp[i - 1][j],  # deleção
                        dp[i][j - 1],  # inserção
                        dp[i - 1][j - 1],  # substituição
                    )

        return dp[m][n]


# =============================================================================
# FACTORY PARA CRIAÇÃO DE DISTANCE CALCULATORS
# =============================================================================


def create_distance_calculator(
    distance_method: str, strings: list[str] | None = None, use_cache: bool = True
) -> DistanceCalculator:
    """
    Factory para criar DistanceCalculator baseado no método especificado.

    Args:
        distance_method: Método de cálculo de distância ("hamming" ou "levenshtein")
        strings: Pool de strings para cálculos (opcional)
        use_cache: Se deve usar cache para otimizar cálculos repetidos

    Returns:
        DistanceCalculator: Instância concreta do calculador apropriado

    Raises:
        ValueError: Se o método de distância não for suportado

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
            f"Método de distância não suportado: '{distance_method}'. "
            f"Métodos disponíveis: 'hamming', 'levenshtein'"
        )
