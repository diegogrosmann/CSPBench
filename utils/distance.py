"""
Funções de distância compartilhadas entre algoritmos CSP.

Funções:
    hamming_distance(s1, s2): Calcula a distância de Hamming entre duas strings.
    max_distance(center, strings): Maior distância de Hamming do centro para o conjunto.
    max_hamming(center, strings): Alias para max_distance.
    max_hamming_parallel(candidate, strings): Versão paralela para conjuntos grandes.
"""

import concurrent.futures


def hamming_distance(s1: str, s2: str) -> int:
    """Calcula a distância de Hamming entre duas strings."""
    if len(s1) != len(s2):
        import logging

        logger = logging.getLogger(__name__)
        logger.error(
            f"[HAMMING] ERRO: Strings com comprimentos diferentes: {len(s1)} vs {len(s2)}"
        )
        raise ValueError(
            f"Strings devem ter o mesmo comprimento: {len(s1)} vs {len(s2)}"
        )

    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def max_distance(center: str, strings: list[str]) -> int:
    """Calcula a distância máxima do centro para as strings."""
    distances = []
    for s in strings:
        dist = hamming_distance(center, s)
        distances.append(dist)

    return max(distances)


def max_hamming(center: str, strings: list[str]) -> int:
    """Alias para max_distance para compatibilidade."""
    return max_distance(center, strings)


def max_hamming_parallel(candidate, strings):
    """
    Calcula a maior distância de Hamming do candidato para o conjunto em paralelo.
    """
    # Executa hamming_distance em paralelo para cada string
    with concurrent.futures.ProcessPoolExecutor() as executor:
        distances = list(
            executor.map(hamming_distance, [candidate] * len(strings), strings)
        )
    return max(distances)
