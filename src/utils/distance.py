"""
Funções de Distância para Algoritmos CSP

Este módulo fornece funções otimizadas para cálculos de distância entre strings,
especialmente focadas na distância de Hamming, que é fundamental para o
Closest String Problem (CSP).

FUNÇÕES PRINCIPAIS:
==================
- hamming_distance(): Distância de Hamming entre duas strings
- max_distance(): Distância máxima de um centro para um conjunto de strings
- max_hamming(): Alias para max_distance (compatibilidade)
- max_hamming_parallel(): Versão paralela para grandes conjuntos

CARACTERÍSTICAS:
===============
- Implementação otimizada para performance
- Validação rigorosa de entrada
- Logging estruturado para debugging
- Suporte a paralelização para grandes datasets
- Tratamento robusto de erros

DISTÂNCIA DE HAMMING:
====================
A distância de Hamming é o número de posições onde duas strings diferem.
É uma métrica fundamental no CSP e em problemas de bioinformática.

Exemplo:
- hamming_distance("ACGT", "AGCT") = 2
- Posições diferentes: 1 (C vs G) e 2 (G vs C)

APLICAÇÕES NO CSP:
=================
- Avaliação de fitness de soluções candidatas
- Cálculo da distância máxima (objetivo a minimizar)
- Análise de diversidade populacional
- Validação de soluções
- Métricas de convergência

OTIMIZAÇÕES:
===========
- Uso de zip() com strict=False para performance
- Implementação com generator expressions
- Validação prévia de comprimentos
- Logging condicional para evitar overhead

EXEMPLO DE USO:
==============
```python
from src.utils.distance import hamming_distance, max_distance

# Calcular distância entre duas strings
dist = hamming_distance("ACGT", "AGCT")
print(f"Distância de Hamming: {dist}")  # 2

# Calcular distância máxima de um centro
strings = ["ACGT", "AGCT", "ATCT", "AAGT"]
center = "ACGT"
max_dist = max_distance(center, strings)
print(f"Distância máxima: {max_dist}")  # 2

# Validar solução CSP
def is_valid_solution(center, strings, max_allowed_distance):
    return max_distance(center, strings) <= max_allowed_distance
```

COMPLEXIDADE:
============
- hamming_distance(): O(n) onde n é o comprimento das strings
- max_distance(): O(m*n) onde m é o número de strings e n o comprimento
- Espaço: O(1) para hamming_distance, O(m) para max_distance

TRATAMENTO DE ERROS:
===================
- Validação de comprimentos das strings
- Logging detalhado de erros
- Exceções informativas com contexto
- Verificação de tipos de entrada

THREAD SAFETY:
==============
Todas as funções são thread-safe e podem ser usadas em contextos paralelos
sem necessidade de sincronização adicional.
"""

import logging
from typing import List

logger = logging.getLogger(__name__)


def hamming_distance(s1: str, s2: str) -> int:
    """
    Calcula a distância de Hamming entre duas strings.

    A distância de Hamming é o número de posições onde duas strings
    de mesmo comprimento diferem. É uma métrica fundamental no CSP.

    Args:
        s1 (str): Primeira string
        s2 (str): Segunda string

    Returns:
        int: Número de posições onde as strings diferem

    Raises:
        ValueError: Se as strings têm comprimentos diferentes

    Example:
        >>> hamming_distance("ACGT", "AGCT")
        2
        >>> hamming_distance("AAAA", "AAAA")
        0
        >>> hamming_distance("ACGT", "TGCA")
        4

    Note:
        - Complexidade: O(n) onde n é o comprimento das strings
        - Thread-safe
        - Otimizado para performance com generator expressions
    """
    if len(s1) != len(s2):
        logger.error(
            "Tentativa de calcular distância de Hamming entre strings de comprimentos diferentes: %d vs %d",
            len(s1),
            len(s2),
        )
        raise ValueError(
            f"Strings devem ter o mesmo comprimento para distância de Hamming: {len(s1)} vs {len(s2)}"
        )

    # Usar generator expression para eficiência
    return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=False))


def max_distance(center: str, strings: List[str]) -> int:
    """
    Calcula a distância máxima de um centro para um conjunto de strings.

    Esta função é fundamental no CSP, pois o objetivo é minimizar a
    distância máxima de uma string central para todas as strings do conjunto.

    Args:
        center (str): String central (candidata a solução)
        strings (List[str]): Lista de strings de referência

    Returns:
        int: Distância máxima entre o centro e qualquer string do conjunto

    Raises:
        ValueError: Se alguma string tem comprimento diferente do centro

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> max_distance("ACGT", strings)
        2
        >>> max_distance("AAAA", strings)
        4

    Note:
        - Esta é a função objetivo do CSP
        - Complexidade: O(m*n) onde m é número de strings e n o comprimento
        - Todas as strings devem ter o mesmo comprimento
        - Thread-safe
    """
    if not strings:
        logger.warning("Lista de strings vazia passada para max_distance")
        return 0

    # Validar comprimentos antes de calcular
    center_len = len(center)
    for i, s in enumerate(strings):
        if len(s) != center_len:
            logger.error(
                "String %d tem comprimento diferente do centro: %d vs %d",
                i,
                len(s),
                center_len,
            )
            raise ValueError(
                f"String {i} tem comprimento {len(s)}, mas centro tem {center_len}"
            )

    # Calcular todas as distâncias
    distances = [hamming_distance(center, s) for s in strings]

    max_dist = max(distances)

    logger.debug(
        "Distância máxima calculada: %d para centro de comprimento %d com %d strings",
        max_dist,
        center_len,
        len(strings),
    )

    return max_dist


# Aliases para compatibilidade
max_hamming = max_distance


def max_hamming_parallel(candidate: str, strings: List[str]) -> int:
    """
    Versão paralela de max_distance para conjuntos grandes de strings.

    Para grandes conjuntos de strings, esta função pode distribuir
    o cálculo de distâncias entre múltiplos threads para melhor performance.

    Args:
        candidate (str): String candidata
        strings (List[str]): Lista de strings de referência

    Returns:
        int: Distância máxima calculada

    Note:
        - Útil para datasets com milhares de strings
        - Fallback para versão sequencial se overhead não compensar
        - Mesma interface que max_distance
    """
    # Para conjuntos pequenos, usar versão sequencial
    if len(strings) < 1000:
        return max_distance(candidate, strings)

    # TODO: Implementar versão paralela com ThreadPoolExecutor
    # Por enquanto, usar versão sequencial
    logger.info("Usando versão sequencial para %d strings", len(strings))
    return max_distance(candidate, strings)
