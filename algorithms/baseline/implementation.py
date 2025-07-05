"""
Implementação do algoritmo de consenso ganancioso (Baseline) para CSP.

Funções:
    greedy_consensus(strings, alphabet): Gera string consenso pelo método guloso.
    max_distance(center, strings): Calcula a maior distância de Hamming.
"""

import logging

logger = logging.getLogger(__name__)


def greedy_consensus(strings: list[str], alphabet: str) -> str:
    """
    Constrói uma string consenso usando a estratégia greedy posição por posição.
    Para cada posição, escolhe o símbolo que minimiza a distância máxima.

    Args:
        strings (List[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
    Returns:
        str: String consenso gerada.
    """

    if not strings:
        return ""

    L = len(strings[0])
    consensus = []

    for pos in range(L):
        best_char = None
        best_max_dist = float("inf")

        # Testar cada símbolo do alfabeto
        for char in alphabet:
            # Calcular distância máxima se usarmos este char na posição
            max_dist = 0
            for s in strings:
                dist = sum(1 for i in range(pos + 1) if (consensus + [char])[i] != s[i])
                max_dist = max(max_dist, dist)

            if max_dist < best_max_dist:
                best_max_dist = max_dist
                best_char = char

        consensus.append(best_char)

    result = "".join(consensus)

    # VALIDAÇÃO: Verificar a distância do consenso
    from src.utils.distance import max_distance

    final_distance = max_distance(result, strings)
    logger.info(f"[CONSENSUS] Consenso: {result}, distância: {final_distance}")

    return result


def max_distance(center: str, strings: list[str]) -> int:
    """
    Calcula a distância máxima de Hamming entre o centro e as strings.

    Args:
        center (str): String central.
        strings (List[str]): Lista de strings de entrada.
    Returns:
        int: Maior distância de Hamming encontrada.
    """
    from src.utils.distance import hamming_distance

    distances = [hamming_distance(center, s) for s in strings]
    return max(distances)
