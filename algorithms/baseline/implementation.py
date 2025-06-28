from collections import Counter
from typing import List
import logging
from utils.distance import max_hamming_parallel

logger = logging.getLogger(__name__)

def greedy_consensus(strings: List[str], alphabet: str) -> str:
    """
    Obtém a string consenso (símbolo mais frequente por posição).

    Args:
        strings (List[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.

    Returns:
        str: String consenso.
    """
    logger.debug("Iniciando greedy_consensus")
    L = len(strings[0])
    consensus_chars = []
    for i in range(L):
        freq = Counter(s[i] for s in strings if i < len(s))
        most, count = freq.most_common(1)[0]
        consensus_chars.append(most)
    result = "".join(consensus_chars)
    logger.debug(f"Consenso final: {result}")
    return result

def max_distance(center: str, strings: List[str]) -> int:
    """
    Calcula a maior distância de Hamming do centro para o conjunto.

    Args:
        center (str): String centro.
        strings (List[str]): Lista de strings de entrada.

    Returns:
        int: Maior distância de Hamming.
    """
    # Usa implementação paralela para performance em conjuntos grandes
    return max_hamming_parallel(center, strings)
