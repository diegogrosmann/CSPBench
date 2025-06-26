# baseline.py
"""
baseline.py
===========

Implementa o algoritmo de consenso guloso usado como baseline
(Berman & Karpinski, 2001). Mantém assinatura simples para comparação.
"""

from collections import Counter
from typing import List
import logging

logger = logging.getLogger(__name__)

def greedy_consensus(strings: List[str], alphabet: str) -> str:
    """Obtém a string consenso (símbolo mais frequente por posição)."""
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
