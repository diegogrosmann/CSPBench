"""
Funções de distância compartilhadas entre algoritmos.

Funções:
    hamming_distance(s1, s2): Calcula a distância de Hamming entre duas strings.
    max_hamming(candidate, strings): Maior distância de Hamming do candidato para o conjunto.
    max_hamming_parallel(candidate, strings): Versão paralela para conjuntos grandes.
"""

import concurrent.futures

def hamming_distance(s1, s2):
    """Calcula a distância de Hamming entre duas strings do mesmo tamanho."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def max_hamming(candidate, strings):
    """Calcula a maior distância de Hamming do candidato para o conjunto."""
    return max(hamming_distance(candidate, s) for s in strings)

def max_hamming_parallel(candidate, strings):
    """
    Calcula a maior distância de Hamming do candidato para o conjunto em paralelo.
    """
    # Executa hamming_distance em paralelo para cada string
    with concurrent.futures.ProcessPoolExecutor() as executor:
        distances = list(executor.map(hamming_distance, [candidate] * len(strings), strings))
    return max(distances)
