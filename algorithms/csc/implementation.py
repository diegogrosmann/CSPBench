import numpy as np
from sklearn.cluster import DBSCAN
from collections import Counter
import logging
from itertools import combinations, product
from utils.distance import hamming_distance, max_hamming

logger = logging.getLogger(__name__)

# ---------- Funções Auxiliares ----------

def consensus_string(strings):
    """Calcula o consenso por maioria para cada posição."""
    L = len(strings[0])
    consensus = ''
    for i in range(L):
        chars = [s[i] for s in strings]
        most_common = Counter(chars).most_common(1)[0][0]
        consensus += most_common
    return consensus

def split_blocks(s, n_blocks):
    """Divide a string s em n_blocks blocos aproximadamente iguais, cobrindo todo o comprimento."""
    L = len(s)
    block_sizes = [L // n_blocks] * n_blocks
    for i in range(L % n_blocks):
        block_sizes[i] += 1
    blocks = []
    idx = 0
    for size in block_sizes:
        blocks.append(s[idx:idx+size])
        idx += size
    return blocks

def recombine_blocks(blocks_list):
    """Combina blocos de diferentes consensos em uma nova string candidata."""
    return ''.join(blocks_list)

# ---------- Etapa 1: Leitura/conversão do conjunto de strings ----------

def strings_to_array(strings):
    """Converte uma lista de strings em um array numpy (letras -> números)."""
    if not strings:
        raise ValueError("A lista de strings está vazia.")
    L = len(strings[0])
    for idx, s in enumerate(strings):
        if len(s) != L:
            raise ValueError(
                f"Todas as strings devem ter o mesmo comprimento. "
                f"Encontrada string de tamanho {len(s)} na posição {idx}, esperado {L}."
            )
    char_map = {c: i for i, c in enumerate(sorted(set(''.join(strings))))}
    arr = np.array([[char_map[c] for c in s] for s in strings])
    return arr, char_map

# ---------- Parâmetros automáticos ----------

def auto_parameters(strings):
    """
    Define automaticamente:
      - d: raio para DBSCAN (eps)
      - n_blocks: número de blocos para recombinação
    """
    # Calcula todas as distâncias de Hamming entre pares
    distancias = [hamming_distance(s1, s2) for s1, s2 in combinations(strings, 2)]
    media = np.mean(distancias)
    mediana = np.median(distancias)
    minimo = np.min(distancias)
    maximo = np.max(distancias)
    
    # Critério prático para DBSCAN:
    # - Valor um pouco abaixo da média para formar clusters não triviais
    d = max(2, int(np.floor(media * 0.8)))
    
    # Critério prático para n_blocks:
    # - Proporcional ao log(n) e ao tamanho do cluster (min(n, L//25, 4)), mas sempre pelo menos 2
    n = len(strings)
    L = len(strings[0])
    n_blocks = max(2, min(4, n // 6, L // 25))
    
    logger.info(f"[auto_parameters] Hamming média={media:.2f}, min={minimo}, max={maximo} | d(auto)={d} | n_blocks(auto)={n_blocks}")
    return d, n_blocks

# ---------- Etapa 2: Clusterização ----------

def cluster_strings(strings, d, min_samples=2):
    logger.debug(f"Clusterizando strings com d={d}, min_samples={min_samples}")
    arr, char_map = strings_to_array(strings)
    def hamming_metric(x, y):
        return np.sum(x != y)
    clustering = DBSCAN(eps=d, min_samples=min_samples, metric=hamming_metric)
    labels = clustering.fit_predict(arr)
    clusters = {}
    for idx, label in enumerate(labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(strings[idx])
    # Remove outliers (label == -1)
    clusters = {k: v for k, v in clusters.items() if k != -1}
    logger.info(f"Clusters encontrados: {len(clusters)}")
    return list(clusters.values())

# ---------- Etapa 3: Consenso local e recombinação ----------

def heuristic_closest_string(strings, d=None, n_blocks=None):
    if d is None or n_blocks is None:
        d_auto, n_blocks_auto = auto_parameters(strings)
        if d is None:
            d = d_auto
        if n_blocks is None:
            n_blocks = n_blocks_auto

    logger.info(f"Iniciando heuristic_closest_string com d={d}, n_blocks={n_blocks}")
    clusters = cluster_strings(strings, d)
    if not clusters:
        logger.warning("Nenhum cluster encontrado (tente aumentar o parâmetro d)")
        return None

    # Consenso local de cada cluster
    consensos = [consensus_string(cluster) for cluster in clusters]
    logger.debug(f"Consensos locais: {consensos}")

    # Recomposição de blocos
    L = len(strings[0])
    blocks_per_consenso = [split_blocks(cons, n_blocks) for cons in consensos]

    # Monta candidatos recombinando blocos (um bloco de cada consenso)
    candidates = []
    for block_tuple in product(*blocks_per_consenso):
        candidate = recombine_blocks(block_tuple)
        # Garante que o candidato tenha o mesmo tamanho das strings originais
        if len(candidate) > L:
            candidate = candidate[:L]
        elif len(candidate) < L:
            candidate += consensos[0][len(candidate):L]
        candidates.append(candidate)

    logger.info(f"{len(candidates)} candidatos gerados por recombinação de blocos")

    # Avalia os candidatos e escolhe o melhor (menor raio máximo)
    best_candidate = min(candidates, key=lambda cand: max_hamming(cand, strings))
    logger.info(f"Melhor candidato antes da busca local: {best_candidate}")

    # Busca local simples: tenta melhorar cada posição
    best_candidate = local_search(best_candidate, strings)
    logger.info(f"Melhor candidato após busca local: {best_candidate}")

    return best_candidate

def local_search(candidate, strings):
    candidate = list(candidate)
    improved = True
    iterations = 0
    while improved and iterations < 50:  # Limite para evitar loops longos
        improved = False
        iterations += 1
        for i in range(len(candidate)):
            current_char = candidate[i]
            chars = set(s[i] for s in strings)
            for alt in chars:
                if alt == current_char:
                    continue
                new_candidate = candidate.copy()
                new_candidate[i] = alt
                new_candidate_str = ''.join(new_candidate)
                if max_hamming(new_candidate_str, strings) < max_hamming(''.join(candidate), strings):
                    candidate[i] = alt
                    improved = True
    logger.debug(f"Resultado da busca local: {''.join(candidate)}")
    return ''.join(candidate)