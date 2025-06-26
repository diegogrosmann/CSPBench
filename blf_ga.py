# blf_ga.py
"""
blf_ga.py
=========

Implementação em Python 3 da heurística BLF-GA (Blockwise Learning Fusion
+ Genetic Algorithm) para o Closest String Problem (CSP) com refinamento
corrigido para evitar loops infinitos.
"""

import random
import sys
import time
from collections import Counter
from typing import List, Tuple
import logging

import numpy as np

logger = logging.getLogger(__name__)

String     = str
Population = List[String]

def hamming_dist(a: String, b: String) -> int:
    """Distância de Hamming entre strings de mesmo tamanho."""
    d = sum(ch1 != ch2 for ch1, ch2 in zip(a, b))
    return d

def max_distance(center: String, strings: List[String]) -> int:
    """Maior Hamming entre o centro e cada string."""
    vals = [hamming_dist(center, s) for s in strings]
    m = max(vals) if vals else 0
    return m

class BLFGA:
    def __init__(
        self,
        strings: List[String],
        alphabet: str,
        pop_size: int     = 50,
        initial_blocks: int = 5,
        min_block_len: int  = 3,
        cross_prob: float   = 0.9,
        mut_prob: float     = 0.5,
        elite_rate: float   = 0.05,
        rediv_freq: int     = 10,
        max_gens: int       = 300,
        max_time: float     = 60.0,
        seed: int | None    = None,
    ):
        logger.debug("__init__ BLFGA")
        self.strings        = strings
        self.n              = len(strings)
        self.L              = len(strings[0])
        self.alphabet       = alphabet
        self.pop_size       = pop_size
        self.initial_blocks = initial_blocks
        self.min_block_len  = min_block_len
        self.cross_prob     = cross_prob
        self.mut_prob       = mut_prob
        self.elite_rate     = elite_rate
        self.rediv_freq     = rediv_freq
        self.max_gens       = max_gens
        self.max_time       = max_time
        self.rng            = random.Random(seed)

        self.blocks = self._initial_blocking()

    def run(self) -> Tuple[String,int]:
        logger.debug("run() iniciado")
        start = time.time()
        pop = self._init_population()
        best     = min(pop, key=lambda s: max_distance(s, self.strings))
        best_val = max_distance(best, self.strings)
        print(f"Avaliação inicial: melhor dist={best_val}")

        for gen in range(1, self.max_gens+1):
            elapsed = time.time() - start
            if elapsed >= self.max_time:
                print("\nTempo esgotado.")
                break

            print(f"\rGeração {gen}/{self.max_gens} | Melhor dist: {best_val}   ", end="")
            sys.stdout.flush()

            repo = self._learn_blocks(pop)
            pop = self._next_generation(pop, repo)
            pop.sort(key=lambda s: max_distance(s, self.strings))

            k = max(1, int(self.elite_rate * self.pop_size))
            pop[:k] = self._refine_elites(pop[:k])

            cur_best = pop[0]
            cur_val  = max_distance(cur_best, self.strings)
            if cur_val < best_val:
                print(f"\nMelhoria na geração {gen}: {best_val} -> {cur_val}")
                best, best_val = cur_best, cur_val
            
            logger.debug(f"Fim da Geração {gen}: melhor_dist_geral={best_val}, melhor_dist_atual={cur_val}")

            if best_val == 0:
                print("\nSolução ótima encontrada!")
                break
            if gen % self.rediv_freq == 0:
                self.blocks = self._adaptive_blocking(pop)
                logger.debug(f"Geração {gen}: Blocos redivididos para {len(self.blocks)} blocos.")

        return best, best_val

    def _init_population(self) -> Population:
        consensus = "".join(Counter(pos).most_common(1)[0][0] for pos in zip(*self.strings))
        pop = [consensus]
        for _ in range(self.pop_size // 3):
            s = list(consensus)
            for idx, (l, r) in enumerate(self.blocks):
                if self.rng.random() < 0.5:
                    src = self.rng.choice(self.strings)
                    s[l:r] = src[l:r]
            pop.append("".join(s))
        while len(pop) < self.pop_size:
            rand_s = "".join(self.rng.choice(self.alphabet) for _ in range(self.L))
            pop.append(rand_s)
        return pop

    def _initial_blocking(self) -> List[Tuple[int,int]]:
        size   = max(self.min_block_len, self.L // self.initial_blocks)
        blocks = [(i, min(i+size, self.L)) for i in range(0, self.L, size)]
        return blocks

    def _learn_blocks(self, pop: Population) -> List[String]:
        repo = []
        for l, r in self.blocks:
            symbols = [ind[l:r] for ind in pop if len(ind)>=r]
            block   = "".join(Counter(chars).most_common(1)[0][0] for chars in zip(*symbols))
            repo.append(block)
        return repo

    def _tournament_selection(self, pop: Population, k: int = 3) -> String:
        """Seleciona o melhor indivíduo de uma amostra aleatória de tamanho k."""
        tournament_pool = self.rng.sample(pop, k)
        winner = min(tournament_pool, key=lambda s: max_distance(s, self.strings))
        return winner

    def _next_generation(self, pop: Population, repo: List[String]) -> Population:
        pop_sorted = sorted(pop, key=lambda s: max_distance(s, self.strings))
        elite_n     = max(1, int(self.elite_rate * self.pop_size))
        new_pop     = pop_sorted[:elite_n]
        while len(new_pop) < self.pop_size:
            p1 = self._tournament_selection(pop_sorted, k=3)
            p2 = self._tournament_selection(pop_sorted, k=3)
            child  = list(p1)
            if self.rng.random() < self.cross_prob:
                for idx, (l, r) in enumerate(self.blocks):
                    if self.rng.random() < 0.5:
                        # Garante índices dentro dos limites
                        l_adj = min(l, len(child))
                        r_adj = min(r, len(child))
                        if idx < len(repo) and l_adj < r_adj:
                            child[l_adj:r_adj] = repo[idx][:r_adj-l_adj]
                    else:
                        # Garante índices dentro dos limites
                        l_adj = min(l, len(child))
                        r_adj = min(r, len(child))
                        if l_adj < r_adj:
                            r_p2 = min(r_adj, len(p2))
                            if l_adj < r_p2:
                                # Copie apenas a parte disponível de p2
                                segment_length = r_p2 - l_adj
                                child[l_adj:l_adj+segment_length] = list(p2[l_adj:l_adj+segment_length])
            for (l, r) in self.blocks:
                if self.rng.random() < self.mut_prob:
                    # Garante que pos está dentro dos limites de child
                    l_adj = min(l, len(child) - 1)
                    r_adj = min(r, len(child))
                    if l_adj < r_adj:
                        pos = self.rng.randint(l_adj, r_adj - 1)
                        if pos < len(child):
                            alternatives = [c for c in self.alphabet if c != child[pos]]
                            if alternatives:  # Evita lista vazia
                                child[pos] = self.rng.choice(alternatives)
            new_pop.append("".join(child))
        return new_pop

    def _refine_elites(self, pop: Population) -> Population:
        refined = []
        for ind in pop:
            ind_list = list(ind)
            old_val = max_distance(ind, self.strings)
            improved = True
            while improved:
                improved = False
                dists = [hamming_dist("".join(ind_list), s) for s in self.strings]
                worst_str = self.strings[int(np.argmax(dists))]
                for i, (a, b) in enumerate(zip(ind_list, worst_str)):
                    if a != b:
                        candidate = ind_list.copy()
                        candidate[i] = b
                        cand_str = "".join(candidate)
                        new_val = max_distance(cand_str, self.strings)
                        if new_val < old_val:
                            ind_list = candidate
                            old_val = new_val
                            improved = True
                            break
            refined.append("".join(ind_list))
        return refined

    def _adaptive_blocking(self, pop: Population) -> List[Tuple[int,int]]:
        ent = np.zeros(self.L)
        for pos in range(self.L):
            # Considere apenas strings que são longas o suficiente
            valid_chars = [ind[pos] for ind in pop if pos < len(ind)]
            if valid_chars:  # Se temos pelo menos um caractere válido
                cnt = Counter(valid_chars)
                probs = np.array(list(cnt.values()))/len(valid_chars)
                ent[pos] = -np.sum(probs * np.log2(probs))
        thr = 0.7 * ent.max() if ent.max()>0 else 0.0

        blocks = []
        cur    = 0
        while cur < self.L:
            length = self.min_block_len if ent[cur]>thr else self.min_block_len*2
            blocks.append((cur, min(self.L, cur+length)))
            cur += length
        return blocks
