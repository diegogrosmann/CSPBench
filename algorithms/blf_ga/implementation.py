"""
Implementação da heurística BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) para CSP.

Classes:
    BLFGA: Implementa o algoritmo BLF-GA.

Funções auxiliares:
    hamming_dist(a, b): Wrapper para distância de Hamming.
    max_distance(center, strings): Wrapper para distância máxima.
"""

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
from typing import List, Tuple, Callable, Optional
import logging

import numpy as np
from utils.distance import hamming_distance, max_hamming
from .config import BLF_GA_DEFAULTS

logger = logging.getLogger(__name__)

String     = str
Population = List[String]

def hamming_dist(a: String, b: String) -> int:
    """Wrapper para manter compatibilidade."""
    return hamming_distance(a, b)

def max_distance(center: String, strings: List[String]) -> int:
    """Wrapper para manter compatibilidade."""
    return max_hamming(center, strings)

class BLFGA:
    def __init__(
        self,
        strings: List[String],
        alphabet: str,
        pop_size: Optional[int]       = None,
        initial_blocks: Optional[int] = None,
        min_block_len: Optional[int]  = None,
        cross_prob: Optional[float]   = None,
        mut_prob: Optional[float]     = None,
        elite_rate: Optional[float]   = None,
        rediv_freq: Optional[int]     = None,
        max_gens: Optional[int]       = None,
        max_time: Optional[float]     = None,
        seed: Optional[int]           = None,
    ):

        self.strings        = strings
        self.n              = len(strings)
        self.L              = len(strings[0])
        self.alphabet       = alphabet

        # Carrega parâmetros do dicionário de defaults, permite sobrescrever via argumentos
        params = {**BLF_GA_DEFAULTS}
        if pop_size is not None: params['pop_size'] = pop_size
        if initial_blocks is not None: params['initial_blocks'] = initial_blocks
        if min_block_len is not None: params['min_block_len'] = min_block_len
        if cross_prob is not None: params['cross_prob'] = cross_prob
        if mut_prob is not None: params['mut_prob'] = mut_prob
        if elite_rate is not None: params['elite_rate'] = elite_rate
        if rediv_freq is not None: params['rediv_freq'] = rediv_freq
        if max_gens is not None: params['max_gens'] = max_gens
        if max_time is not None: params['max_time'] = max_time
        if seed is not None: params['seed'] = seed

        self.pop_size       = params['pop_size']
        self.initial_blocks = params['initial_blocks']
        self.min_block_len  = params['min_block_len']
        self.cross_prob     = params['cross_prob']
        self.mut_prob       = params['mut_prob']
        self.elite_rate     = params['elite_rate']
        self.rediv_freq     = params['rediv_freq']
        self.max_gens       = params['max_gens']
        self.max_time       = params['max_time']
        self.rng            = random.Random(params['seed'])
        self.progress_callback: Optional[Callable[[str], None]] = None

        self.blocks = self._initial_blocking()

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Define um callback para relatar o progresso."""
        self.progress_callback = callback

    def run(self) -> Tuple[String,int]:

        start = time.time()
        
        if self.progress_callback:
            self.progress_callback("Criando população inicial...")
        
        pop = self._init_population()
        best     = min(pop, key=lambda s: max_distance(s, self.strings))
        best_val = max_distance(best, self.strings)

        for gen in range(1, self.max_gens+1):
            elapsed = time.time() - start
            if elapsed >= self.max_time:
                if self.progress_callback:
                    self.progress_callback(f"Timeout após {elapsed:.1f}s")
                break

            # Progresso via callback (que pode lançar exceção se cancelado)
            if self.progress_callback:
                progress_msg = f"Geração {gen}/{self.max_gens}, melhor={best_val}"
                self.progress_callback(progress_msg)

            repo = self._learn_blocks(pop)
            pop = self._next_generation(pop, repo)
            pop.sort(key=lambda s: max_distance(s, self.strings))

            k = max(1, int(self.elite_rate * self.pop_size))
            pop[:k] = self._refine_elites(pop[:k])

            cur_best = pop[0]
            cur_val  = max_distance(cur_best, self.strings)
            if cur_val < best_val:
                best, best_val = cur_best, cur_val
            
            # Log apenas a cada 50 gerações ou na última
            if gen % 50 == 0 or gen == self.max_gens or best_val == 0:
                logger.debug(f"Geração {gen}: melhor_dist={best_val}")

            if best_val == 0:
                if self.progress_callback:
                    self.progress_callback("Solução ótima encontrada!")
                break
                
            if gen % self.rediv_freq == 0:
                # Redivisão adaptativa (log apenas quando necessário)
                self.blocks = self._adaptive_blocking(pop)

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
