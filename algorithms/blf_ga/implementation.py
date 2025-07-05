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

import logging
import random
import time
from collections import Counter
from collections.abc import Callable

import numpy as np

from src.utils.distance import hamming_distance, max_distance

from .config import BLF_GA_DEFAULTS
from .ops import genetic_ops

logger = logging.getLogger(__name__)

String = str
Population = list[String]


def hamming_dist(a: String, b: String) -> int:
    """Wrapper para manter compatibilidade."""
    return hamming_distance(a, b)


class BLFGA:
    def __init__(
        self,
        strings: list[String],
        alphabet: str,
        pop_size: int | None = None,
        initial_blocks: int | None = None,
        min_block_len: int | None = None,
        cross_prob: float | None = None,
        mut_prob: float | None = None,
        elite_rate: float | None = None,
        rediv_freq: int | None = None,
        max_gens: int | None = None,
        max_time: float | None = None,
        seed: int | None = None,
        # Novos parâmetros
        immigrant_freq: int | None = None,
        immigrant_ratio: float | None = None,
        diversity_threshold: float | None = None,
        mutation_adapt_N: int | None = None,
        mutation_adapt_factor: float | None = None,
        mutation_adapt_duration: int | None = None,
        mutation_type: str | None = None,
        mutation_multi_n: int | None = None,
        tournament_k: int | None = None,
        crossover_type: str | None = None,
        niching: bool | None = None,
        niching_radius: int | None = None,
        refinement_type: str | None = None,
        refine_elites: str | None = None,
        refine_iter_limit: int | None = None,
        restart_patience: int | None = None,
        restart_ratio: float | None = None,
        disable_elitism_gens: int | None = None,
        no_improve_patience: int | None = None,
    ):
        self.strings = strings
        self.n = len(strings)
        self.L = len(strings[0])
        self.alphabet = alphabet

        # Carrega parâmetros do dicionário de defaults, permite sobrescrever via argumentos
        params = {**BLF_GA_DEFAULTS}
        if pop_size is not None:
            params["pop_size"] = pop_size
        if initial_blocks is not None:
            params["initial_blocks"] = initial_blocks
        if min_block_len is not None:
            params["min_block_len"] = min_block_len
        if cross_prob is not None:
            params["cross_prob"] = cross_prob
        if mut_prob is not None:
            params["mut_prob"] = mut_prob
        if elite_rate is not None:
            params["elite_rate"] = elite_rate
        if rediv_freq is not None:
            params["rediv_freq"] = rediv_freq
        if max_gens is not None:
            params["max_gens"] = max_gens
        if max_time is not None:
            params["max_time"] = max_time
        if seed is not None:
            params["seed"] = seed

        # Novos parâmetros
        self.immigrant_freq = params["immigrant_freq"]
        self.immigrant_ratio = params["immigrant_ratio"]
        self.diversity_threshold = params["diversity_threshold"]
        self.mutation_adapt_N = params["mutation_adapt_N"]
        self.mutation_adapt_factor = params["mutation_adapt_factor"]
        self.mutation_adapt_duration = params["mutation_adapt_duration"]
        self.mutation_type = params["mutation_type"]
        self.mutation_multi_n = params["mutation_multi_n"]
        self.tournament_k = params["tournament_k"]
        self.crossover_type = params["crossover_type"]
        self.niching = params["niching"]
        self.niching_radius = params["niching_radius"]
        self.refinement_type = params["refinement_type"]
        self.refine_elites = params["refine_elites"]
        self.refine_iter_limit = params["refine_iter_limit"]
        self.restart_patience = params["restart_patience"]
        self.restart_ratio = params["restart_ratio"]
        self.disable_elitism_gens = params["disable_elitism_gens"]

        # Conversão dinâmica de pop_size e initial_blocks se forem proporções
        self.pop_size = self._resolve_dynamic_param(params["pop_size"], self.n, mode="multiplier")
        self.initial_blocks = self._resolve_dynamic_param(params["initial_blocks"], self.L, mode="proportion")
        self.min_block_len = params["min_block_len"]
        self.cross_prob = params["cross_prob"]
        self.mut_prob = params["mut_prob"]
        self.elite_rate = params["elite_rate"]
        self.rediv_freq = params["rediv_freq"]
        self.max_gens = params["max_gens"]
        self.max_time = params["max_time"]
        self.rng = random.Random(params["seed"])
        self.progress_callback: Callable[[str], None] | None = None

        # Inicializa os blocos após todos os parâmetros necessários
        self.blocks = self._initial_blocking()
        self.history = []  # Histórico de distâncias por geração

        # Early stopping: converte no_improve_patience para int se for proporção
        patience_param = params.get("no_improve_patience", 0)
        if isinstance(patience_param, float) and 0 < patience_param < 1:
            self.no_improve_patience = max(1, int(patience_param * self.max_gens))
        else:
            self.no_improve_patience = int(patience_param)

    @staticmethod
    def _resolve_dynamic_param(param, ref_value, mode="proportion"):
        """
        Para pop_size (mode='multiplier'):
            - Se param for float >= 1, retorna int(param * ref_value)
            - Se param for int, retorna param
        Para initial_blocks (mode='proportion'):
            - Se param for float entre 0 e 1, retorna int(param * ref_value)
            - Se param for int, retorna param
        """
        if mode == "multiplier":
            if isinstance(param, float) and param >= 1:
                return max(1, int(param * ref_value))
            return int(param)
        else:
            if isinstance(param, float) and 0 < param <= 1:
                return max(1, int(param * ref_value))
            return int(param)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Define um callback para relatar o progresso."""
        self.progress_callback = callback

    def run(self) -> tuple[String, int, list]:
        start = time.time()

        if self.progress_callback:
            self.progress_callback("Criando população inicial...")

        pop = self._init_population()
        best = min(pop, key=lambda s: max_distance(s, self.strings))
        best_val = max_distance(best, self.strings)
        self.history = [best_val]
        no_improve = 0
        mut_prob_backup = self.mut_prob
        mut_adapt_timer = 0
        for gen in range(1, self.max_gens + 1):
            elapsed = time.time() - start
            if elapsed >= self.max_time:
                if self.progress_callback:
                    self.progress_callback(f"Timeout após {elapsed:.1f}s")
                break

            # Progresso via callback (que pode lançar exceção se cancelado)
            if self.progress_callback:
                progress_msg = f"Geração {gen}/{self.max_gens}, melhor={best_val}"
                self.progress_callback(progress_msg)

            # --- Imigrantes aleatórios ---
            if self.immigrant_freq and gen % self.immigrant_freq == 0:
                n_imm = int(self.immigrant_ratio * self.pop_size)
                for _ in range(n_imm):
                    rand_s = "".join(self.rng.choice(self.alphabet) for _ in range(self.L))
                    pop[-(_ + 1)] = rand_s

            # --- Diversidade e mutação adaptativa ---
            diversity = genetic_ops.mean_hamming_distance(pop)
            if diversity < self.diversity_threshold * self.L:
                self.mut_prob = mut_prob_backup * self.mutation_adapt_factor
                mut_adapt_timer = self.mutation_adapt_duration
            # Convergência de fitness
            if len(self.history) > self.mutation_adapt_N and all(
                self.history[-i] == self.history[-1] for i in range(1, self.mutation_adapt_N + 1)
            ):
                self.mut_prob = mut_prob_backup * self.mutation_adapt_factor
                mut_adapt_timer = self.mutation_adapt_duration
            if mut_adapt_timer > 0:
                mut_adapt_timer -= 1
                if mut_adapt_timer == 0:
                    self.mut_prob = mut_prob_backup

            repo = self._learn_blocks(pop)
            pop = self._next_generation(pop, repo)
            pop.sort(key=lambda s: max_distance(s, self.strings))

            # --- Elitismo adaptativo ---
            k = max(1, int(self.elite_rate * self.pop_size))
            if self.disable_elitism_gens and (gen % self.disable_elitism_gens == 0):
                elites = []
            else:
                elites = pop[:k]
            if self.refine_elites == "all":
                pop[:k] = self._refine_elites(elites)
            else:
                if elites:
                    pop[0] = self._refine_elites([elites[0]])[0]

            cur_best = pop[0]
            cur_val = max_distance(cur_best, self.strings)
            self.history.append(cur_val)
            if cur_val < best_val:
                best, best_val = cur_best, cur_val
                no_improve = 0
            else:
                no_improve += 1

            # --- Restart ---
            if self.restart_patience and no_improve >= self.restart_patience:
                n_restart = int(self.restart_ratio * self.pop_size)
                for i in range(n_restart):
                    pop[-(i + 1)] = "".join(self.rng.choice(self.alphabet) for _ in range(self.L))
                no_improve = 0

            # --- Early stopping: encerra se não houver melhoria por X gerações ---
            if self.no_improve_patience and no_improve >= self.no_improve_patience:
                if self.progress_callback:
                    self.progress_callback(f"Encerrando por early stopping após {no_improve} gerações sem melhoria.")
                break

            # Log apenas a cada 50 gerações ou na última
            if gen % 50 == 0 or gen == self.max_gens or best_val == 0:
                logger.debug(f"Geração {gen}: melhor_dist={best_val}, diversidade={diversity:.2f}")

            if best_val == 0:
                if self.progress_callback:
                    self.progress_callback("Solução ótima encontrada!")
                break

            if gen % self.rediv_freq == 0:
                # Redivisão adaptativa (log apenas quando necessário)
                self.blocks = self._adaptive_blocking(pop)

        return best, best_val, self.history

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

    def _initial_blocking(self) -> list[tuple[int, int]]:
        size = max(self.min_block_len, self.L // self.initial_blocks)
        blocks = [(i, min(i + size, self.L)) for i in range(0, self.L, size)]
        return blocks

    def _learn_blocks(self, pop: Population) -> list[String]:
        repo = []
        for l, r in self.blocks:
            symbols = [ind[l:r] for ind in pop if len(ind) >= r]
            block = "".join(Counter(chars).most_common(1)[0][0] for chars in zip(*symbols))
            repo.append(block)
        return repo

    def _tournament_selection(self, pop: Population, k: int = 2) -> String:
        """Seleciona o melhor indivíduo de uma amostra aleatória de tamanho k."""
        if k is None:
            k = self.tournament_k if self.tournament_k is not None else 2
        tournament_pool = self.rng.sample(pop, k)
        winner = min(tournament_pool, key=lambda s: max_distance(s, self.strings))
        return winner

    def _apply_crossover(self, p1, p2):
        if self.crossover_type == "one_point":
            c1, c2 = genetic_ops.crossover_one_point(p1, p2, self.rng)
        elif self.crossover_type == "uniform":
            c1, c2 = genetic_ops.crossover_uniform(p1, p2, self.rng)
        elif self.crossover_type == "blend_blocks":
            c1, c2 = genetic_ops.crossover_blend_blocks(p1, p2, self.blocks, self.rng)
        else:
            c1, c2 = p1, p2
        return c1, c2

    def _apply_mutation(self, ind):
        if self.mutation_type == "multi":
            return genetic_ops.mutate_multi(ind, self.alphabet, self.rng, n=self.mutation_multi_n)
        elif self.mutation_type == "inversion":
            return genetic_ops.mutate_inversion(ind, self.rng)
        elif self.mutation_type == "transposition":
            return genetic_ops.mutate_transposition(ind, self.rng)
        else:
            return ind

    def _apply_refinement(self, ind):
        if self.refinement_type == "swap":
            return genetic_ops.refine_swap(ind, self.strings)
        elif self.refinement_type == "insertion":
            return genetic_ops.refine_insertion(ind, self.strings)
        elif self.refinement_type == "2opt":
            return genetic_ops.refine_2opt(ind, self.strings)
        else:
            return ind

    def _next_generation(self, pop: Population, repo: list[String]) -> Population:
        pop_sorted = sorted(pop, key=lambda s: max_distance(s, self.strings))
        elite_n = max(1, int(self.elite_rate * self.pop_size))
        new_pop = pop_sorted[:elite_n]
        while len(new_pop) < self.pop_size:
            p1 = self._tournament_selection(pop_sorted, k=self.tournament_k)
            p2 = self._tournament_selection(pop_sorted, k=self.tournament_k)
            if self.rng.random() < self.cross_prob:
                c1, c2 = self._apply_crossover(p1, p2)
            else:
                c1, c2 = p1, p2
            c1 = self._apply_mutation(c1)
            c2 = self._apply_mutation(c2)
            new_pop.append(c1)
            if len(new_pop) < self.pop_size:
                new_pop.append(c2)
        return new_pop[: self.pop_size]

    def _refine_elites(self, pop: Population) -> Population:
        # Paralelização opcional pode ser feita aqui
        from concurrent.futures import ThreadPoolExecutor

        def refine(ind):
            return self._apply_refinement(ind)

        with ThreadPoolExecutor() as executor:
            refined = list(executor.map(refine, pop))
        return refined

    def _adaptive_blocking(self, pop: Population) -> list[tuple[int, int]]:
        ent = np.zeros(self.L)
        for pos in range(self.L):
            # Considere apenas strings que são longas o suficiente
            valid_chars = [ind[pos] for ind in pop if pos < len(ind)]
            if valid_chars:  # Se temos pelo menos um caractere válido
                cnt = Counter(valid_chars)
                probs = np.array(list(cnt.values())) / len(valid_chars)
                ent[pos] = -np.sum(probs * np.log2(probs))
        thr = 0.7 * ent.max() if ent.max() > 0 else 0.0

        blocks = []
        cur = 0
        while cur < self.L:
            length = self.min_block_len if ent[cur] > thr else self.min_block_len * 2
            blocks.append((cur, min(self.L, cur + length)))
            cur += length
        return blocks
