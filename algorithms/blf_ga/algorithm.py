"""BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) implementation.

Implementation of the BLF-GA heuristic (Blockwise Learning Fusion + Genetic Algorithm) for CSP.

BLF-GA is an advanced hybrid metaheuristic that combines three main strategies:
1. **Blockwise Learning**: Divides strings into blocks and learns local patterns
2. **Genetic Algorithm**: Uses genetic operations for global evolution
3. **Fusion**: Combines local and global knowledge for better convergence

ALGORITHMIC ARCHITECTURE:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                           DETAILED BLF-GA ALGORITHM                             │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. INITIALIZATION                                                              │
│   ├── Create hybrid initial population:                                        │
│   │   ├── Consensus of original strings (initial quality)                     │
│   │   ├── Block-wise consensus variations (intelligent diversity)             │
│   │   └── Random strings (broad exploration)                                  │
│   └── Divide strings into adaptive blocks based on initial_blocks             │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. MAIN LOOP (until stopping criterion)                                        │
│   ├── A) BLOCKWISE LEARNING                                                   │
│   │   ├── For each block, find consensus pattern of current population        │
│   │   ├── Create local knowledge repository                                   │
│   │   └── Use patterns to guide genetic operations                            │
│   ├── B) GENETIC EVOLUTION                                                    │
│   │   ├── Tournament selection (balanced selective pressure)                  │
│   │   ├── Adaptive crossover:                                                 │
│   │   │   ├── one-point: preserves large segments                             │
│   │   │   ├── uniform: maximum recombination                                  │
│   │   │   └── blend-blocks: preserves local patterns                          │
│   │   ├── Adaptive mutation:                                                  │
│   │   │   ├── multi: point alterations                                        │
│   │   │   ├── inversion: structural reorganization                            │
│   │   │   └── transposition: segment relocation                               │
│   │   └── Adaptive elitism (can be disabled for diversity)                   │
│   ├── C) ADAPTIVE MECHANISMS                                                  │
│   │   ├── Block redivision based on positional entropy                       │
│   │   ├── Random immigrants to maintain diversity                             │
│   │   ├── Adaptive mutation based on convergence and diversity               │
│   │   ├── Intensive local refinement of best individuals                     │
│   │   ├── Niching to preserve multiple solutions                             │
│   │   └── Population restart in case of stagnation                           │
│   └── D) STOPPING CRITERIA                                                    │
│       ├── Optimal solution found (distance = 0)                               │
│       ├── Time/generation limit reached                                       │
│       └── Early stopping due to stagnation (no_improve_patience)             │
└─────────────────────────────────────────────────────────────────────────────────┘

DISTINCTIVE CHARACTERISTICS:

• **INTELLIGENT HYBRIDIZATION**: Combines global search (GA) with local learning (blocks)
• **DYNAMIC ADAPTIVITY**: Parameters automatically adjust during evolution
• **CONTROLLED DIVERSITY**: Multiple mechanisms to avoid premature convergence
• **EFFICIENT PARALLELIZATION**: Leverages multiple cores for evaluation and refinement
• **ADVANCED CONFIGURABILITY**: Dynamic parameters and interchangeable strategies

CSP APPLICATION:
BLF-GA is especially suitable for the Closest String Problem because:
- Blocks capture common local patterns between strings
- Learning adapts to the specific problem structure
- Hybridization balances global exploration with local intensification
- Adaptive mechanisms avoid premature convergence traps

EXAMPLE FLOW:
Strings: ["ACGT", "AGCT", "ATCT"]
1. Initial consensus: "ACCT"
2. Blocks: [(0,2), (2,4)] based on entropy
3. Learning: block1="AC", block2="CT"
4. Evolution: crossover_blend_blocks uses learned patterns
5. Adaptation: redivides blocks if entropy changes
6. Convergence: refines solutions until optimal or stopping criterion


Refactored to conform to new CSPAlgorithm contract (AlgorithmResult, progress
reporting via store callbacks, no separate implementation module).
"""

from __future__ import annotations

import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from src.domain.algorithms import AlgorithmResult, CSPAlgorithm, register_algorithm

from .config import BLF_GA_DEFAULTS
from .ops import genetic_ops


@register_algorithm
class BLFGAAlgorithm(CSPAlgorithm):
    name = "BLF-GA"
    default_params: dict = BLF_GA_DEFAULTS
    supports_internal_parallel = True

    def __init__(
        self,
        strings: list[str],
        alphabet: str,
        distance_calculator,
        store=None,
        seed: int | None = None,
        internal_jobs: int = 1,
        **params: Any,
    ) -> None:
        if strings:
            lens = {len(s) for s in strings}
            if len(lens) > 1:
                raise ValueError(
                    f"BLF-GA requires strings of equal length. Found lengths: {sorted(lens)}"
                )
        super().__init__(
            strings,
            alphabet,
            distance_calculator=distance_calculator,
            store=store,
            seed=seed,
            internal_jobs=internal_jobs,
            **params,
        )
        self._normalize_parameters()
        self.blocks: list[tuple[int, int]] = self._initial_blocking()
        # Counters
        self.generations_executed = 0
        self.improvement_generations = 0
        self.immigrants_injections = 0
        self.mutation_adaptations = 0
        self.mutation_resets = 0
        self.restarts_executed = 0
        self.block_redivisions = 0
        self.elites_refined_total = 0
        self._no_improve = 0
        self._mut_prob_base = self.params["mut_prob"]
        self._mut_boost_timer = 0
        self.initial_fitness: int | None = None

    # -------- parameter normalization --------
    def _resolve_dynamic(
        self, value, ref: int, mode: str, min_value: int | None = None
    ) -> int:
        if mode == "multiplier":
            if isinstance(value, float) and value >= 1:
                resolved = int(value * ref)
            else:
                resolved = int(value)
        else:
            if isinstance(value, float) and 0 < value <= 1:
                resolved = int(value * ref)
            else:
                resolved = int(value)
        if min_value is not None:
            resolved = max(min_value, resolved)
        return max(1, resolved)

    def _normalize_parameters(self) -> None:
        p = self.params
        n = len(self.strings) if self.strings else 0
        L = len(self.strings[0]) if self.strings else 0
        p["pop_size"] = self._resolve_dynamic(
            p.get("pop_size", BLF_GA_DEFAULTS["pop_size"]), n if n > 0 else 1, "multiplier", p.get("min_pop_size", BLF_GA_DEFAULTS["min_pop_size"])  # type: ignore[arg-type]
        )
        p["initial_blocks"] = self._resolve_dynamic(
            p.get("initial_blocks", BLF_GA_DEFAULTS["initial_blocks"]),
            L if L > 0 else 1,
            "proportion",
        )
        nip = p.get("no_improve_patience", 0)
        max_g = int(p.get("max_gens", BLF_GA_DEFAULTS["max_gens"]))
        if isinstance(nip, float) and 0 < nip < 1:
            p["no_improve_patience"] = max(1, int(nip * max_g))
        else:
            p["no_improve_patience"] = int(nip)
        p["max_gens"] = max_g
        p["mut_prob"] = float(p.get("mut_prob", BLF_GA_DEFAULTS["mut_prob"]))
        p["cross_prob"] = float(p.get("cross_prob", BLF_GA_DEFAULTS["cross_prob"]))
        p["elite_rate"] = float(p.get("elite_rate", BLF_GA_DEFAULTS["elite_rate"]))
        if p["pop_size"] < 2:
            p["pop_size"] = 2
        if L and p["initial_blocks"] > L:
            p["initial_blocks"] = L
        self._elite_count = max(1, int(p["elite_rate"] * p["pop_size"]))

    # -------- blocks & population --------
    def _initial_blocking(self) -> list[tuple[int, int]]:
        L = len(self.strings[0]) if self.strings else 0
        if L == 0:
            return []
        size = max(
            self.params.get("min_block_len", 1), L // self.params["initial_blocks"]
        )
        return [(i, min(i + size, L)) for i in range(0, L, size)]

    def _adaptive_blocking(self, pop: list[str]) -> list[tuple[int, int]]:
        import math

        import numpy as np

        L = len(self.strings[0]) if self.strings else 0
        if L == 0:
            return []
        ent = np.zeros(L)
        for pos in range(L):
            chars = [ind[pos] for ind in pop]
            cnt = Counter(chars)
            probs = [c / len(chars) for c in cnt.values()]
            ent[pos] = -sum(p * math.log2(p) for p in probs if p > 0)
        thr = 0.7 * ent.max() if ent.max() > 0 else 0.0
        min_len = self.params.get("min_block_len", 1)
        blocks: list[tuple[int, int]] = []
        cur = 0
        while cur < L:
            length = min_len if ent[cur] > thr else min_len * 2
            blocks.append((cur, min(L, cur + length)))
            cur += length
        return blocks

    def _build_initial_population(self) -> list[str]:
        if not self.strings:
            return []
        L = len(self.strings[0])
        pop_size = self.params["pop_size"]
        consensus = "".join(
            Counter(chars).most_common(1)[0][0] for chars in zip(*self.strings)
        )
        pop: list[str] = [consensus]
        for _ in range(pop_size // 3):
            seq = list(consensus)
            for l, r in self.blocks:
                if self.rng.random() < 0.5:
                    src = self.rng.choice(self.strings)
                    seq[l:r] = src[l:r]
            pop.append("".join(seq))
        while len(pop) < pop_size:
            pop.append("".join(self.rng.choice(self.alphabet) for _ in range(L)))
        return pop

    # -------- evaluation --------
    def _evaluate_population(self, pop: list[str]):
        if self.internal_jobs <= 1:
            return [(s, self.max_distance(s)) for s in pop]
        with ThreadPoolExecutor(max_workers=self.internal_jobs) as ex:
            return list(ex.map(lambda s: (s, self.max_distance(s)), pop))

    def _sort_population(self, pop: list[str]) -> list[str]:
        evaluated = self._evaluate_population(pop)
        evaluated.sort(key=lambda x: x[1])
        return [s for s, _ in evaluated]

    def _population_diversity(self, pop: list[str]) -> float:
        return float(genetic_ops.mean_hamming_distance(pop)) if pop else 0.0

    # -------- operators --------
    def _select_parent(self, pop_sorted: list[str]) -> str:
        k = self.params.get("tournament_k", 2)
        sample = self.rng.sample(pop_sorted, min(k, len(pop_sorted)))
        return min(sample, key=lambda s: self.max_distance(s))

    def _crossover(self, p1: str, p2: str) -> tuple[str, str]:
        ctype = self.params.get("crossover_type")
        if ctype == "one_point":
            return genetic_ops.crossover_one_point(p1, p2, self.rng)
        if ctype == "uniform":
            return genetic_ops.crossover_uniform(p1, p2, self.rng)
        if ctype == "blend_blocks":
            return genetic_ops.crossover_blend_blocks(p1, p2, self.blocks, self.rng)
        return p1, p2

    def _mutate(self, ind: str) -> str:
        mtype = self.params.get("mutation_type")
        if mtype == "multi":
            return genetic_ops.mutate_multi(
                ind, self.alphabet, self.rng, n=self.params.get("mutation_multi_n", 2)
            )
        if mtype == "inversion":
            return genetic_ops.mutate_inversion(ind, self.rng)
        if mtype == "transposition":
            return genetic_ops.mutate_transposition(ind, self.rng)
        return ind

    def _refine_individuals(self, elites: list[str]) -> list[str]:
        rtype = self.params.get("refinement_type")
        if rtype == "greedy":
            func = genetic_ops.refine_greedy
        elif rtype == "swap":
            func = genetic_ops.refine_swap
        elif rtype == "insertion":
            func = genetic_ops.refine_insertion
        elif rtype == "2opt":
            func = genetic_ops.refine_2opt
        else:
            return elites
        if self.internal_jobs > 1 and len(elites) > 1:
            with ThreadPoolExecutor(
                max_workers=min(self.internal_jobs, len(elites))
            ) as ex:
                return list(ex.map(lambda e: func(e, self.strings), elites))
        return [func(e, self.strings) for e in elites]

    def _apply_niching(self, pop: list[str]) -> list[str]:
        if not self.params.get("niching", False) or len(pop) < 2:
            return pop
        radius = self.params.get("niching_radius", 3)
        fit = [(s, self.max_distance(s)) for s in pop]
        fit.sort(key=lambda x: x[1])
        selected: list[str] = []
        for s, _ in fit:
            if all(
                genetic_ops.mean_hamming_distance([s, t]) >= radius for t in selected
            ):
                selected.append(s)
            if len(selected) >= self.params["pop_size"]:
                break
        while len(selected) < self.params["pop_size"]:
            L = len(self.strings[0])
            selected.append("".join(self.rng.choice(self.alphabet) for _ in range(L)))
        return selected

    # -------- learning placeholder --------
    def _learn_blocks(self, pop: list[str]) -> None:
        # Placeholder to maintain structure (consensus per block if needed)
        return None

    # -------- adaptive mechanisms --------
    def _immigration(self, pop: list[str]) -> int:
        freq = self.params.get("immigrant_freq", 0)
        if not freq or self.generations_executed % freq != 0:
            return 0
        ratio = self.params.get("immigrant_ratio", 0.0)
        count = max(0, int(ratio * self.params["pop_size"]))
        if count == 0:
            return 0
        L = len(self.strings[0])
        for i in range(1, count + 1):
            pop[-i] = "".join(self.rng.choice(self.alphabet) for _ in range(L))
        self.immigrants_injections += 1
        return count

    def _adaptive_mutation(self) -> str:
        diversity = self._population_diversity(self._current_population)
        threshold = self.params.get("diversity_threshold", 0.0) * (
            len(self.strings[0]) if self.strings else 0
        )
        triggered = diversity < threshold if threshold > 0 else False
        N = self.params.get("mutation_adapt_N", 0)
        if N and self._no_improve >= N:
            triggered = True
        if triggered and self._mut_boost_timer == 0:
            self.params["mut_prob"] = self._mut_prob_base * self.params.get(
                "mutation_adapt_factor", 2.0
            )
            self._mut_boost_timer = self.params.get("mutation_adapt_duration", 5)
            self.mutation_adaptations += 1
        if self._mut_boost_timer > 0:
            self._mut_boost_timer -= 1
            if self._mut_boost_timer == 0:
                self.params["mut_prob"] = self._mut_prob_base
                self.mutation_resets += 1
        return "boosted" if self.params["mut_prob"] != self._mut_prob_base else "normal"

    def _restart(self, pop: list[str]) -> bool:
        patience = self.params.get("restart_patience", 0)
        if patience and self._no_improve >= patience:
            ratio = self.params.get("restart_ratio", 0.0)
            count = int(ratio * self.params["pop_size"])
            L = len(self.strings[0])
            for i in range(count):
                pop[-(i + 1)] = "".join(
                    self.rng.choice(self.alphabet) for _ in range(L)
                )
            self._no_improve = 0
            self.restarts_executed += 1
            return True
        return False

    def _maybe_redivide_blocks(self, pop: list[str]) -> bool:
        rf = self.params.get("rediv_freq", 0)
        if rf and self.generations_executed % rf == 0:
            self.blocks = self._adaptive_blocking(pop)
            self.block_redivisions += 1
            return True
        return False

    # -------- main run --------
    def run(self) -> AlgorithmResult:  # type: ignore[override]
        start_time = time.time()
        try:
            if not self.strings:
                raise ValueError("String list cannot be empty")
            if not self.alphabet:
                raise ValueError("Alphabet cannot be empty")
            L = len(self.strings[0])
            if self._monitor:
                self._monitor.on_progress(
                    0.0,
                    "Starting BLF-GA",
                    phase="init",
                    pop_size=self.params["pop_size"],
                    blocks=len(self.blocks),
                )

            # Verificação inicial de cancelamento
            if self._monitor and self._monitor.is_cancelled():
                return self._build_cancelled_result(start_time)

            population = self._build_initial_population()
            population = self._sort_population(population)
            self._current_population = population  # for adaptive mutation diversity
            best = population[0]
            best_val = self.max_distance(best)
            self.initial_fitness = best_val
            if self._monitor:
                self._monitor.on_progress(
                    0.05,
                    "Initial population built",
                    phase="population_init",
                    initial_fitness=best_val,
                )
            max_gens = self.params["max_gens"]
            no_improve_patience = self.params.get("no_improve_patience", 0)
            termination_reason = "max_generations"
            for gen in range(1, max_gens + 1):
                # Verificar cancelamento a cada geração
                if self._monitor and self._monitor.is_cancelled():
                    termination_reason = "cancelled"
                    best = population[0] if population else ""
                    best_val = self.max_distance(best) if best else -1
                    break

                self.generations_executed = gen
                elapsed = time.time() - start_time
                if elapsed >= self.params.get("max_time", float("inf")):
                    termination_reason = "time_limit"
                    break
                immigrants = self._immigration(population)
                mutation_mode = self._adaptive_mutation()
                self._learn_blocks(population)
                sorted_pop = self._sort_population(population)
                new_pop: list[str] = []
                disable_elitism_gens = self.params.get("disable_elitism_gens", 0)
                elitism_enabled = not (
                    disable_elitism_gens and gen % disable_elitism_gens == 0
                )
                elites_refined = 0
                if elitism_enabled:
                    elites = sorted_pop[: self._elite_count]
                    new_pop.extend(elites)
                while len(new_pop) < self.params["pop_size"]:
                    p1 = self._select_parent(sorted_pop)
                    p2 = self._select_parent(sorted_pop)
                    if self.rng.random() < self.params["cross_prob"]:
                        c1, c2 = self._crossover(p1, p2)
                    else:
                        c1, c2 = p1, p2
                    if self.rng.random() < self.params["mut_prob"]:
                        c1 = self._mutate(c1)
                    if (
                        len(new_pop) + 1 < self.params["pop_size"]
                        and self.rng.random() < self.params["mut_prob"]
                    ):
                        c2 = self._mutate(c2)
                    new_pop.append(c1)
                    if len(new_pop) < self.params["pop_size"]:
                        new_pop.append(c2)
                new_pop = self._apply_niching(new_pop)
                population = self._sort_population(new_pop)
                self._current_population = population
                refine_policy = self.params.get("refine_elites", "best")
                if refine_policy in {"best", "all"} and population:
                    if refine_policy == "best":
                        refined = self._refine_individuals([population[0]])
                        population[0] = refined[0]
                        elites_refined = 1
                    else:
                        refined = self._refine_individuals(
                            population[: self._elite_count]
                        )
                        population[: self._elite_count] = refined
                        elites_refined = self._elite_count
                    self.elites_refined_total += elites_refined
                current_best = population[0]
                current_val = self.max_distance(current_best)
                improved = current_val < best_val
                if improved:
                    best, best_val = current_best, current_val
                    self._no_improve = 0
                    self.improvement_generations += 1
                else:
                    self._no_improve += 1
                restarted = self._restart(population)
                block_rediv = self._maybe_redivide_blocks(population)
                progress = 0.05 + (gen / max_gens) * 0.80
                diversity = self._population_diversity(population)
                if self._monitor:
                    self._monitor.on_progress(
                        progress,
                        f"Gen {gen} best={best_val}",
                        phase="evolution",
                        generation=gen,
                        generations_total=max_gens,
                        best_fitness=best_val,
                        current_fitness=current_val,
                        improvement=improved,
                        stagnation_count=self._no_improve,
                        immigrants_injected=immigrants,
                        mutation_prob=self.params["mut_prob"],
                        mutation_mode=mutation_mode,
                        restarted=restarted,
                        block_redivided=block_rediv,
                        elites_refined=elites_refined,
                        elitism_enabled=elitism_enabled,
                        diversity=diversity,
                        elapsed_seconds=time.time() - start_time,
                    )
                if best_val == 0:
                    termination_reason = "optimal_found"
                    break
                if no_improve_patience and self._no_improve >= no_improve_patience:
                    termination_reason = "early_stopping"
                    break

            # Verificar cancelamento antes do refinamento final
            if (
                termination_reason != "cancelled"
                and best_val > 0
                and self.params.get("refinement_type")
                and self.params.get("refine_elites")
            ):
                if self._monitor and self._monitor.is_cancelled():
                    termination_reason = "cancelled"
                else:
                    if self._monitor:
                        self._monitor.on_progress(
                            0.85,
                            "Final refinement pass",
                            phase="refinement_start",
                            best_fitness_pre=best_val,
                        )
                    refined = self._refine_individuals([best])
                    refined_best = refined[0]
                    refined_val = self.max_distance(refined_best)
                    if refined_val < best_val:
                        best, best_val = refined_best, refined_val
                        self.improvement_generations += 1
                    if self._monitor:
                        self._monitor.on_progress(
                            0.95,
                            "Refinement complete",
                            phase="refinement_end",
                            best_fitness_post=best_val,
                        )
            else:
                if self._monitor:
                    self._monitor.on_progress(
                        0.90,
                        "Skipping final refinement",
                        phase="refinement_skip",
                        best_fitness=best_val,
                    )

            end_time = time.time()

            # Para execução cancelada, assegurar que temos pelo menos algum resultado
            if termination_reason == "cancelled":
                if not best:
                    return self._build_cancelled_result(start_time)

            avg_dist = self.average_distance(best)
            total_dist = self.total_distance(best)
            if self._monitor:
                self._monitor.on_progress(
                    1.0,
                    "BLF-GA finished",
                    phase="finish",
                    termination_reason=termination_reason,
                    best_fitness=best_val,
                )
            metadata = {
                "algorithm_name": self.name,
                "execution_time": end_time - start_time,
                "num_strings": len(self.strings),
                "string_length": L,
                "alphabet_size": len(self.alphabet),
                "population_size": self.params["pop_size"],
                "generations_executed": self.generations_executed,
                "best_fitness": best_val,
                "initial_fitness": self.initial_fitness,
                "improvement_generations": self.improvement_generations,
                "immigrants_injections": self.immigrants_injections,
                "mutation_adaptations": self.mutation_adaptations,
                "mutation_resets": self.mutation_resets,
                "restarts_executed": self.restarts_executed,
                "block_redivisions": self.block_redivisions,
                "elites_refined_total": self.elites_refined_total,
                "refinement_type": self.params.get("refinement_type"),
                "refine_elites": self.params.get("refine_elites"),
                "termination_reason": termination_reason,
                "seed": self.seed,
                "internal_jobs": self.internal_jobs,
                "deterministic": False,
                "final_mutation_prob": self.params.get("mut_prob"),
                "diversity_last": self._population_diversity(self._current_population),
                "avg_distance": avg_dist,
                "total_distance": total_dist,
            }

            # Se foi cancelado, retornar como falha
            if termination_reason == "cancelled":
                return AlgorithmResult(
                    success=False,
                    center_string=best,
                    max_distance=best_val,
                    parameters=self.get_actual_params(),
                    error="Algorithm execution was cancelled",
                    metadata=metadata,
                )

            return AlgorithmResult(
                success=True,
                center_string=best,
                max_distance=best_val,
                parameters=self.get_actual_params(),
                error=None,
                metadata=metadata,
            )
        except Exception as e:  # noqa: BLE001
            end_time = time.time()
            msg = f"Error executing BLF-GA: {e}"
            if self._monitor:
                self._monitor.on_warning(msg)
                self._monitor.on_progress(
                    0.0, msg, phase="error", error_type=type(e).__name__
                )
            return AlgorithmResult(
                success=False,
                center_string="",
                max_distance=-1,
                parameters=self.get_actual_params(),
                error=msg,
                metadata={
                    "algorithm_name": self.name,
                    "execution_time": end_time - start_time,
                    "num_strings": len(self.strings),
                    "string_length": len(self.strings[0]) if self.strings else 0,
                    "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                    "deterministic": False,
                    "error_type": type(e).__name__,
                },
            )

    def _build_cancelled_result(self, start_time: float) -> AlgorithmResult:
        """Constrói resultado para execução cancelada."""
        end_time = time.time()
        return AlgorithmResult(
            success=False,
            center_string="",
            max_distance=-1,
            parameters=self.get_actual_params(),
            error="Algorithm execution was cancelled",
            metadata={
                "algorithm_name": self.name,
                "execution_time": end_time - start_time,
                "num_strings": len(self.strings),
                "string_length": len(self.strings[0]) if self.strings else 0,
                "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                "seed": self.seed,
                "internal_jobs": self.internal_jobs,
                "deterministic": False,
                "termination_reason": "cancelled",
            },
        )
