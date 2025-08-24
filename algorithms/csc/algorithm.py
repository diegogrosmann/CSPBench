"""CSC (Consensus String Clustering) algorithm for the Closest String Problem.

High‑level idea
================
1. Distance analysis → derive automatic parameters (``d`` for DBSCAN radius, ``n_blocks`` for recombination granularity).
2. Density based clustering (DBSCAN over the injected distance metric, using a precomputed matrix) to uncover local structure.
3. Per‑cluster consensus construction (majority vote per position) → local representative strings.
4. Block recombination: each consensus is split into ``n_blocks``; Cartesian product of block choices explores hybrid candidates.
5. Candidate evaluation by maximum distance to the input set (metric‑consistent).
6. Local hill‑climbing refinement (position wise substitutions) on the best candidate.

Why it can work
---------------
Local similarity pockets often contain stronger positional signal than the full heterogeneous instance. Recombining blocks from independent consensuses allows controlled diversification, while the final local search performs intensification.

Key characteristics
-------------------
* Automatic parameter inference (``d`` & ``n_blocks``) from empirical pairwise distance statistics (metric agnostic).
* Deterministic execution (no random sampling) – reproducible given same input.
* Structured progress reporting across 8 macro steps (start → parameters → clustering → consensus → candidates → evaluation → local_search → finish) or fallback path.
* Graceful degradation: if clustering fails or yields zero clusters, falls back to a single global consensus + local search.

Complexity (informal)
---------------------
* Clustering: O(n² · L) worst case for pairwise comparisons within DBSCAN (n strings, length L).
* Consensus build: O(c · s_c · L) (c clusters, s_c average cluster size).
* Recombination: O(c^{n_blocks}) candidates (bounded by ``max_candidates`` safeguard).
* Local search: O(iter · L · |Σ|) with early stopping on no improvement.

Critical parameters (see ``config.py``)
--------------------------------------
* d (int): DBSCAN radius (auto if not provided).
* n_blocks (int): number of blocks used in recombination (auto if not provided).
* max_candidates: hard cap to mitigate combinatorial explosion.
* local_search_max_iterations: limit for refinement phase.
* candidate_eval_batch: granularity of progress events during evaluation.

Returned result (``AlgorithmResult``)
-------------------------------------
``center_string`` (str), ``max_distance`` (int), ``parameters`` (dict of effective parameters), ``metadata`` (timings, flags, counters, diagnostics), and ``success`` flag.

Fallback behaviour
------------------
If no valid clusters are produced: build one global consensus, run local search, mark ``fallback_used`` & possibly ``degraded_mode`` in metadata.

Notes
-----
Metric agnostic: all distance computations (parameter inference, clustering, evaluation, local search) use the injected ``DistanceCalculator``. DBSCAN is fed with a precomputed pairwise distance matrix, so any supported metric is applied consistently end‑to‑end.

This file intentionally consolidates implementation + orchestration for clarity inside the plugin boundary. All domain‑level purity concerns are preserved by only using external libraries inside the algorithm adapter layer.
"""

from __future__ import annotations

import time
import logging
from collections import Counter
from itertools import combinations, product
from typing import Any, List

import numpy as np
from sklearn.cluster import DBSCAN

from src.domain.algorithms import AlgorithmResult, CSPAlgorithm, register_algorithm

from .config import CSC_DEFAULTS

logger = logging.getLogger(__name__)


@register_algorithm
class CSCAlgorithm(CSPAlgorithm):
    name = "CSC"
    default_params = CSC_DEFAULTS
    supports_internal_parallel = False
    deterministic = True

    # Macro phases for progress tracking
    _TOTAL_STEPS = 8

    def __init__(
        self,
        strings: list[str],
        alphabet: str,
        distance_calculator,
        store=None,
        seed: int | None = None,
        internal_jobs: int = 1,
        **params,
    ):
        # Validation: all lengths equal
        if strings:
            lens = {len(s) for s in strings}
            if len(lens) > 1:
                raise ValueError(
                    f"CSC requires strings of equal length. Found lengths: {sorted(lens)}"
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

    # ------------------------- Internal helper methods -------------------------
    def _determine_parameters(self) -> tuple[int, int, bool, bool]:
        d = self.params.get("d")
        n_blocks = self.params.get("n_blocks")
        d_auto_flag = False
        n_blocks_auto_flag = False
        if d is None or n_blocks is None:
            auto_d, auto_n_blocks = self._auto_parameters(self.strings)
            if d is None:
                d = auto_d
                d_auto_flag = True
            if n_blocks is None:
                n_blocks = auto_n_blocks
                n_blocks_auto_flag = True
        return d, n_blocks, d_auto_flag, n_blocks_auto_flag

    def _auto_parameters(self, strings: list[str]) -> tuple[int, int]:
        if len(strings) < 2:
            return CSC_DEFAULTS["min_d"], CSC_DEFAULTS["min_blocks"]
        # Use central distance calculator instead of local duplicate implementation
        pairwise_distances = [self.distance(a, b) for a, b in combinations(strings, 2)]
        mean_distance = float(np.mean(pairwise_distances))
        min_distance = int(np.min(pairwise_distances))
        max_distance = int(np.max(pairwise_distances))
        d = max(
            CSC_DEFAULTS["min_d"],
            int(np.floor(mean_distance * CSC_DEFAULTS["d_factor"])),
        )
        n = len(strings)
        L = len(strings[0])
        n_blocks = max(
            CSC_DEFAULTS["min_blocks"],
            min(
                CSC_DEFAULTS["max_blocks"],
                n // CSC_DEFAULTS["n_div"],
                L // CSC_DEFAULTS["l_div"],
            ),
        )
        logger.info(
            "Automatic parameters: mean=%.2f, min=%d, max=%d | d=%d | n_blocks=%d",
            mean_distance,
            min_distance,
            max_distance,
            d,
            n_blocks,
        )
        return d, n_blocks

    @staticmethod
    def _split_blocks(s: str, n_blocks: int) -> list[str]:
        L = len(s)
        if n_blocks <= 0:
            return [s]
        block_sizes = [L // n_blocks] * n_blocks
        for i in range(L % n_blocks):
            block_sizes[i] += 1
        blocks = []
        idx = 0
        for size in block_sizes:
            blocks.append(s[idx : idx + size])
            idx += size
        return blocks

    @staticmethod
    def _recombine_blocks(blocks_list) -> str:
        return "".join(blocks_list)

    @staticmethod
    def _consensus_string(strings: list[str]) -> str:
        if not strings:
            return ""
        L = len(strings[0])
        result = []
        for i in range(L):
            chars = [s[i] for s in strings]
            most_common = Counter(chars).most_common(1)[0][0]
            result.append(most_common)
        return "".join(result)

    def _compute_distance_matrix(self, strings: list[str]) -> np.ndarray:
        """Build full pairwise distance matrix using injected metric.

        Complexity: O(n^2 * cost(distance))
        """
        n = len(strings)
        mat = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i + 1, n):
                d = self.distance(strings[i], strings[j])
                mat[i, j] = d
                mat[j, i] = d
        return mat

    def _cluster_strings(
        self, strings: list[str], d: int, min_samples: int = 2
    ) -> list[list[str]]:
        if not strings:
            return []
        # Precompute matrix (metric agnostic)
        dist_matrix = self._compute_distance_matrix(strings)
        clustering = DBSCAN(eps=d, min_samples=min_samples, metric="precomputed")
        labels = clustering.fit_predict(dist_matrix)
        clusters: dict[int, list[str]] = {}
        for idx, label in enumerate(labels):
            if label == -1:
                continue  # ignore noise
            clusters.setdefault(label, []).append(strings[idx])
        logger.info("Clusters found: %d", len(clusters))
        return list(clusters.values())

    def _generate_candidates(
        self,
        consensus_list: List[str],
        n_blocks: int,
        L: int,
        max_candidates: int,
    ) -> tuple[list[str], bool]:
        # Split consensus strings into blocks
        blocks_per_consenso = [
            self._split_blocks(cons, n_blocks) for cons in consensus_list
        ]
        from itertools import product

        candidates: list[str] = []
        truncated = False
        for block_tuple in product(*blocks_per_consenso):
            candidate = self._recombine_blocks(block_tuple)
            # Normalização tamanho
            if len(candidate) > L:
                candidate = candidate[:L]
            elif len(candidate) < L:
                candidate += consensus_list[0][len(candidate) : L]
            candidates.append(candidate)
            if len(candidates) >= max_candidates:
                truncated = True
                break
        return candidates, truncated

    def _evaluate_candidates(
        self,
        candidates: list[str],
        batch_size: int,
        start_progress: float,
        end_progress: float,
    ) -> tuple[str, int, int]:
        total = len(candidates)
        best = None
        best_dist = None
        for idx, cand in enumerate(candidates, 1):
            dist = self.max_distance(cand)
            if best is None or dist < best_dist:  # smaller max_distance is better
                best = cand
                best_dist = dist
            # batch progress update
            if idx % batch_size == 0 or idx == total:
                frac = idx / total if total else 1.0
                progress = start_progress + (end_progress - start_progress) * frac
                if self._monitor:
                    self._monitor.on_progress(
                        progress,
                        f"Evaluating candidates ({idx}/{total})",
                        phase="evaluation",
                        step_index=6,
                        total_steps=self._TOTAL_STEPS,
                        candidates_evaluated=idx,
                        candidates_total=total,
                        best_distance_so_far=best_dist,
                    )
        return best, best_dist if best_dist is not None else -1, total

    def _refine_local(
        self,
        candidate: str,
        max_iterations: int,
        start_progress: float,
        end_progress: float,
    ) -> tuple[str, int]:
        # Local search using injected distance calculator
        current = list(candidate)
        iterations = 0
        improved = True
        base_dist = self.max_distance(candidate)
        L = len(current)
        while improved and iterations < max_iterations:
            improved = False
            iterations += 1
            # periodic progress report
            if iterations == 1 or iterations % 5 == 0:
                frac = iterations / max_iterations
                progress = start_progress + (end_progress - start_progress) * min(
                    1.0, frac
                )
                if self._monitor:
                    self._monitor.on_progress(
                        progress,
                        f"Local search iteration {iterations}",
                        phase="local_search",
                        step_index=7,
                        total_steps=self._TOTAL_STEPS,
                        iterations=iterations,
                        base_distance=base_dist,
                    )
            current_dist = self.max_distance("".join(current))
            for i in range(L):
                chars = sorted({s[i] for s in self.strings})
                original = current[i]
                for alt in chars:
                    if alt == original:
                        continue
                    current[i] = alt
                    new_str = "".join(current)
                    new_dist = self.max_distance(new_str)
                    if new_dist < current_dist:
                        current_dist = new_dist
                        improved = True
                        break
                    else:
                        current[i] = original
                if improved:
                    break
        return "".join(current), iterations

    # ------------------------------ Main method ------------------------------
    def run(self) -> AlgorithmResult:  # type: ignore[override]
        start_time = time.time()
        try:
            if not self.strings:
                raise ValueError("String list cannot be empty")
            if not self.alphabet:
                raise ValueError("Alphabet cannot be empty")

            L = len(self.strings[0]) if self.strings else 0
            if self._monitor:
                self._monitor.on_progress(
                    0.0,
                    "Starting CSC",
                    phase="start",
                    step_index=1,
                    total_steps=self._TOTAL_STEPS,
                )

            # Parameter determination
            d, n_blocks, d_auto_flag, n_blocks_auto_flag = self._determine_parameters()
            self.params.setdefault("d", d)
            self.params.setdefault("n_blocks", n_blocks)

            if self._monitor:
                self._monitor.on_progress(
                    0.10,
                    f"Parameters set d={d}, n_blocks={n_blocks}",
                    phase="parameters",
                    step_index=2,
                    total_steps=self._TOTAL_STEPS,
                    d=d,
                    n_blocks=n_blocks,
                    d_auto=d_auto_flag,
                    n_blocks_auto=n_blocks_auto_flag,
                )

            # Clustering
            degraded_mode = False
            try:
                clusters = self._cluster_strings(self.strings, d)
            except Exception as ce:
                self._report_warning(
                    f"Clustering failed: {ce}; using global consensus fallback"
                )
                clusters = []
                degraded_mode = True

            n_clusters = len(clusters)
            if self._monitor:
                self._monitor.on_progress(
                    0.25,
                    f"Clustering completed ({n_clusters} clusters)",
                    phase="clustering",
                    step_index=3,
                    total_steps=self._TOTAL_STEPS,
                    n_clusters=n_clusters,
                    degraded_mode=degraded_mode,
                )

            fallback_used = False

            if n_clusters == 0:
                # Fallback: global consensus + local search
                fallback_used = True
                initial = self._consensus_string(self.strings)
                if self._monitor:
                    self._monitor.on_progress(
                        0.35,
                        "Using global consensus fallback",
                        phase="fallback",
                        step_index=4,
                        total_steps=self._TOTAL_STEPS,
                    )
                refined, ls_iters = self._refine_local(
                    initial,
                    self.params.get(
                        "local_search_max_iterations",
                        CSC_DEFAULTS["local_search_max_iterations"],
                    ),
                    0.90,
                    0.97,
                )
                best_center = refined
                best_dist = self.max_distance(best_center)
                execution_time = time.time() - start_time
                avg_distance = self.average_distance(best_center)
                total_distance = self.total_distance(best_center)
                if self._monitor:
                    self._monitor.on_progress(
                        1.0,
                        "CSC finished (fallback)",
                        phase="finish",
                        step_index=8,
                        total_steps=self._TOTAL_STEPS,
                    )
                return AlgorithmResult(
                    success=True,
                    center_string=best_center,
                    max_distance=best_dist,
                    parameters=self.get_actual_params(),
                    error=None,
                    metadata={
                        "algorithm_name": self.name,
                        "execution_time": execution_time,
                        "num_strings": len(self.strings),
                        "string_length": L,
                        "alphabet_size": len(self.alphabet),
                        "d": d,
                        "n_blocks": n_blocks,
                        "d_auto": d_auto_flag,
                        "n_blocks_auto": n_blocks_auto_flag,
                        "n_clusters": n_clusters,
                        "fallback_used": True,
                        "degraded_mode": degraded_mode,
                        "candidates_generated": 0,
                        "candidates_evaluated": 0,
                        "candidates_truncated": False,
                        "local_search_iterations": ls_iters,
                        "avg_distance": avg_distance,
                        "total_distance": total_distance,
                        "deterministic": True,
                        "seed": self.seed,
                        "internal_jobs": self.internal_jobs,
                    },
                )

            # Local consensus strings
            consensus_list = [self._consensus_string(c) for c in clusters]
            if self._monitor:
                self._monitor.on_progress(
                    0.35,
                    f"{len(consensus_list)} consensus strings computed",
                    phase="consensus",
                    step_index=4,
                    total_steps=self._TOTAL_STEPS,
                )

            # Candidate generation
            max_candidates = self.params.get(
                "max_candidates", CSC_DEFAULTS["max_candidates"]
            )
            candidates, truncated = self._generate_candidates(
                consensus_list, n_blocks, L, max_candidates
            )
            if self._monitor:
                self._monitor.on_progress(
                    0.55,
                    f"Generated {len(candidates)} candidates"
                    + (" (truncated)" if truncated else ""),
                    phase="candidates",
                    step_index=5,
                    total_steps=self._TOTAL_STEPS,
                    candidates_generated=len(candidates),
                    candidates_truncated=truncated,
                )

            if not candidates:
                # additional unlikely fallback
                fallback_used = True
                candidates = [consensus_list[0]]

            # Evaluation
            batch_size = self.params.get(
                "candidate_eval_batch", CSC_DEFAULTS["candidate_eval_batch"]
            )
            best_initial, best_initial_dist, evaluated = self._evaluate_candidates(
                candidates, batch_size, 0.55, 0.75
            )

            # Local search
            refined, ls_iters = self._refine_local(
                best_initial,
                self.params.get(
                    "local_search_max_iterations",
                    CSC_DEFAULTS["local_search_max_iterations"],
                ),
                0.90,
                0.97,
            )

            best_center = refined
            best_dist = self.max_distance(best_center)
            execution_time = time.time() - start_time
            avg_distance = self.average_distance(best_center)
            total_distance = self.total_distance(best_center)

            if self._monitor:
                self._monitor.on_progress(
                    1.0,
                    "CSC finished",
                    phase="finish",
                    step_index=8,
                    total_steps=self._TOTAL_STEPS,
                )

            return AlgorithmResult(
                success=True,
                center_string=best_center,
                max_distance=best_dist,
                parameters=self.get_actual_params(),
                error=None,
                metadata={
                    "algorithm_name": self.name,
                    "execution_time": execution_time,
                    "num_strings": len(self.strings),
                    "string_length": L,
                    "alphabet_size": len(self.alphabet),
                    "d": d,
                    "n_blocks": n_blocks,
                    "d_auto": d_auto_flag,
                    "n_blocks_auto": n_blocks_auto_flag,
                    "n_clusters": n_clusters,
                    "fallback_used": fallback_used,
                    "degraded_mode": degraded_mode,
                    "candidates_generated": len(candidates),
                    "candidates_evaluated": evaluated,
                    "candidates_truncated": truncated,
                    "best_candidate_initial_distance": best_initial_dist,
                    "local_search_iterations": ls_iters,
                    "avg_distance": avg_distance,
                    "total_distance": total_distance,
                    "deterministic": True,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                },
            )
        except Exception as e:  # general error
            execution_time = time.time() - start_time
            msg = f"Error during CSC execution: {e}"
            if self._monitor:
                self._monitor.on_warning(msg)
                self._monitor.on_progress(
                    0.0,
                    msg,
                    phase="error",
                    step_index=self._TOTAL_STEPS,
                    total_steps=self._TOTAL_STEPS,
                    error_type=type(e).__name__,
                )
            return AlgorithmResult(
                success=False,
                center_string="",
                max_distance=-1,
                parameters=self.get_actual_params(),
                error=msg,
                metadata={
                    "algorithm_name": self.name,
                    "execution_time": execution_time,
                    "num_strings": len(self.strings),
                    "string_length": len(self.strings[0]) if self.strings else 0,
                    "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                    "deterministic": True,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                    "error_type": type(e).__name__,
                },
            )
