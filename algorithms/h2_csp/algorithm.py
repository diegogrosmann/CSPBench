"""H2-CSP Algorithm (Hybrid Hierarchical Search) unified implementation.

H²-CSP Implementation (Hybrid Hierarchical Search) for CSP.

H²-CSP is a sophisticated hybrid algorithm that solves the Closest String Problem
through a two-layer hierarchical approach, combining structural decomposition,
adaptive technique selection, and global refinement.

ALGORITHMIC ARCHITECTURE:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                          DETAILED H²-CSP ALGORITHM                             │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. B-SPLITTER (Hierarchical Division)                                          │
│   ├── Splits strings into ~√L contiguous blocks                                │
│   ├── Each block: size ≈ √L (balancing heuristic)                              │
│   ├── Preserves spatial locality (contiguous blocks)                           │
│   └── Reduces complexity: O(|Σ|^L) → O(B×|Σ|^(L/B))                          │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. SMART-CORE (Adaptive Block Selection)                                       │
│   ├── For each block, computes difficulty d_b (consensus distance)             │
│   ├── EASY BLOCK (d_b ≤ block_small):                                          │
│   │   ├── Exhaustive search: explores the entire space |Σ|^(r-l)               │
│   │   ├── Guarantees local optimality                                          │
│   │   └── Used when |Σ|^(r-l) ≤ 10,000                                        │
│   ├── MEDIUM BLOCK (block_small < d_b ≤ block_medium):                         │
│   │   ├── Reduced beam search (beam_width/2)                                   │
│   │   ├── Balance between quality and efficiency                               │
│   │   └── Incremental position-by-position construction                        │
│   └── HARD BLOCK (d_b > block_medium):                                         │
│       ├── Full beam search (beam_width)                                        │
│       ├── Maximum computational effort                                         │
│       └── Broad exploration of the search space                                │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 3. GLOBAL REFINE (Fusion and Global Refinement)                                │
│   ├── FUSION: Concatenates best candidates from each block                     │
│   ├── REFINEMENT: Position-by-position hill-climbing                           │
│   │   ├── Tests all possible substitutions                                     │
│   │   ├── Accepts only improvements                                            │
│   │   └── Stops at the
│   └── CONTROLE: Timeout e callback de progresso                               │
└─────────────────────────────────────────────────────────────────────────────────┘

FILOSOFIA ALGORÍTMICA:

• **DECOMPOSIÇÃO INTELIGENTE**: Divide o problema em subproblemas menores e
  independentes, reduzindo drasticamente a complexidade computacional.

• **ADAPTATIVIDADE**: Seleciona a técnica mais apropriada para cada bloco
  baseado na dificuldade, otimizando o uso de recursos computacionais.

• **HIERARQUIA**: Combina otimização local (por bloco) com refinamento global,
  balanceando qualidade local e consistência global.

• **DETERMINISMO**: Algoritmo determinístico que garante reprodutibilidade
  de resultados (exceto por timeout).

CARACTERÍSTICAS DISTINTIVAS:

• **EFICIÊNCIA ESCALÁVEL**: Complexidade O(B×|Σ|^(L/B)) vs. O(|Σ|^L)
• **QUALIDADE ADAPTATIVA**: Técnicas sofisticadas apenas onde necessário
• **ROBUSTEZ**: Fallbacks e controle de timeout para casos extremos for extreme cases
• **FLEXIBILITY**: Configurable parameters for different scenarios

APPLICATION TO CSP:
H²-CSP is especially effective for:
- Medium-sized instances (L=50-500)
- Data with local patterns or hierarchical structure
- Scenarios where solution quality is critical
- Problems requiring determinism and reproducibility

ADVANTAGES:
- Exponential complexity reduction through decomposition
- Automatic adaptation to problem difficulty
- Guarantee of local optimality where computationally feasible
- Global refinement for consistency between blocks

LIMITATIONS:
- Reduced efficiency for very short strings (decomposition overhead)
- Possible suboptimality due to block decomposition
- Additional memory to store candidates from all blocks

FLOW EXAMPLE:
For strings ["ACGTACGT", "AGCTACGT", "ATGTACGT"] with L=8:
1. B-Splitter: √8≈3 blocks → [(0,3), (3,6), (6,8)]
2. Smart-Core per block:
     - Block 1: consensus="ACG", d_b=1 → exhaustive search
     - Block 2: consensus="TAC", d_b=0 → exhaustive search
     - Block 3: consensus="GT", d_b=0 → exhaustive search
3. Fusion: "ACG" + "TAC" + "GT" = "ACGTACGT"
4. Refinement: hill-climbing for final optimization

Detailed algorithmic rationale moved to README to keep source concise.
"""

from __future__ import annotations

from collections import Counter
from itertools import product
from typing import Any, Iterable
import math
import time

from src.domain.algorithms import AlgorithmResult, CSPAlgorithm, register_algorithm

from .config import H2_CSP_DEFAULTS


# ---------------------------------------------------------------------------
# Helper (internal) functions
# ---------------------------------------------------------------------------


def _split_blocks(L: int) -> list[tuple[int, int]]:
    if L <= 0:
        return []
    B = math.ceil(math.sqrt(L))
    base = math.ceil(L / B)
    blocks: list[tuple[int, int]] = []
    cur = 0
    while cur < L:
        blocks.append((cur, min(cur + base, L)))
        cur += base
    return blocks


def _consensus_block(strings: list[str], l: int, r: int) -> str:
    res: list[str] = []
    for pos in range(l, r):
        chars = [s[pos] for s in strings]
        most_common = Counter(chars).most_common(1)[0][0]
        res.append(most_common)
    return "".join(res)


def _exhaustive_block(
    strings: list[str], alphabet: str, l: int, r: int, k: int, limit: int, distance_fn
) -> list[str]:
    m = r - l
    space = len(alphabet) ** m if m > 0 else 1
    segs = [s[l:r] for s in strings]
    candidates: list[tuple[int, str]] = []
    if space <= limit:
        for tpl in product(alphabet, repeat=m):
            cand = "".join(tpl)
            d = max(distance_fn(cand, seg) for seg in segs)
            if len(candidates) < k:
                candidates.append((d, cand))
                candidates.sort(key=lambda x: x[0])
            elif d < candidates[-1][0]:
                candidates.append((d, cand))
                candidates.sort(key=lambda x: x[0])
                del candidates[k:]
    else:  # fallback: use dataset segments + consensus
        seen: set[str] = set()
        for seg in segs:
            if seg in seen:
                continue
            seen.add(seg)
            d = max(distance_fn(seg, other) for other in segs)
            candidates.append((d, seg))
        consensus = _consensus_block(strings, l, r)
        dcons = max(distance_fn(consensus, other) for other in segs)
        candidates.append((dcons, consensus))
    # deduplicate and pick best k
    candidates.sort(key=lambda x: x[0])
    out: list[str] = []
    seen_out: set[str] = set()
    for d, s in candidates:
        if s not in seen_out:
            out.append(s)
            seen_out.add(s)
        if len(out) >= k:
            break
    return out


def _beam_search_block(
    strings: list[str],
    alphabet: str,
    l: int,
    r: int,
    beam_width: int,
    k: int,
    distance_fn,
) -> list[str]:
    m = r - l
    beam = [""]
    # incremental build
    for depth in range(m):
        ext: list[tuple[int, str]] = []
        for pref in beam:
            for c in alphabet:
                new_pref = pref + c
                # partial mismatch count vs segments (direct compare)
                partial_scores = []
                for s in strings:
                    mism = 0
                    seg = s[l : l + len(new_pref)]
                    for i, ch in enumerate(new_pref):
                        if ch != seg[i]:
                            mism += 1
                    partial_scores.append(mism)
                score = max(partial_scores)
                ext.append((score, new_pref))
        ext.sort(key=lambda x: x[0])
        beam = [p for _, p in ext[:beam_width]]
    segs = [s[l:r] for s in strings]
    finals: list[tuple[int, str]] = []
    for cand in beam:
        d = max(distance_fn(cand, seg) for seg in segs)
        finals.append((d, cand))
    finals.sort(key=lambda x: x[0])
    return [c for _, c in finals[:k]]


def _fuse_blocks(
    selected: list[str], blocks: list[tuple[int, int]], total_length: int
) -> str:
    if len(selected) != len(blocks):
        raise ValueError("Selected candidates count mismatch blocks count")
    fused = "".join(selected)
    if len(fused) != total_length:
        raise ValueError("Fused length mismatch original length")
    return fused


def _hill_climb(
    candidate: str,
    strings: list[str],
    distance_center_fn,
    alphabet_by_pos: list[set[str]],
    max_iters: int,
    report_iter,
) -> tuple[str, int, int]:
    current = list(candidate)
    best_dist = distance_center_fn(candidate)
    L = len(current)
    improvements = 0
    for it in range(1, max_iters + 1):
        improved = False
        report_iter(it, max_iters, best_distance=best_dist)
        for pos in range(L):
            original = current[pos]
            for alt in alphabet_by_pos[pos]:
                if alt == original:
                    continue
                current[pos] = alt
                new_str = "".join(current)
                d = distance_center_fn(new_str)
                if d < best_dist:
                    best_dist = d
                    improvements += 1
                    improved = True
                    break
                else:
                    current[pos] = original
            if improved:
                break
        if not improved:
            break
    return "".join(current), best_dist, improvements


# ---------------------------------------------------------------------------
# Algorithm class
# ---------------------------------------------------------------------------


@register_algorithm
class H2CSPAlgorithm(CSPAlgorithm):
    """Hybrid Hierarchical Search (H2CSP).

    Adaptive block‑based CSP algorithm with exhaustive / beam search per block
    followed by global hill‑climbing refinement.
    """

    name = "H2-CSP"  # ASCII registry key
    display_name = "H²-CSP"
    default_params = H2_CSP_DEFAULTS
    supports_internal_parallel = False
    deterministic = True

    # Progress phase boundaries (percent suggestions)
    _PHASES = {
        "init": 0.0,
        "analysis_end": 0.18,
        "generation_end": 0.55,
        "fusion_end": 0.60,
        "initial_eval_end": 0.65,
        "refine_end": 0.95,
        "finish": 1.0,
    }

    def __init__(
        self,
        strings: list[str],
        alphabet: str,
        distance_calculator,
        store=None,
        seed: int | None = None,
        internal_jobs: int = 1,
        **params: Any,
    ):
        # Uniform length validation
        if strings:
            lens = {len(s) for s in strings}
            if len(lens) > 1:
                raise ValueError(
                    f"All strings must share the same length. Found lengths: {sorted(lens)}"
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
        # Normalize basics
        self.params["local_iters"] = max(0, int(self.params.get("local_iters", 3)))
        self.params["beam_width"] = max(1, int(self.params.get("beam_width", 32)))
        self.params["k_candidates"] = max(1, int(self.params.get("k_candidates", 5)))
        # Threshold ordering safeguard
        bs = self.params.get("block_small", 2)
        bm = self.params.get("block_medium", 4)
        bl = self.params.get("block_large", 8)
        if not (bs <= bm <= bl):
            self._report_warning("Block thresholds out of order; auto-sorting applied")
            ordered = sorted([bs, bm, bl])
            (
                self.params["block_small"],
                self.params["block_medium"],
                self.params["block_large"],
            ) = ordered
        # Precompute blocks
        self.blocks = self._prepare_blocks(len(strings[0]) if strings else 0)

    # ------------------------- block preparation -------------------------
    def _prepare_blocks(self, L: int) -> list[tuple[int, int]]:
        if L == 0:
            return []
        if self.params.get("auto_blocks", True):
            blocks = _split_blocks(L)
        else:
            size = max(1, int(self.params.get("block_size", 2)))
            blocks = [(i, min(i + size, L)) for i in range(0, L, size)]
        mb = self.params.get("max_blocks")
        if mb is not None:
            blocks = blocks[: int(mb)]
        # Enforce min block size by merging trailing tiny blocks if necessary
        min_size = self.params.get("min_block_size", 1)
        merged: list[tuple[int, int]] = []
        for b in blocks:
            if not merged:
                merged.append(b)
                continue
            if b[1] - b[0] < min_size:
                last_l, last_r = merged[-1]
                merged[-1] = (last_l, b[1])
            else:
                merged.append(b)
        return merged

    # ------------------------- block analysis ----------------------------
    def _analyze_blocks(self) -> tuple[list[int], list[str], dict[str, float]]:
        difficulties: list[int] = []
        consensuses: list[str] = []
        for l, r in self.blocks:
            cons = _consensus_block(self.strings, l, r)
            segs = [s[l:r] for s in self.strings]
            # distance between candidate block string and block segments via raw mismatch
            d = max(self._block_distance(cons, seg) for seg in segs)
            consensuses.append(cons)
            difficulties.append(d)
        if difficulties:
            stats = {
                "min": min(difficulties),
                "max": max(difficulties),
                "mean": sum(difficulties) / len(difficulties),
            }
        else:
            stats = {"min": 0.0, "max": 0.0, "mean": 0.0}
        return difficulties, consensuses, stats

    def _block_distance(self, a: str, b: str) -> int:
        return sum(1 for x, y in zip(a, b) if x != y)

    # -------------------- per-block candidate generation -----------------
    def _generate_block_candidates(
        self,
        block_index: int,
        l: int,
        r: int,
        difficulty: int,
        seg_consensus: str,
    ) -> tuple[list[str], str]:
        k = self.params["k_candidates"]
        bw = self.params["beam_width"]
        limit = self.params.get("exhaustive_limit", 10_000)
        block_small = self.params.get("block_small")
        block_medium = self.params.get("block_medium")
        block_large = self.params.get("block_large")
        # technique decision
        if difficulty <= block_small:
            tech = "exhaustive"
            cands = _exhaustive_block(
                self.strings,
                self.alphabet,
                l,
                r,
                k,
                limit,
                self._block_distance,
            )
        else:
            if difficulty <= block_medium:
                tech = "beam_reduced"
                bw_local = max(1, bw // 2)
            else:
                tech = "beam_full"
                bw_local = bw
            cands = _beam_search_block(
                self.strings,
                self.alphabet,
                l,
                r,
                bw_local,
                k,
                self._block_distance,
            )
        if not cands:
            cands = [seg_consensus]
            tech += "+fallback"
        return cands, tech

    # ------------------------- global refinement ------------------------
    def _refine(
        self, center: str, start_pct: float, end_pct: float
    ) -> tuple[str, int, int]:
        # alphabet by position
        if not self.strings:
            return center, 0, 0
        L = len(center)
        alph_pos: list[set[str]] = []
        for i in range(L):
            alph_pos.append({s[i] for s in self.strings})
        max_iters = self.params["local_iters"]

        def report_iter(current: int, total: int, **extra):
            if total == 0:
                pct = end_pct
            else:
                pct = start_pct + (end_pct - start_pct) * (current / total)
            if self._monitor:
                pct = (current / total) if total else 1.0
                self._monitor.on_progress(
                    pct,
                    f"refinement iteration {current}/{total}",
                    phase="refinement",
                    **extra,
                )
            # Force progress override (iteration helper already reports pct based on main formula)
            # but we keep separate to enrich metadata.

        refined, best_dist, improvements = _hill_climb(
            center,
            self.strings,
            self.max_distance,
            alph_pos,
            max_iters,
            report_iter,
        )
        return refined, best_dist, improvements

    # ------------------------------ run ---------------------------------
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
                    self._PHASES["init"], "Starting H2CSP", phase="init"
                )

            # Phase 1: analysis
            difficulties, consensuses, stats = self._analyze_blocks()
            if self._monitor:
                self._monitor.on_progress(
                    self._PHASES["analysis_end"],
                    "Block analysis complete",
                    phase="analysis",
                    blocks=len(self.blocks),
                    block_sizes=[r - l for (l, r) in self.blocks],
                    difficulty_min=stats["min"],
                    difficulty_max=stats["max"],
                    difficulty_mean=stats["mean"],
                )

            # Phase 2: per-block generation
            all_candidates: list[list[str]] = []
            techniques: list[str] = []
            for idx, ((l, r), diff, cons) in enumerate(
                zip(self.blocks, difficulties, consensuses), start=1
            ):
                cands, tech = self._generate_block_candidates(idx - 1, l, r, diff, cons)
                all_candidates.append(cands)
                techniques.append(tech)
                # progress within generation segment
                frac = idx / len(self.blocks) if self.blocks else 1.0
                pct = (
                    self._PHASES["analysis_end"]
                    + (self._PHASES["generation_end"] - self._PHASES["analysis_end"])
                    * frac
                )
                if self._monitor:
                    self._monitor.on_progress(
                        pct,
                        f"Block {idx}/{len(self.blocks)} candidates generated",
                        phase="generation",
                        block_index=idx - 1,
                        difficulty=diff,
                        technique=tech,
                        candidates=len(cands),
                    )

            # Phase 3: fusion (choose first candidate per block)
            chosen = [cands[0] for cands in all_candidates] if all_candidates else []
            center = _fuse_blocks(chosen, self.blocks, L)
            initial_dist = self.max_distance(center)
            if self._monitor:
                self._monitor.on_progress(
                    self._PHASES["fusion_end"],
                    "Fusion complete",
                    phase="fusion",
                    initial_distance=initial_dist,
                )

            # Early exit: perfect
            if initial_dist == 0:
                end_time = time.time()
                avg_d = self.average_distance(center)
                total_d = self.total_distance(center)
                if self._monitor:
                    self._monitor.on_progress(
                        self._PHASES["finish"], "Perfect solution", phase="finish"
                    )
                return AlgorithmResult(
                    success=True,
                    center_string=center,
                    max_distance=0,
                    parameters=self.get_actual_params(),
                    error=None,
                    metadata={
                        "algorithm_name": self.name,
                        "display_name": self.display_name,
                        "execution_time": end_time - start_time,
                        "num_strings": len(self.strings),
                        "string_length": L,
                        "alphabet_size": len(self.alphabet),
                        "blocks": len(self.blocks),
                        "techniques": techniques,
                        "refinement_iterations": 0,
                        "refinement_improvements": 0,
                        "avg_distance": avg_d,
                        "total_distance": total_d,
                        "seed": self.seed,
                        "internal_jobs": self.internal_jobs,
                        "deterministic": True,
                    },
                )

            if self._monitor:
                self._monitor.on_progress(
                    self._PHASES["initial_eval_end"],
                    "Initial evaluation done",
                    phase="initial_eval",
                    initial_distance=initial_dist,
                )

            # Phase 4: refinement
            refined, best_dist, improvements = self._refine(
                center, self._PHASES["initial_eval_end"], self._PHASES["refine_end"]
            )
            final_center = refined
            final_dist = best_dist
            end_time = time.time()
            avg_d = self.average_distance(final_center)
            total_d = self.total_distance(final_center)
            if self._monitor:
                self._monitor.on_progress(
                    self._PHASES["finish"], "H2CSP finished", phase="finish"
                )

            return AlgorithmResult(
                success=True,
                center_string=final_center,
                max_distance=final_dist,
                parameters=self.get_actual_params(),
                error=None,
                metadata={
                    "algorithm_name": self.name,
                    "display_name": self.display_name,
                    "execution_time": end_time - start_time,
                    "num_strings": len(self.strings),
                    "string_length": L,
                    "alphabet_size": len(self.alphabet),
                    "blocks": len(self.blocks),
                    "block_sizes": [r - l for (l, r) in self.blocks],
                    "block_difficulties": difficulties,
                    "techniques": techniques,
                    "initial_distance": initial_dist,
                    "refinement_iterations": self.params.get("local_iters", 0),
                    "refinement_improvements": improvements,
                    "avg_distance": avg_d,
                    "total_distance": total_d,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                    "deterministic": True,
                },
            )
        except Exception as e:  # error path
            end_time = time.time()
            msg = f"Error executing H2CSP: {e}"
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
                    "display_name": getattr(self, "display_name", self.name),
                    "execution_time": end_time - start_time,
                    "num_strings": len(self.strings),
                    "string_length": len(self.strings[0]) if self.strings else 0,
                    "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                    "deterministic": True,
                    "error_type": type(e).__name__,
                },
            )
