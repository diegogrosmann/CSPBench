"""Tests for the H2CSP algorithm (Hybrid Hierarchical Search).

Coverage objectives (>80% of `algorithm.py`):
- Basic success path with refinement
- Block classification into exhaustive / beam_reduced / beam_full
- Reordering of out-of-order thresholds + warning
- Manual block splitting with merge via `min_block_size`
- Early exit when initial distance = 0 (no refinement)
- Error handling: empty list of strings / empty alphabet
- Main metadata fields and distance consistency
- Parameter diversity (beam_width, k_candidates, local_iters)
"""

from __future__ import annotations

import pytest

from src.domain.algorithms import global_registry
from src.domain.distance import create_distance_calculator


def _ensure_h2_registered():
    if "H2CSP" not in global_registry:
        try:
            import algorithms.h2_csp.algorithm as _h2_mod  # noqa: F401
        except Exception:
            import importlib

            try:
                importlib.import_module("algorithms")
            except Exception:
                pass


class MemoryStore:
    def __init__(self):
        self.progress_events: list[dict] = []
        self.warnings: list[str] = []

    def report_algorithm_progress(
        self, progress: float, message: str, **extra
    ):  # pragma: no cover - simple
        self.progress_events.append({"progress": progress, "message": message, **extra})

    def warning(self, message: str):  # pragma: no cover - simple
        self.warnings.append(message)


def build(strings, alphabet, store=None, **params):
    _ensure_h2_registered()
    calc = create_distance_calculator("hamming", strings)
    Alg = global_registry["H2-CSP"]
    return Alg(strings, alphabet, distance_calculator=calc, store=store, **params)


# ---------------------------------------------------------------------------
# Basic with refinement
# ---------------------------------------------------------------------------


def test_h2csp_basic_success_with_refinement():
    # Data that yields distance > 0 to trigger refinement
    strings = [
        "AAAAB",  # probable center adjusted during refinement
        "AAAAC",
        "BAAAC",
        "AAAAC",
    ]
    alphabet = "ABC"  # diversity
    store = MemoryStore()
    alg = build(
        strings, alphabet, store=store, local_iters=2, beam_width=8, k_candidates=3
    )
    res = alg.run()
    assert res["success"] is True
    center = res["center_string"]
    assert len(center) == len(strings[0])
    # Consistent max_distance
    assert res["max_distance"] == alg.max_distance(center)
    md = res["metadata"]
    assert md["blocks"] >= 1
    assert "techniques" in md
    assert md["refinement_iterations"] == alg.params["local_iters"]
    # At least one final phase progress event
    assert any(ev.get("phase") == "finish" for ev in store.progress_events)


# ---------------------------------------------------------------------------
# Technique classification by difficulty (exhaustive / beam_reduced / beam_full)
# ---------------------------------------------------------------------------


def test_h2csp_block_technique_classification():
    # 3 blocks (L=9) with difficulties 0,1,2 using custom thresholds
    strings = [
        # Block1 (0..2): all AAA => difficulty 0
        # Block2 (3..5): max distance 1 (one string differs in 1 pos)
        # Block3 (6..8): max distance 2 (one string differs in 2 pos)
        "AAABBBCCC",
        "AAABBBCCC",
        "AAABCBCAC",  # changes 1 pos in block2 and 1 pos in block3
        "AAABBBACA",  # changes 2 positions in block3
    ]
    # Ensure consistent length = 9
    assert all(len(s) == 9 for s in strings)
    alphabet = "ABC"
    store = MemoryStore()
    alg = build(
        strings,
        alphabet,
        store=store,
        block_small=0,
        block_medium=1,
        block_large=3,
        local_iters=0,  # remover refinamento para focar em geração
        beam_width=6,
        k_candidates=2,
    )
    res = alg.run()
    assert res["success"]
    techniques = res["metadata"]["techniques"]
    # Ensure the three categories appear across blocks
    assert any(t.startswith("exhaustive") for t in techniques)
    assert any(t.startswith("beam_reduced") for t in techniques)
    assert any(t.startswith("beam_full") for t in techniques)


# ---------------------------------------------------------------------------
# Automatic reordering of out-of-order thresholds
# ---------------------------------------------------------------------------


def test_h2csp_threshold_reordering_warning():
    strings = ["AAAA", "AAAT", "AATT", "TTTT"]
    alphabet = "AT"
    store = MemoryStore()
    alg = build(
        strings,
        alphabet,
        store=store,
        # Purposely invalid order
        block_small=5,
        block_medium=2,
        block_large=3,
        local_iters=0,
        beam_width=4,
    )
    res = alg.run()
    assert res["success"]
    # Warning generated
    assert any("Block thresholds out of order" in w for w in store.warnings)
    # Params reordered
    params = res["parameters"]
    assert params["block_small"] <= params["block_medium"] <= params["block_large"]


# ---------------------------------------------------------------------------
# Manual block splitting + merge of small trailing block
# ---------------------------------------------------------------------------


def test_h2csp_manual_blocks_and_merge():
    # Length 5 with block_size=2 => blocks (0,2),(2,4),(4,5) -> last (< min=2) should be merged
    strings = ["ABCDE", "ABCXE", "ABCYE"]
    alphabet = "ABCDEXY"  # contains all used letters
    store = MemoryStore()
    alg = build(
        strings,
        alphabet,
        store=store,
        auto_blocks=False,
        block_size=2,
        min_block_size=2,
        local_iters=0,
        beam_width=4,
        k_candidates=2,
    )
    # Validate prepared blocks before running
    assert alg.blocks == [(0, 2), (2, 5)], f"Unexpected blocks {alg.blocks}"
    res = alg.run()
    assert res["success"]
    assert len(res["metadata"]["block_sizes"]) == 2


# ---------------------------------------------------------------------------
# Early exit when initial distance = 0
# ---------------------------------------------------------------------------


def test_h2csp_early_exit_perfect_solution():
    strings = ["AAAAAA", "AAAAAA", "AAAAAA"]
    alphabet = "A"
    store = MemoryStore()
    alg = build(strings, alphabet, store=store, local_iters=5)
    res = alg.run()
    assert res["success"]
    assert res["max_distance"] == 0
    # No refinement (improvements 0) in metadata
    assert res["metadata"]["refinement_improvements"] == 0
    # Center identical
    assert res["center_string"] == "AAAAAA"
    # Final progress recorded
    assert any(ev.get("phase") == "finish" for ev in store.progress_events)


# ---------------------------------------------------------------------------
# Input errors (empty strings / empty alphabet)
# ---------------------------------------------------------------------------


def test_h2csp_error_empty_strings():
    alphabet = "AC"
    alg = build([], alphabet)
    res = alg.run()
    assert res["success"] is False
    assert "cannot be empty" in (res["error"] or "")


def test_h2csp_error_empty_alphabet():
    strings = ["AA", "AA"]
    alg = build(strings, "")
    res = alg.run()
    assert res["success"] is False
    assert "Alphabet" in (res["error"] or "")


# ---------------------------------------------------------------------------
# Light parametrization for parameter diversity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "beam_width,k_candidates,local_iters",
    [
        (4, 1, 0),
        (6, 2, 1),
        (8, 3, 2),
    ],
)
def test_h2csp_parameter_diversity(beam_width, k_candidates, local_iters):
    strings = ["ACGT", "ACGA", "ACGG", "ACGC"]
    alphabet = "ACGT"
    alg = build(
        strings,
        alphabet,
        beam_width=beam_width,
        k_candidates=k_candidates,
        local_iters=local_iters,
    )
    res = alg.run()
    assert res["success"]
    assert res["metadata"]["blocks"] >= 1
    assert res["parameters"]["beam_width"] == beam_width
    assert res["parameters"]["k_candidates"] == k_candidates
    assert res["parameters"]["local_iters"] == local_iters


# ---------------------------------------------------------------------------
# Average and total distance consistency
# ---------------------------------------------------------------------------


def test_h2csp_distance_consistency():
    strings = ["AC", "AG", "AT"]
    alphabet = "ACGT"
    alg = build(strings, alphabet, local_iters=1, beam_width=4)
    res = alg.run()
    assert res["success"]
    center = res["center_string"]
    md = res["metadata"]
    # Recompute and compare

    assert md["avg_distance"] == pytest.approx(
        sum(alg.distances_to_all(center)) / len(strings)
    )
    assert md["total_distance"] == sum(alg.distances_to_all(center))
