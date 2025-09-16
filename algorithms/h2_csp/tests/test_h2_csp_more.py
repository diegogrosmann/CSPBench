"""Additional tests focused on internal branches of H2CSP to raise coverage >80%.

(File moved from the global suite to the plugin directory.)
"""

from __future__ import annotations

import algorithms.h2_csp.algorithm as h2_mod  # ensures registration
from src.domain.algorithms import global_registry
from src.domain.distance import create_distance_calculator


def ensure_reg():  # small helper
    if "H2CSP" not in global_registry:  # pragma: no cover - defensive redundancy
        import importlib

        importlib.import_module("algorithms.h2_csp")


def build(strings, alphabet, **params):
    ensure_reg()
    calc = create_distance_calculator("hamming", strings)
    Alg = global_registry["H2-CSP"]
    return Alg(strings, alphabet, distance_calculator=calc, **params)


def test_split_blocks_zero_length():
    assert h2_mod._split_blocks(0) == []
    assert h2_mod._split_blocks(-5) == []


def test_exhaustive_block_pruning():
    strings = ["AA", "BB"]
    alphabet = "AB"
    alg = build(
        strings,
        alphabet,
        auto_blocks=False,
        block_size=2,
        block_small=5,
        block_medium=6,
        block_large=7,
        k_candidates=2,
        local_iters=0,
    )
    res = alg.run()
    assert res["success"]
    assert any(t.startswith("exhaustive") for t in res["metadata"]["techniques"])


def test_exhaustive_block_fallback_duplicates():
    strings = ["ABCD", "ABCD", "ABXD"]
    alphabet = "ABCDX"
    alg = build(
        strings,
        alphabet,
        exhaustive_limit=1,
        local_iters=0,
        k_candidates=3,
    )
    res = alg.run()
    assert res["success"]
    assert any("exhaustive" in t for t in res["metadata"]["techniques"])


def test_run_internal_exception(monkeypatch):
    strings = ["AAAA", "AAAT"]
    alphabet = "AT"
    alg = build(strings, alphabet, local_iters=0)

    def boom(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(alg, "_generate_block_candidates", boom)
    res = alg.run()
    assert res["success"] is False
    assert "boom" in (res["error"] or "")
    assert res["metadata"].get("error_type") == "RuntimeError"


def test_beam_search_block_direct():
    strings = ["ABC", "ABD", "ACC"]
    alphabet = "ABCD"
    seg = h2_mod._beam_search_block(
        strings,
        alphabet,
        0,
        3,
        beam_width=2,
        k=2,
        distance_fn=lambda a, b: sum(x != y for x, y in zip(a, b)),
    )
    assert 1 <= len(seg) <= 2


def test_hill_climb_no_improvement():
    strings = ["AAAA", "AAAA"]
    calc = create_distance_calculator("hamming", strings)
    center = "AAAA"
    alph_pos = [{"A"} for _ in range(4)]

    def rpt(*a, **k):  # pragma: no cover
        pass

    new_center, dist, improvements = h2_mod._hill_climb(
        center, strings, calc.max_distance, alph_pos, 5, rpt
    )
    assert new_center == center
    assert dist == 0
    assert improvements == 0


def test_hill_climb_with_improvement():
    strings = ["AAAA", "AAAT"]
    calc = create_distance_calculator("hamming", strings)
    center = "TTTT"
    alph_pos = [{"A", "T"} for _ in range(4)]
    it_calls = []

    def rpt(*a, **k):
        it_calls.append(1)

    new_center, dist, improvements = h2_mod._hill_climb(
        center, strings, calc.max_distance, alph_pos, 10, rpt
    )
    assert dist == calc.max_distance(new_center)
    assert improvements >= 1
    assert len(it_calls) >= 1
