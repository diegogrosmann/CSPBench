"""Testes adicionais H2CSP (movido para diretório do plugin).

Cobre fallback exhaustive e erro de fuse mismatch.
"""

from __future__ import annotations

import pytest

from src.domain.algorithms import global_registry
from src.domain.distance import create_distance_calculator


def _ensure_h2_registered():
    if "H2CSP" not in global_registry:  # pragma: no cover
        import importlib

        importlib.import_module("algorithms.h2_csp.algorithm")


def build(strings, alphabet, **params):
    _ensure_h2_registered()
    calc = create_distance_calculator("hamming", strings)
    Alg = global_registry["H2-CSP"]
    return Alg(strings, alphabet, distance_calculator=calc, **params)


def test_h2csp_exhaustive_fallback_large_space():
    strings = ["ABCDEF", "ABXDYF", "ABCDQF"]
    alphabet = "ABCDEFGXY"
    alg = build(
        strings,
        alphabet,
        exhaustive_limit=2,
        local_iters=0,
        beam_width=4,
        k_candidates=2,
    )
    res = alg.run()
    assert res["success"]
    assert any("exhaustive" in t for t in res["metadata"]["techniques"])


def test_h2csp_fuse_blocks_mismatch_error():
    strings = ["AAAA", "AAAT"]
    alphabet = "AT"
    _ = build(strings, alphabet, local_iters=0)
    from algorithms.h2_csp.algorithm import _fuse_blocks

    with pytest.raises(ValueError):
        _fuse_blocks(["AA"], [(0, 2), (2, 4)], 4)
