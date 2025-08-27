"""Cobertura adicional de caminhos H2CSP.

Inclui branches específicos não exercitados pelos testes principais:
- Substituição do pior candidato em _exhaustive_block (k=1)
- Comportamento de *erro* quando ``max_blocks`` trunca cobertura (gera ValueError via _fuse_blocks)
- Caso strings de comprimento zero (blocos vazios -> stats zeradas)
- Erro de _fuse_blocks (mismatch quantidade de blocos)
"""

from __future__ import annotations

import pytest

from src.domain.distance import create_distance_calculator
from src.domain.algorithms import global_registry

import algorithms.h2_csp.algorithm as h2_mod  # garante registro


def _ensure_reg():  # pragma: no cover - simples salvaguarda
    if "H2CSP" not in global_registry:
        import importlib

        importlib.import_module("algorithms.h2_csp.algorithm")


def _build(strings, alphabet, **params):
    _ensure_reg()
    calc = create_distance_calculator("hamming", strings)
    Alg = global_registry["H2-CSP"]
    return Alg(strings, alphabet, distance_calculator=calc, **params)


# ---------------------------------------------------------------------------
# Substituição do pior candidato (branch elif d < candidates[-1][0])
# ---------------------------------------------------------------------------


def test_exhaustive_block_replacement_branch():
    strings = ["AA", "BB"]
    alphabet = "AB"
    alg = _build(
        strings,
        alphabet,
        auto_blocks=False,
        block_size=2,
        block_small=5,
        block_medium=6,
        block_large=7,
        k_candidates=1,  # força substituições
        local_iters=0,
    )
    res = alg.run()
    assert res["success"]
    assert any(t.startswith("exhaustive") for t in res["metadata"]["techniques"])


# ---------------------------------------------------------------------------
# max_blocks reduzindo blocos
# ---------------------------------------------------------------------------


def test_max_blocks_truncation_error():
    """Truncar blocos antes de cobrir todo o comprimento resulta em erro.

    Justificativa: ``_prepare_blocks`` gera blocos que cobrem todo L; aplicar
    ``max_blocks`` com um valor menor que o número necessário causa perda de
    cobertura e, na fusão, ``_fuse_blocks`` detecta comprimento agregado menor
    que L, levantando ``ValueError``. O algoritmo então retorna ``success=False``.
    Este teste valida esse comportamento (em vez de assumir sucesso inválido).
    """
    strings = ["AAAAAAAAAAAA", "AAAABAAAACAA"]  # L=12 -> 4 blocos auto (3,3,3,3)
    alphabet = "ABC"
    alg = _build(strings, alphabet, max_blocks=2, local_iters=0, beam_width=4)
    assert len(alg.blocks) == 2  # truncado para dois blocos parciais (0..3,3..6)
    res = alg.run()
    assert res["success"] is False
    assert "Fused length mismatch" in (res["error"] or "")


# ---------------------------------------------------------------------------
# Strings de comprimento zero (lista não vazia) => nenhum bloco
# ---------------------------------------------------------------------------


def test_zero_length_strings_blocks_empty():
    strings = ["", "", ""]
    alphabet = "A"
    alg = _build(strings, alphabet, local_iters=0)
    res = alg.run()
    assert res["success"]
    md = res["metadata"]
    assert md["blocks"] == 0
    # Distância inicial/final deve ser 0 (centro string vazia)
    assert res["max_distance"] == 0
    assert res["center_string"] == ""


# ---------------------------------------------------------------------------
# _fuse_blocks erro por quantidade (primeiro if)
# ---------------------------------------------------------------------------


def test_fuse_blocks_selected_count_mismatch():
    with pytest.raises(ValueError):
        h2_mod._fuse_blocks(["AA"], [(0, 2), (2, 4)], 4)
