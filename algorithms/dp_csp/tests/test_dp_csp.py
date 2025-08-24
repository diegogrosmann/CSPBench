"""Testes para o algoritmo DP-CSP.

Objetivos:
- Verificar caminho básico de sucesso (linear)
- Verificar estratégia binary retorna mesmo resultado que linear
- Testar caso trivial d=0 (todas iguais)
- Testar abort por complexidade
- Testar erro de comprimentos inconsistentes
"""

from __future__ import annotations

import pytest

from algorithms.dp_csp import DPCSPAlgorithm
from algorithms.baseline import BaselineAlg
from src.domain.distance import create_distance_calculator


class DummyStore:
    def __init__(self):
        self.progress_events: list[tuple[float, str]] = []
        self.warnings: list[str] = []

    def report_algorithm_progress(
        self, progress: float, message: str, **_
    ):  # pragma: no cover
        self.progress_events.append((progress, message))

    def warning(self, message: str):  # pragma: no cover
        self.warnings.append(message)


def build_dp(strings, alphabet, **params) -> DPCSPAlgorithm:
    calc = create_distance_calculator("hamming", strings)
    store = params.pop("store", DummyStore())
    return DPCSPAlgorithm(
        strings=strings,
        alphabet=alphabet,
        distance_calculator=calc,
        store=store,
        seed=params.pop("seed", 123),
        **params,
    )


def test_dp_csp_basic_linear_success():
    strings = ["ACGT", "AGGT", "ACGA"]
    alphabet = "ACGT"
    alg = build_dp(strings, alphabet, search_strategy="linear")
    result = alg.run()
    assert result["success"] is True
    assert result["metadata"]["search_strategy"] == "linear"
    assert result["metadata"]["optimal_d"] == result["max_distance"]


def test_dp_csp_trivial_d_zero():
    strings = ["AAAA", "AAAA", "AAAA"]
    alphabet = "A"
    alg = build_dp(strings, alphabet)
    result = alg.run()
    assert result["success"] is True
    assert result["max_distance"] == 0
    assert result["center_string"] == "AAAA"


def test_dp_csp_binary_matches_linear():
    strings = ["ACGT", "AGGT", "ACGA"]
    alphabet = "ACGT"
    linear = build_dp(strings, alphabet, search_strategy="linear")
    binary = build_dp(strings, alphabet, search_strategy="binary")
    r_lin = linear.run()
    r_bin = binary.run()
    assert r_lin["success"] and r_bin["success"]
    assert r_lin["max_distance"] == r_bin["max_distance"]
    # Centro pode ser diferente mas deve validar distâncias iguais
    assert r_lin["max_distance"] == r_lin["metadata"]["validation_max_distance"]
    assert r_bin["max_distance"] == r_bin["metadata"]["validation_max_distance"]


def test_dp_csp_complexity_abort():
    # Força abort definindo threshold abaixo de (1+1)^n = 2^n
    strings = ["AAAA", "AAAT", "AATA", "TAAA"]
    alphabet = "AT"
    # d=0 => (0+1)^4 = 1 não aborta, d=1 => 2^4=16 > 10 aborta
    alg = build_dp(strings, alphabet, complexity_abort_threshold=10, max_d=5)
    result = alg.run()
    assert result["success"] is False
    assert result["metadata"]["termination_reason"] == "complexity_abort"


def test_dp_csp_inconsistent_lengths_error():
    strings = ["AAA", "AA", "AAAA"]
    alphabet = "A"
    with pytest.raises(ValueError):
        build_dp(strings, alphabet)


def test_dp_csp_matches_baseline_distance_small():
    strings = ["AC", "AG", "AT"]
    alphabet = "ACGT"
    # Baseline center & max distance
    base = BaselineAlg(
        strings=strings,
        alphabet=alphabet,
        distance_calculator=create_distance_calculator("hamming", strings),
        store=DummyStore(),
        seed=42,
    )
    b_res = base.run()

    dp = build_dp(strings, alphabet, search_strategy="binary")
    dp_res = dp.run()
    assert dp_res["success"] is True
    assert (
        dp_res["max_distance"] <= b_res["max_distance"]
    ), "DP exato não deve ser pior que baseline"
