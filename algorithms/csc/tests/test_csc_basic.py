import pytest

from src.domain.algorithms import global_registry
from src.domain.distance import create_distance_calculator


def build_algorithm(strings, alphabet, **params):
    calc = create_distance_calculator("hamming", strings)
    AlgCls = global_registry["CSC"]
    return AlgCls(strings, alphabet, distance_calculator=calc, **params)


def test_returns_algorithm_result_structure():
    strings = ["ACGT", "ACGA", "ACGG"]
    alg = build_algorithm(strings, "ACGT")
    result = alg.run()
    assert result["success"] is True
    assert isinstance(result["center_string"], str)
    assert isinstance(result["max_distance"], int)
    assert isinstance(result["parameters"], dict)
    assert isinstance(result["metadata"], dict)


def test_determinism_same_seed_same_result():
    strings = ["ACGT", "ACGA", "ACGG", "ACGC"]
    alg1 = build_algorithm(strings, "ACGT", seed=42)
    alg2 = build_algorithm(strings, "ACGT", seed=42)
    r1 = alg1.run()
    r2 = alg2.run()
    assert r1["center_string"] == r2["center_string"]
    assert r1["max_distance"] == r2["max_distance"]


def test_string_length_mismatch_error():
    strings = ["AAA", "AAAA"]
    calc = create_distance_calculator(
        "hamming", ["AAA", "AAAA"]
    )  # still for constructor
    from algorithms.csc.algorithm import CSCAlgorithm

    with pytest.raises(ValueError):
        CSCAlgorithm(strings, "ACGT", distance_calculator=calc)
