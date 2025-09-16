from src.domain.algorithms import global_registry
from src.domain.distance import create_distance_calculator


def build_algorithm(strings, alphabet, **params):
    calc = create_distance_calculator("hamming", strings)
    AlgCls = global_registry["CSC"]
    return AlgCls(strings, alphabet, distance_calculator=calc, **params)


def test_no_clusters_fallback():
    # Forçar d muito pequeno para não formar clusters (depende dos dados)
    strings = ["AAAA", "TTTT", "CCCC", "GGGG"]
    alg = build_algorithm(strings, "ATCG", d=0, n_blocks=2)
    result = alg.run()
    assert result["success"]
    assert result["metadata"]["fallback_used"] is True


def test_candidate_truncation_flag():
    # Forçar muitos candidatos: usar 3 consensos artificiais com n_blocks alto
    strings = [
        "AAAAAA",
        "AAAACC",
        "CCCCCC",
        "CCCCAA",
        "GGGGGG",
        "GGGGTT",
    ]
    # Parâmetros para tentar gerar mais combinações e limitar rapidamente
    alg = build_algorithm(
        strings,
        "ACGT",
        d=10,
        n_blocks=4,
        max_candidates=2,  # limite bem baixo
    )
    result = alg.run()
    assert result["success"]
    assert result["metadata"]["candidates_truncated"] is True


def test_auto_parameters_flags():
    strings = ["ACGT", "ACGA", "ACGG", "ACGC"]
    alg = build_algorithm(strings, "ACGT")  # sem d/n_blocks explícitos
    result = alg.run()
    assert (
        result["metadata"]["d_auto"] is True
        or result["metadata"]["n_blocks_auto"] is True
    )
