"""Tests for Baseline (Greedy Consensus) algorithm.

Covers:
- Basic successful execution & distance correctness
- Determinism with same seed & tie_break
- Tie-breaking strategies (lex, first, random)
- Error handling: empty strings list, empty alphabet
- Metadata integrity
- Known example consensus
"""

from __future__ import annotations

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


def build(strings, alphabet, **params) -> BaselineAlg:
    calc = create_distance_calculator("hamming", strings)
    store = params.pop("store", DummyStore())
    return BaselineAlg(
        strings=strings,
        alphabet=alphabet,
        distance_calculator=calc,
        store=store,
        seed=params.pop("seed", 123),
        **params,
    )


def test_baseline_basic_success():
    strings = ["ACGT", "AGGT", "ACGA"]
    alphabet = "ACGT"
    alg = build(strings, alphabet)
    result = alg.run()
    assert result["success"] is True
    center = result["center_string"]
    assert len(center) == 4
    # Validate reported max distance matches calculator
    assert result["max_distance"] == alg.max_distance(center)
    # Basic metadata sanity
    md = result["metadata"]
    assert md["algorithm_name"] == "Baseline"
    assert md["num_strings"] == len(strings)
    assert md["alphabet_size"] == len(alphabet)


def test_baseline_known_example():
    strings = ["ACGT", "AGCT", "ATCT"]
    alphabet = "ACGT"
    alg = build(strings, alphabet, tie_break="lex", seed=1)
    res = alg.run()
    assert res["success"]
    # Example from documentation: expected consensus ACCT distance 1
    assert res["center_string"] == "ACCT"
    assert res["max_distance"] == 1


def test_baseline_determinism_same_seed():
    strings = ["AC", "AG", "AT"]
    alphabet = "ACGT"
    r1 = build(strings, alphabet, seed=99, tie_break="lex").run()
    r2 = build(strings, alphabet, seed=99, tie_break="lex").run()
    assert r1["center_string"] == r2["center_string"]
    assert r1["max_distance"] == r2["max_distance"]


def test_baseline_tie_break_strategies():
    # At position 1 there is a tie deliberately constructed
    strings = ["AA", "BA", "CA"]  # position 0 has A,B,C; position 1 all 'A'
    alphabet = "ABC"
    # For position 0 all letters produce same max prefix distance (=1), tie -> lex => 'A'
    lex = build(strings, alphabet, tie_break="lex", seed=7).run()["center_string"]
    first = build(strings, alphabet, tie_break="first", seed=7).run()["center_string"]
    rnd1 = build(strings, alphabet, tie_break="random", seed=1).run()["center_string"]
    rnd2 = build(strings, alphabet, tie_break="random", seed=1).run()["center_string"]
    assert lex[0] == "A"  # lexicographic min
    assert first[0] == "A"  # first encountered order is alphabet order
    # Random with same seed must be deterministic
    assert rnd1 == rnd2
    assert rnd1[0] in {"A", "B", "C"}


def test_baseline_empty_strings_error():
    alphabet = "ACGT"
    alg = build([], alphabet)
    res = alg.run()
    assert res["success"] is False
    assert "cannot be empty" in (res["error"] or "")


def test_baseline_empty_alphabet_error():
    strings = ["AC", "AG"]
    alg = build(strings, "")
    res = alg.run()
    assert res["success"] is False
    assert "Alphabet" in (res["error"] or "")


def test_baseline_metadata_fields_presence():
    strings = ["AC", "AG", "AT"]
    alphabet = "ACGT"
    res = build(strings, alphabet, seed=5).run()
    assert res["success"]
    md = res["metadata"]
    for key in [
        "execution_time",
        "avg_distance",
        "total_distance",
        "num_strings",
        "string_length",
        "alphabet_size",
        "seed",
        "internal_jobs",
    ]:
        assert key in md, f"Missing metadata key {key}"
