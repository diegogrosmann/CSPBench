import pytest

from src.domain.dataset import Dataset


def test_dataset_basic_stats_and_properties():
    ds = Dataset(["ACGT", "AGGT", "ACCT"])  # infer alphabet
    assert ds.size == 3
    assert ds.length == 4
    assert set(ds.alphabet) == set("ACGT")

    stats = ds.get_statistics()
    assert stats["size"] == 3
    assert stats["min_length"] == 4
    assert stats["max_length"] == 4
    assert stats["alphabet_size"] == 4
    assert stats["total_characters"] == 12


def test_dataset_add_remove_and_sample():
    ds = Dataset(["AAAA", "AAAT", "AATA", "ATAA", "TAAA"], alphabet="AT")
    ds.add_sequence("TTTT")
    assert ds.size == 6
    with pytest.raises(IndexError):
        ds.remove_sequence(99)
    removed = ds.remove_sequence(0)
    assert removed == "AAAA"
    assert ds.size == 5

    sampled = ds.sample(3, seed=123)
    assert sampled.size == 3
    assert all(len(s) == 4 for s in sampled.get_sequences())


def test_dataset_filter_and_non_uniform_length():
    ds = Dataset(["ACG", "ATG", "AAG", "AC"])
    filtered = ds.filter_by_pattern("C", 1)
    assert filtered.get_sequences() == ["ACG", "AC"]
    # non-uniform -> representative length is max
    assert ds.length == 3


def test_dataset_validation_with_alphabet():
    with pytest.raises(ValueError):
        Dataset(["AX"], alphabet="ACGT")
    # valid when alphabet is None (inferred)
    ds = Dataset(["AX"])  # X becomes part of inferred alphabet
    assert "X" in ds.alphabet
