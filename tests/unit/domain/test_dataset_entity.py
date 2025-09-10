import pytest

from src.domain.dataset import Dataset


def test_dataset_requires_id():
    """Test that Dataset requires an id parameter."""
    with pytest.raises(ValueError, match="Dataset id is required"):
        Dataset(id="", name="test_dataset", sequences=["ACGT"])

    # Valid dataset with id
    ds = Dataset(id="test_id", name="test_dataset", sequences=["ACGT"])
    assert ds.id == "test_id"


def test_dataset_basic_stats_and_properties():
    ds = Dataset(
        id="test_dataset_id", name="test_dataset", sequences=["ACGT", "AGGT", "ACCT"]
    )  # infer alphabet
    assert ds.size == 3
    assert ds.max_length == 4
    assert ds.min_length == 4
    assert set(ds.alphabet) == set("ACGT")

    stats = ds.get_statistics()
    assert stats["size"] == 3
    assert stats["min_length"] == 4
    assert stats["max_length"] == 4
    assert stats["alphabet_size"] == 4
    assert stats["total_characters"] == 12


def test_dataset_add_remove_operations():
    ds = Dataset(
        id="test_dataset_id",
        name="test_dataset",
        sequences=["AAAA", "AAAT", "AATA", "ATAA", "TAAA"],
        alphabet="AT",
    )
    ds.add_sequence("TTTT")
    assert ds.size == 6
    with pytest.raises(IndexError):
        ds.remove_sequence(99)
    removed = ds.remove_sequence(0)
    assert removed == "AAAA"
    assert ds.size == 5


def test_dataset_filter_and_non_uniform_length():
    ds = Dataset(
        id="test_dataset_id", name="test_dataset", sequences=["ACG", "ATG", "AAG", "AC"]
    )
    filtered = ds.filter_by_pattern("C", 1)
    assert filtered.get_sequences() == ["ACG", "AC"]
    # non-uniform -> max length is 3
    assert ds.max_length == 3
    assert ds.min_length == 2


def test_dataset_validation_with_alphabet():
    with pytest.raises(ValueError):
        Dataset(
            id="test_dataset_id", name="test_dataset", sequences=["AX"], alphabet="ACGT"
        )
    # valid when alphabet is None (inferred)
    ds = Dataset(
        id="test_dataset_id", name="test_dataset", sequences=["AX"]
    )  # X becomes part of inferred alphabet
    assert "X" in ds.alphabet


def test_dataset_requires_id():
    """Test that Dataset requires an id parameter."""
    with pytest.raises(ValueError, match="Dataset id is required"):
        Dataset(id="", name="test_dataset", sequences=["ACGT"])

    # Valid dataset with id
    ds = Dataset(id="test_id", name="test_dataset", sequences=["ACGT"])
    assert ds.id == "test_id"
