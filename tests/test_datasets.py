import pytest

from datasets.dataset_file import _load_fasta, _load_text, load_dataset_with_params
from datasets.dataset_utils import ask_save_dataset


def test_load_fasta_and_text(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nACGTAC\n>seq2\nTGCAAA\n")
    seqs = _load_fasta(fasta)
    assert seqs == ["ACGTAC", "TGCAAA"]

    txt = tmp_path / "test.txt"
    txt.write_text("acgtac\ntgcaaa\n#comentario\n\n")
    seqs2 = _load_text(txt)
    assert seqs2 == ["ACGTAC", "TGCAAA"]

    txt2 = tmp_path / "test2.txt"
    txt2.write_text("1234\nACGT\n!!@@\n")
    seqs3 = _load_text(txt2)
    assert seqs3 == ["ACGT"]


def test_load_dataset_with_params(tmp_path):
    fasta = tmp_path / "valid.fasta"
    fasta.write_text(">seq1\nAAAA\n>seq2\nAAAT\n")
    params = {"filepath": str(fasta)}
    seqs, used = load_dataset_with_params(params)
    assert seqs == ["AAAA", "AAAT"]
    assert used["n"] == 2 and used["L"] == 4

    txt = tmp_path / "inv.txt"
    txt.write_text("AAAA\nAAA\n")
    params2 = {"filepath": str(txt)}
    seqs2, used2 = load_dataset_with_params(params2)
    assert seqs2 == ["AAAA"]
    assert used2["n"] == 1 and used2["L"] == 4

    txt2 = tmp_path / "empty.txt"
    txt2.write_text("")
    params3 = {"filepath": str(txt2)}
    with pytest.raises(ValueError):
        load_dataset_with_params(params3)


def test_ask_save_dataset(monkeypatch, tmp_path):
    monkeypatch.setattr("builtins.input", lambda _: "s")
    seqs = ["ACGT", "TGCA"]
    params = {
        "n": 2,
        "L": 4,
        "alphabet": "ACGT",
        "noise": 0.1,
        "term": "test",
        "db": "nucleotide",
    }
    assert ask_save_dataset(seqs, "synthetic", params)
    assert ask_save_dataset(seqs, "entrez", params)
    assert ask_save_dataset(seqs, "custom", params)

    monkeypatch.setattr("builtins.input", lambda _: "n")
    assert not ask_save_dataset(seqs, "synthetic", params)

    monkeypatch.setattr("builtins.input", lambda _: "s")

    def fail_save(*a, **kw):
        raise Exception("fail")

    monkeypatch.setattr("datasets.dataset_utils.save_dataset_fasta", fail_save)
    assert not ask_save_dataset(seqs, "synthetic", params)
