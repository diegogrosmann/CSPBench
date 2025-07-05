import pytest

from src.core.io import load_fasta, load_txt
from src.datasets.dataset_file import load_dataset_with_params
from src.ui.cli.save_wizard import ask_save_dataset


def test_load_fasta_and_text(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nACGTAC\n>seq2\nTGCAAA\n")
    seqs = load_fasta(fasta)
    assert seqs == ["ACGTAC", "TGCAAA"]

    txt = tmp_path / "test.txt"
    txt.write_text("acgtac\ntgcaaa\n#comentario\n\n")
    seqs2 = load_txt(txt)
    assert seqs2 == ["acgtac", "tgcaaa"]  # Mudança: aceitar minúsculas

    # Teste com dados mistos - load_txt pega tudo
    txt2 = tmp_path / "test2.txt"
    txt2.write_text("1234\nACGT\n!!@@\n")
    seqs3 = load_txt(txt2)
    assert seqs3 == ["1234", "ACGT", "!!@@"]  # load_txt carrega tudo

    # Teste com load_dataset_with_params - validação atual só filtra por comprimento
    # Todas têm comprimento 4, então todas passam (mas gera warning)
    params = {"filepath": str(txt2)}
    seqs4, used = load_dataset_with_params(params)
    assert seqs4 == ["1234", "ACGT", "!!@@"]  # Filtro atual é só por comprimento
    assert used["n"] == 3 and used["L"] == 4


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


def test_ask_save_dataset(monkeypatch, tmp_path):  # pylint: disable=unused-argument
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
        raise ValueError("fail")

    monkeypatch.setattr("src.ui.cli.save_wizard.save_dataset_fasta", fail_save)
    assert not ask_save_dataset(seqs, "synthetic", params)
