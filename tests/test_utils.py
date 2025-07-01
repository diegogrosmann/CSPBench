import pytest
from datasets.dataset_synthetic import generate_dataset
from utils.config import SYNTHETIC_DEFAULTS, safe_input
from datasets.dataset_utils import ensure_datasets_folder
import builtins

def test_generate_synthetic_dataset(monkeypatch):
    inputs = iter(["5", "6", "ACGT", "0.2", "n", "42"])
    monkeypatch.setattr(builtins, "input", lambda _: next(inputs))
    seqs, params = generate_dataset()
    assert len(seqs) == 5
    assert all(len(s) == 6 for s in seqs)
    assert set("".join(seqs)).issubset(set("ACGT"))
    assert params["noise"] == 0.2
    assert params["seed"] == 42

def test_generate_synthetic_fully_random(monkeypatch):
    inputs = iter(["3", "4", "AB", "0.0", "s", "123"])
    monkeypatch.setattr(builtins, "input", lambda _: next(inputs))
    seqs, params = generate_dataset()
    assert len(seqs) == 3
    assert all(len(s) == 4 for s in seqs)
    assert set("".join(seqs)).issubset(set("AB"))
    assert params["fully_random"] is True
    assert params["seed"] == 123

def test_synthetic_defaults():
    assert SYNTHETIC_DEFAULTS["n"] > 0
    assert SYNTHETIC_DEFAULTS["L"] > 0
    assert set(SYNTHETIC_DEFAULTS["alphabet"]).issuperset({"A", "C", "G", "T"})

def test_ensure_datasets_folder():
    path = ensure_datasets_folder()
    assert path.exists()
    assert path.is_dir()

def test_safe_input_keyboard_interrupt(monkeypatch):
    def raise_interrupt(prompt, default=""): raise KeyboardInterrupt
    monkeypatch.setattr(builtins, "input", lambda _: (_ for _ in ()).throw(KeyboardInterrupt))
    try:
        safe_input("Teste: ")
    except SystemExit:
        assert True
    else:
        assert False, "safe_input deveria encerrar o programa em KeyboardInterrupt"
