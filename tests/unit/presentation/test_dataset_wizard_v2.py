from unittest.mock import patch

from src.presentation.cli.dataset_wizard import DatasetWizard


def test_dataset_wizard_synthetic_defaults_v2():
    wiz = DatasetWizard()
    # Choose synthetic (1), then accept all defaults: method(1), n, L, alphabet, seed empty
    # inputs: menu choice + method + n + length + alphabet + seed
    inputs = ["1", "", "", "", "", ""]
    with patch("builtins.input", side_effect=inputs):
        with patch("builtins.print"):
            kind = wiz.show_main_menu()
            assert kind == "synthetic"
            params = wiz.collect_synthetic_params()
            assert params["method"] in {"random", "noise", "clustered", "mutations"}
            assert isinstance(params["n"], int)
            assert isinstance(params["length"], int)
            assert isinstance(params["alphabet"], str)


def test_dataset_wizard_real_file_v2():
    wiz = DatasetWizard()
    # Choose real (2), select file (2), provide path
    inputs = ["2", "2", "datasets/example.fasta", ""]
    with patch("builtins.input", side_effect=inputs):
        with patch("builtins.print"):
            kind = wiz.show_main_menu()
            assert kind == "real"
            params = wiz.collect_real_params()
            assert params["source"] == "file"
            assert params["file_path"].endswith(".fasta")


def test_dataset_wizard_default_filename_rules_v2():
    wiz = DatasetWizard()
    sf = wiz.generate_default_filename(
        "synthetic", {"method": "random", "n": 5, "length": 8}
    )
    assert sf.startswith("synthetic_random_n5_L8")
    rf = wiz.generate_default_filename(
        "real", {"source": "ncbi", "query": "gene[ALL] AND species:Human"}
    )
    assert rf.startswith("real_gene_ALL_AND_species")
    # truncation to 30 chars for query part
    long = wiz.generate_default_filename("real", {"source": "ncbi", "query": "a" * 80})
    assert len(long) <= len("real_" + "a" * 30 + ".fasta")
