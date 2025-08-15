from unittest.mock import patch, Mock

from src.application.services.dataset_service import load_dataset
from src.domain.config import (
    SyntheticDatasetConfig,
    FileDatasetConfig,
    EntrezDatasetConfig,
)
from src.domain.dataset import Dataset


def test_load_dataset_synthetic_v2(tmp_path):
    cfg = SyntheticDatasetConfig(
        id="s1", name="Synthetic 1", mode="random", n=3, L=5, alphabet="ACGT", seed=123
    )
    with patch(
        "src.application.services.dataset_service.SyntheticDatasetGenerator.generate_from_config"
    ) as gen:
        gen.return_value = (
            Dataset(["AAAAA", "AAAAT", "AAATA"], alphabet="AT"),
            {"mode": "random"},
        )
        ds, params = load_dataset(cfg)
        assert isinstance(ds, Dataset)
        assert ds.size == 3
        assert params.get("mode") == "random"


def test_load_dataset_file_v2(tmp_path):
    fasta = tmp_path / "d.fasta"
    fasta.write_text(">s1\nACGT\n>s2\nATGT\n")
    cfg = FileDatasetConfig(id="f1", name="File 1", filename=str(fasta))
    with patch(
        "src.application.services.dataset_service.FileDatasetRepository.load"
    ) as loader:
        loader.return_value = (Dataset(["ACGT", "ATGT"]), {"file_path": str(fasta)})
        ds, params = load_dataset(cfg)
        assert ds.size == 2
        assert "file_path" in params


def test_load_dataset_entrez_v2():
    cfg = EntrezDatasetConfig(
        id="e1", name="Entrez 1", query="cox1", db="nucleotide", retmax=2
    )
    with patch(
        "src.application.services.dataset_service.EntrezDatasetDownloader.download"
    ) as dl:
        dl.return_value = (Dataset(["ACGT", "ATGT"]), {"term": "cox1", "n": 2})
        ds, params = load_dataset(cfg)
        assert ds.size == 2
        assert params.get("n") == 2
