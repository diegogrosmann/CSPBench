from unittest.mock import patch, Mock

from src.infrastructure.external.dataset_entrez import EntrezDatasetDownloader


def test_entrez_downloader_happy_path_v2():
    # Mock Entrez and SeqIO
    with patch("src.infrastructure.external.dataset_entrez.Entrez") as Entrez:
        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as SeqIO:
            Entrez.esearch.return_value = Mock()
            Entrez.read.return_value = {"IdList": ["1", "2", "3"]}
            Entrez.efetch.return_value = Mock()

            class Rec:
                def __init__(self, s):
                    self.seq = s

            SeqIO.parse.return_value = [Rec("ACGT"), Rec("ATGT")]
            ds, used = EntrezDatasetDownloader.download({"term": "q", "n": 2})
            assert ds.size == 2
            assert used["n_obtained"] == 2
