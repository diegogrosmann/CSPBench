"""
Unit tests for EntrezDatasetDownloader.

Tests NCBI dataset download functionality including mocking
of external API calls and error handling.
"""

import os
from unittest.mock import Mock, patch, MagicMock
from urllib.error import HTTPError

import pytest

from src.domain.config import EntrezDatasetConfig
from src.domain.dataset import Dataset
from src.infrastructure.external.dataset_entrez import (
    EntrezDatasetDownloader,
    ENTREZ_DEFAULTS,
)


class MockEntrezResult:
    """Mock for Entrez search result."""

    def __init__(self, ids):
        self.data = {"IdList": ids}

    def get(self, key, default=None):
        return self.data.get(key, default)


class MockSeqRecord:
    """Mock for BioPython SeqRecord."""

    def __init__(self, sequence):
        self.seq = sequence


class TestEntrezDatasetDownloaderBasic:
    """Basic tests for Entrez dataset downloader."""

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_download_with_dict_params(self, mock_entrez):
        """Test download with dictionary parameters."""
        # Mock Entrez search
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456", "789"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        # Mock Entrez fetch
        mock_fetch_handle = Mock()
        mock_records = [
            MockSeqRecord("ACGTACGT"),
            MockSeqRecord("ATGTACGT"),
            MockSeqRecord("GCGTACGT"),
        ]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            params = {
                "term": "COX1 AND mitochondrion",
                "db": "nucleotide",
                "n": 3,
                "email": "test@example.com",
            }

            dataset, used_params = EntrezDatasetDownloader.download(params)

            assert isinstance(dataset, Dataset)
            assert dataset.size == 3
            assert dataset.sequences == ["ACGTACGT", "ATGTACGT", "GCGTACGT"]
            assert isinstance(used_params, dict)

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_download_with_config_object(self, mock_entrez):
        """Test download with EntrezDatasetConfig object."""
        # Mock Entrez operations
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        mock_records = [MockSeqRecord("ACGTACGT"), MockSeqRecord("ATGTACGT")]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            config = EntrezDatasetConfig(query="test query", retmax=2, db="nucleotide")

            dataset, used_params = EntrezDatasetDownloader.download(config)

            assert dataset.size == 2
            assert "term" in used_params
            assert used_params["term"] == "test query"

    def test_config_to_params_conversion(self):
        """Test conversion from EntrezDatasetConfig to parameters."""
        config = EntrezDatasetConfig(
            query="COX1",
            retmax=10,
            db="protein",
            min_length=100,
            max_length=500,
            uniform_policy="pad",
        )

        params = EntrezDatasetDownloader._config_to_params(config)

        assert params["term"] == "COX1"
        assert params["n"] == 10
        assert params["db"] == "protein"
        assert params["min_length"] == 100
        assert params["max_length"] == 500
        assert params["uniform_policy"] == "pad"
        # Should include defaults
        assert "email" in params
        assert "rettype" in params


class TestEntrezDatasetDownloaderErrorHandling:
    """Tests for error handling scenarios."""

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_http_error_handling(self, mock_entrez):
        """Test handling of HTTP errors from NCBI."""
        mock_entrez.esearch.side_effect = HTTPError(
            url="", code=429, msg="Too Many Requests", hdrs={}, fp=None
        )

        params = {"term": "test", "n": 5}

        with pytest.raises(ValueError, match="Erro ao acessar NCBI: HTTP 429"):
            EntrezDatasetDownloader.download(params)

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_general_exception_handling(self, mock_entrez):
        """Test handling of general exceptions."""
        mock_entrez.esearch.side_effect = Exception("Connection timeout")

        params = {"term": "test", "n": 5}

        with pytest.raises(ValueError, match="Erro na busca no NCBI"):
            EntrezDatasetDownloader.download(params)

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_no_results_error(self, mock_entrez):
        """Test error when no results are found."""
        mock_search_handle = Mock()
        mock_search_result = {"IdList": []}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.return_value = mock_search_result

        params = {"term": "nonexistent_query", "n": 5}

        with pytest.raises(ValueError, match="Nenhum resultado encontrado"):
            EntrezDatasetDownloader.download(params)

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_invalid_search_result_type(self, mock_entrez):
        """Test error when search result is not a dict."""
        mock_search_handle = Mock()
        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.return_value = "invalid_result"  # Not a dict

        params = {"term": "test", "n": 5}

        with pytest.raises(
            TypeError, match="Resultado da busca Entrez não é um dicionário"
        ):
            EntrezDatasetDownloader.download(params)


class TestEntrezDatasetDownloaderSequenceHandling:
    """Tests for sequence processing and filtering."""

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_fewer_results_than_requested(self, mock_entrez):
        """Test when fewer results are found than requested."""
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456"]}  # Only 2 IDs

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        mock_records = [MockSeqRecord("ACGTACGT"), MockSeqRecord("ATGTACGT")]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            params = {"term": "rare_sequence", "n": 5}  # Request 5 but only 2 available

            dataset, used_params = EntrezDatasetDownloader.download(params)

            assert dataset.size == 2  # Should get 2, not 5

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_length_filtering(self, mock_entrez):
        """Test sequence length filtering."""
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456", "789"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        # Different length sequences
        mock_records = [
            MockSeqRecord("ACGT"),  # length 4
            MockSeqRecord("ACGTACGTACGT"),  # length 12
            MockSeqRecord("ACGTACGT"),  # length 8
        ]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            params = {"term": "test", "n": 3, "min_length": 6, "max_length": 10}

            dataset, used_params = EntrezDatasetDownloader.download(params)

            # Should only include sequence of length 8 (within 6-10 range)
            assert dataset.size == 1
            assert dataset.sequences == ["ACGTACGT"]

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_uniform_policy_strict(self, mock_entrez):
        """Test strict uniform policy filtering."""
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456", "789"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        mock_records = [
            MockSeqRecord("ACGTACGT"),  # length 8
            MockSeqRecord("ACGTACGTAA"),  # length 10
            MockSeqRecord("ACGTACGA"),  # length 8
        ]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            params = {"term": "test", "n": 3, "uniform_policy": "strict"}

            dataset, used_params = EntrezDatasetDownloader.download(params)

            # Should only include sequences of the same length
            assert all(
                len(seq) == len(dataset.sequences[0]) for seq in dataset.sequences
            )

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_uniform_policy_pad(self, mock_entrez):
        """Test padding uniform policy."""
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        mock_records = [
            MockSeqRecord("ACGT"),  # length 4
            MockSeqRecord("ACGTACGT"),  # length 8
        ]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            params = {"term": "test", "n": 2, "uniform_policy": "pad", "pad_char": "N"}

            dataset, used_params = EntrezDatasetDownloader.download(params)

            # All sequences should be padded to same length
            assert len(set(len(seq) for seq in dataset.sequences)) == 1
            assert "N" in dataset.sequences[0]  # First sequence should be padded

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_uniform_policy_truncate(self, mock_entrez):
        """Test truncation uniform policy."""
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        mock_records = [
            MockSeqRecord("ACGT"),  # length 4
            MockSeqRecord("ACGTACGTAA"),  # length 10
        ]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            params = {"term": "test", "n": 2, "uniform_policy": "truncate"}

            dataset, used_params = EntrezDatasetDownloader.download(params)

            # All sequences should be truncated to same length (shortest)
            assert all(len(seq) == 4 for seq in dataset.sequences)


class TestEntrezDatasetDownloaderDefaults:
    """Tests for default value handling."""

    def test_entrez_defaults_exist(self):
        """Test that ENTREZ_DEFAULTS are properly defined."""
        assert "email" in ENTREZ_DEFAULTS
        assert "db" in ENTREZ_DEFAULTS
        assert "n" in ENTREZ_DEFAULTS
        assert "rettype" in ENTREZ_DEFAULTS
        assert "retmode" in ENTREZ_DEFAULTS

    @patch.dict(os.environ, {"NCBI_EMAIL": "test@env.com", "NCBI_API_KEY": "test_key"})
    def test_environment_defaults(self):
        """Test that defaults are loaded from environment."""
        # Reload module to pick up env vars
        import importlib
        from src.infrastructure.external import dataset_entrez

        importlib.reload(dataset_entrez)

        assert dataset_entrez.ENTREZ_DEFAULTS["email"] == "test@env.com"
        assert dataset_entrez.ENTREZ_DEFAULTS["api_key"] == "test_key"

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_parameter_merging_with_defaults(self, mock_entrez):
        """Test that parameters are correctly merged with defaults."""
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        mock_records = [MockSeqRecord("ACGT")]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            # Only provide minimal params
            params = {"term": "test"}

            dataset, used_params = EntrezDatasetDownloader.download(params)

            # Should have merged with defaults
            assert used_params["db"] == ENTREZ_DEFAULTS["db"]
            assert used_params["n"] == ENTREZ_DEFAULTS["n"]
            assert used_params["email"] == ENTREZ_DEFAULTS["email"]


class TestEntrezDatasetDownloaderIntegration:
    """Integration-style tests (still mocked but more comprehensive)."""

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_complete_download_workflow(self, mock_entrez):
        """Test complete download workflow from search to dataset."""
        # Setup mocks for full workflow
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456", "789", "101", "102"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()
        mock_records = [
            MockSeqRecord("ACGTACGTACGT"),  # 12bp
            MockSeqRecord("ATGTACGTACGT"),  # 12bp
            MockSeqRecord("GCGTACGTACGT"),  # 12bp
        ]

        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            params = {
                "term": "COX1 AND mitochondrion",
                "db": "nucleotide",
                "n": 3,
                "email": "researcher@university.edu",
                "min_length": 10,
                "max_length": 15,
                "uniform_policy": "strict",
            }

            dataset, used_params = EntrezDatasetDownloader.download(params)

            # Verify final dataset
            assert isinstance(dataset, Dataset)
            assert dataset.size == 3
            assert all(len(seq) == 12 for seq in dataset.sequences)
            assert set(dataset.alphabet).issubset(set("ACGT"))

            # Verify API was called correctly
            mock_entrez.esearch.assert_called_once()
            mock_entrez.efetch.assert_called_once()

            # Verify used parameters
            assert used_params["term"] == "COX1 AND mitochondrion"
            assert used_params["n"] == 3
