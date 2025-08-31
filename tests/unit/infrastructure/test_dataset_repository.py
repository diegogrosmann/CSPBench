"""
Unit tests for FileDatasetRepository.

Tests file-based persistence for datasets including loading,
saving, and various file format handling scenarios.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from src.domain.dataset import Dataset
from src.domain.errors import DatasetNotFoundError, DatasetValidationError
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository


class TestFileDatasetRepositoryBasic:
    """Basic tests for file repository operations."""

    def test_save_and_load_dataset(self, tmp_path):
        """Test basic save and load operations."""
        # Create test dataset
        sequences = ["ACGTACGT", "ATGTACGT", "ACGTCCGT"]
        dataset = Dataset(id="test_id", name="test_dataset", sequences=sequences)

        # Save dataset
        saved_path = FileDatasetRepository.save(dataset, "test_dataset", str(tmp_path))

        # Verify file was created
        assert Path(saved_path).exists()
        assert saved_path.endswith("test_dataset.fasta")

        # Load dataset back
        loaded_dataset, params = FileDatasetRepository.load(
            "test_dataset.fasta", str(tmp_path)
        )

        # Verify data integrity
        assert loaded_dataset.sequences == sequences
        assert loaded_dataset.size == 3
        assert loaded_dataset.max_length == 8
        assert "file_path" in params

    def test_save_dataset_creates_directories(self, tmp_path):
        """Test that save creates necessary directories."""
        nested_path = tmp_path / "datasets" / "subfolder"
        sequences = ["ACGT"]
        dataset = Dataset(id="test_id", name="test_dataset", sequences=sequences)

        saved_path = FileDatasetRepository.save(dataset, "test", str(nested_path))

        assert Path(saved_path).exists()
        assert Path(saved_path).parent == nested_path

    def test_load_nonexistent_file(self, tmp_path):
        """Test loading non-existent file raises error."""
        with pytest.raises(
            DatasetNotFoundError, match="Dataset not found: nonexistent"
        ):
            FileDatasetRepository.load("nonexistent", str(tmp_path))

    def test_save_empty_dataset(self, tmp_path):
        """Test saving empty dataset."""
        dataset = Dataset(id="test_id", name="test_dataset", sequences=[])

        saved_path = FileDatasetRepository.save(dataset, "empty", str(tmp_path))

        # File should be created even if empty
        assert Path(saved_path).exists()

        # Content should be empty
        with open(saved_path) as f:
            content = f.read().strip()
            assert content == ""


class TestFileDatasetRepositoryFastaFormat:
    """Tests for FASTA format handling."""

    def test_parse_simple_fasta(self, tmp_path):
        """Test parsing simple FASTA format."""
        fasta_content = """>seq1
ACGTACGT
>seq2
ATGTACGT
>seq3
ACGTCCGT
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        dataset, params = FileDatasetRepository.load(str(fasta_file))

        assert dataset.sequences == ["ACGTACGT", "ATGTACGT", "ACGTCCGT"]
        assert dataset.size == 3

    def test_parse_multiline_fasta(self, tmp_path):
        """Test parsing multi-line FASTA sequences."""
        fasta_content = """>seq1
ACGTACGT
GGCCTTAA
>seq2
ATGTACGT
CCGGAATT
"""
        fasta_file = tmp_path / "multiline.fasta"
        fasta_file.write_text(fasta_content)

        dataset, params = FileDatasetRepository.load(str(fasta_file))

        expected_sequences = ["ACGTACGTGGCCTTAA", "ATGTACGTCCGGAATT"]
        assert dataset.sequences == expected_sequences

    def test_parse_fasta_with_empty_lines(self, tmp_path):
        """Test parsing FASTA with empty lines."""
        fasta_content = """>seq1
ACGTACGT

>seq2

ATGTACGT

>seq3
GCGTACGT
"""
        fasta_file = tmp_path / "empty_lines.fasta"
        fasta_file.write_text(fasta_content)

        dataset, params = FileDatasetRepository.load(str(fasta_file))

        assert dataset.sequences == ["ACGTACGT", "ATGTACGT", "GCGTACGT"]

    def test_parse_fasta_no_sequences(self, tmp_path):
        """Test parsing FASTA with no sequences raises error."""
        fasta_content = """>header1
>header2
"""
        fasta_file = tmp_path / "no_seqs.fasta"
        fasta_file.write_text(fasta_content)

        with pytest.raises(DatasetValidationError, match="No sequences found"):
            FileDatasetRepository.load(str(fasta_file))

    def test_parse_fasta_only_headers(self, tmp_path):
        """Test parsing FASTA with only headers (no sequence data)."""
        fasta_content = """>seq1
>seq2
>seq3
"""
        fasta_file = tmp_path / "headers_only.fasta"
        fasta_file.write_text(fasta_content)

        with pytest.raises(DatasetValidationError, match="No sequences found"):
            FileDatasetRepository.load(str(fasta_file))

    def test_save_fasta_format(self, tmp_path):
        """Test that saved FASTA format is correct."""
        sequences = ["ACGTACGT", "ATGTACGT"]
        dataset = Dataset(id="test_id", name="test_dataset", sequences=sequences)

        saved_path = FileDatasetRepository.save(dataset, "format_test", str(tmp_path))

        with open(saved_path) as f:
            content = f.read()

        expected_content = ">seq_0\nACGTACGT\n>seq_1\nATGTACGT\n"
        assert content == expected_content


class TestFileDatasetRepositoryPathHandling:
    """Tests for path resolution and handling."""

    def test_resolve_path_absolute(self, tmp_path):
        """Test resolving absolute paths."""
        absolute_path = tmp_path / "test.fasta"
        absolute_path.write_text(">seq1\nACGT\n")

        resolved = FileDatasetRepository._resolve_path(str(absolute_path), tmp_path)

        assert resolved == absolute_path

    def test_resolve_path_relative_with_extension(self, tmp_path):
        """Test resolving relative path with .fasta extension."""
        resolved = FileDatasetRepository._resolve_path("test.fasta", tmp_path)

        assert resolved == tmp_path / "test.fasta"

    def test_resolve_path_relative_without_extension(self, tmp_path):
        """Test resolving relative path without extension."""
        resolved = FileDatasetRepository._resolve_path("test", tmp_path)

        assert resolved == tmp_path / "test.fasta"

    def test_load_with_absolute_path(self, tmp_path):
        """Test loading using absolute path."""
        fasta_content = ">seq1\nACGT\n"
        absolute_file = tmp_path / "absolute.fasta"
        absolute_file.write_text(fasta_content)

        dataset, params = FileDatasetRepository.load(str(absolute_file))

        assert dataset.sequences == ["ACGT"]

    def test_load_with_relative_path(self, tmp_path):
        """Test loading using relative path."""
        fasta_content = ">seq1\nACGT\n"
        fasta_file = tmp_path / "relative.fasta"
        fasta_file.write_text(fasta_content)

        dataset, params = FileDatasetRepository.load("relative.fasta", str(tmp_path))

        assert dataset.sequences == ["ACGT"]


class TestFileDatasetRepositoryEnvironment:
    """Tests for environment variable handling."""

    def test_get_base_path_default(self, tmp_path):
        """Test getting base path with default."""
        with patch.dict(os.environ, {}, clear=True):
            # Remove DATASET_DIRECTORY if it exists
            if "DATASET_DIRECTORY" in os.environ:
                del os.environ["DATASET_DIRECTORY"]

            base_path = FileDatasetRepository._get_base_path()

            # The implementation returns an absolute path, so we should compare with the resolved path
            expected_path = Path("./datasets").resolve()
            assert base_path == expected_path

    def test_get_base_path_from_env(self, tmp_path):
        """Test getting base path from environment variable."""
        with patch.dict(os.environ, {"DATASET_DIRECTORY": str(tmp_path)}):
            base_path = FileDatasetRepository._get_base_path()

            assert base_path == tmp_path

    def test_get_base_path_explicit(self, tmp_path):
        """Test getting base path with explicit parameter."""
        base_path = FileDatasetRepository._get_base_path(str(tmp_path))

        assert base_path == tmp_path

    def test_get_base_path_creates_directory(self, tmp_path):
        """Test that _get_base_path creates directory if it doesn't exist."""
        nonexistent_path = tmp_path / "nonexistent"

        base_path = FileDatasetRepository._get_base_path(str(nonexistent_path))

        assert base_path.exists()
        assert base_path.is_dir()


class TestFileDatasetRepositoryListOperations:
    """Tests for listing and checking datasets."""

    def test_list_available_empty(self, tmp_path):
        """Test listing available datasets in empty directory."""
        available = FileDatasetRepository.list_available(str(tmp_path))

        assert available == []

    def test_list_available_with_files(self, tmp_path):
        """Test listing available datasets."""
        # Create test files
        (tmp_path / "dataset1.fasta").write_text(">seq1\nACGT\n")
        (tmp_path / "dataset2.fasta").write_text(">seq2\nATGC\n")
        (tmp_path / "not_fasta.txt").write_text("not a fasta file")

        available = FileDatasetRepository.list_available(str(tmp_path))

        # Should only list .fasta files without extension
        assert set(available) == {"dataset1", "dataset2"}

    def test_exists_true(self, tmp_path):
        """Test exists returns True for existing file."""
        fasta_file = tmp_path / "exists.fasta"
        fasta_file.write_text(">seq1\nACGT\n")

        assert FileDatasetRepository.exists("exists", str(tmp_path)) is True

    def test_exists_false(self, tmp_path):
        """Test exists returns False for non-existing file."""
        assert FileDatasetRepository.exists("nonexistent", str(tmp_path)) is False

    def test_delete_existing(self, tmp_path):
        """Test deleting existing dataset."""
        fasta_file = tmp_path / "to_delete.fasta"
        fasta_file.write_text(">seq1\nACGT\n")

        result = FileDatasetRepository.delete("to_delete", str(tmp_path))

        assert result is True
        assert not fasta_file.exists()

    def test_delete_nonexistent(self, tmp_path):
        """Test deleting non-existent dataset."""
        result = FileDatasetRepository.delete("nonexistent", str(tmp_path))

        assert result is False


class TestFileDatasetRepositoryEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_parse_fasta_empty_file(self, tmp_path):
        """Test parsing completely empty FASTA file."""
        empty_file = tmp_path / "empty.fasta"
        empty_file.write_text("")

        with pytest.raises(DatasetValidationError, match="No sequences found"):
            FileDatasetRepository.load(str(empty_file))

    def test_parse_fasta_only_whitespace(self, tmp_path):
        """Test parsing FASTA with only whitespace."""
        whitespace_file = tmp_path / "whitespace.fasta"
        whitespace_file.write_text("   \n\n   \t\t\n   ")

        with pytest.raises(DatasetValidationError, match="No sequences found"):
            FileDatasetRepository.load(str(whitespace_file))

    def test_parse_fasta_with_special_characters(self, tmp_path):
        """Test parsing FASTA with special characters in sequences."""
        fasta_content = """>seq1
ACGT-N*X
>seq2
WWSSKKM
"""
        fasta_file = tmp_path / "special.fasta"
        fasta_file.write_text(fasta_content)

        dataset, params = FileDatasetRepository.load(str(fasta_file))

        assert dataset.sequences == ["ACGT-N*X", "WWSSKKM"]

    def test_save_dataset_with_special_sequences(self, tmp_path):
        """Test saving dataset with special characters."""
        sequences = ["ACGT-N", "WSKMR*"]
        dataset = Dataset(id="test_id", name="test_dataset", sequences=sequences)

        saved_path = FileDatasetRepository.save(dataset, "special", str(tmp_path))
        loaded_dataset, params = FileDatasetRepository.load("special", str(tmp_path))

        assert loaded_dataset.sequences == sequences

    def test_concurrent_access(self, tmp_path):
        """Test concurrent save/load operations."""
        sequences1 = ["ACGT"]
        sequences2 = ["TGCA"]
        dataset1 = Dataset(id="test_id", name="test_dataset", sequences=sequences1)
        dataset2 = Dataset(id="test_id", name="test_dataset", sequences=sequences2)

        # Save multiple datasets
        path1 = FileDatasetRepository.save(dataset1, "concurrent1", str(tmp_path))
        path2 = FileDatasetRepository.save(dataset2, "concurrent2", str(tmp_path))

        # Load both back
        loaded1, _ = FileDatasetRepository.load("concurrent1", str(tmp_path))
        loaded2, _ = FileDatasetRepository.load("concurrent2", str(tmp_path))

        assert loaded1.sequences == sequences1
        assert loaded2.sequences == sequences2

    def test_unicode_file_handling(self, tmp_path):
        """Test handling files with unicode content."""
        fasta_content = """>séquence1
ACGT
>séquence2 with ñoñ-ASCII
TGCA
"""
        fasta_file = tmp_path / "unicode.fasta"
        fasta_file.write_text(fasta_content, encoding="utf-8")

        dataset, params = FileDatasetRepository.load(str(fasta_file))

        assert dataset.sequences == ["ACGT", "TGCA"]


class TestFileDatasetRepositoryBackwardCompatibility:
    """Tests for backward compatibility features."""

    def test_fasta_dataset_repository_alias(self):
        """Test that FastaDatasetRepository alias works."""
        from src.infrastructure.persistence.dataset_repository import (
            FastaDatasetRepository,
        )

        # Should be the same class
        assert FastaDatasetRepository is FileDatasetRepository

    def test_legacy_import_compatibility(self, tmp_path):
        """Test that legacy imports still work."""
        from src.infrastructure.persistence.dataset_repository import (
            FastaDatasetRepository,
        )

        sequences = ["ACGT"]
        dataset = Dataset(id="test_id", name="test_dataset", sequences=sequences)

        # Should work with alias
        saved_path = FastaDatasetRepository.save(dataset, "legacy", str(tmp_path))
        loaded_dataset, params = FastaDatasetRepository.load("legacy", str(tmp_path))

        assert loaded_dataset.sequences == sequences
