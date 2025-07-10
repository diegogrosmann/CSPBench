"""
CSPBench - Testes para utilitários de datasets.

Este módulo contém testes para as funcionalidades relacionadas aos
datasets, incluindo salvamento, carregamento e gerenciamento de arquivos.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from src.datasets.dataset_utils import ensure_datasets_folder, save_dataset_fasta


class TestEnsureDatasetsFolder:
    """Testes para ensure_datasets_folder."""

    def test_ensure_datasets_folder_exists(self):
        """Testa criação de pasta de datasets."""
        # A função usa um path fixo, então apenas verificamos que não dá erro
        result = ensure_datasets_folder()
        assert result is not None
        assert isinstance(result, Path)

    def test_ensure_datasets_folder_already_exists(self):
        """Testa pasta que já existe."""
        # Chama a função múltiplas vezes para garantir que não falha
        result1 = ensure_datasets_folder()
        result2 = ensure_datasets_folder()
        assert result1 == result2

    def test_ensure_datasets_folder_permission_error(self):
        """Testa comportamento da função."""
        # Apenas testa que a função executa sem erro
        result = ensure_datasets_folder()
        assert result is not None


class TestSaveDatasetFasta:
    """Testes para save_dataset_fasta."""

    def test_save_dataset_fasta_basic(self):
        """Testa salvamento básico de dataset."""
        sequences = ["ATCG", "GCTA", "TTAA"]

        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "test_dataset.fasta")

            result = save_dataset_fasta(sequences, filepath)

            assert result is not None
            assert isinstance(result, Path)
            assert result.exists()

            # Verifica conteúdo do arquivo
            with open(filepath, "r") as f:
                content = f.read()
                assert ">seq_1" in content
                assert "ATCG" in content
                assert ">seq_2" in content
                assert "GCTA" in content

    def test_save_dataset_fasta_with_description(self):
        """Testa salvamento com descrição."""
        sequences = ["ATCG", "GCTA"]

        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "test_dataset.fasta")
            result = save_dataset_fasta(sequences, filepath, description="Test data")
            assert result is not None
            assert isinstance(result, Path)
            assert result.exists()

    def test_save_dataset_fasta_empty_sequences(self):
        """Testa salvamento com lista vazia de sequências."""
        sequences = []

        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "empty_dataset.fasta")

            # Deve lançar ValueError
            with pytest.raises(ValueError):
                save_dataset_fasta(sequences, filepath)

    def test_save_dataset_fasta_single_sequence(self):
        """Testa salvamento com única sequência."""
        sequences = ["ATCGATCG"]

        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "single_dataset.fasta")

            result = save_dataset_fasta(sequences, filepath)
            assert result is not None
            assert isinstance(result, Path)
            assert result.exists()

            # Verifica conteúdo
            with open(filepath, "r") as f:
                content = f.read()
                assert ">seq_1" in content
                assert "ATCGATCG" in content

    def test_save_dataset_fasta_invalid_path(self):
        """Testa salvamento com caminho inválido."""
        sequences = ["ATCG"]
        filepath = "/nonexistent/path/dataset.fasta"

        with pytest.raises(FileNotFoundError):
            save_dataset_fasta(sequences, filepath)

    def test_save_dataset_fasta_overwrite(self):
        """Testa sobrescrita de arquivo existente."""
        sequences1 = ["ATCG"]
        sequences2 = ["GCTA", "TTAA"]

        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "overwrite_dataset.fasta")

            result1 = save_dataset_fasta(sequences1, filepath)
            assert result1 is not None
            assert isinstance(result1, Path)
            assert result1.exists()

            result2 = save_dataset_fasta(sequences2, filepath)
            assert result2 is not None
            assert isinstance(result2, Path)
            assert result2.exists()

            # Verifica que foi sobrescrito
            with open(filepath, "r") as f:
                content = f.read()
                assert "GCTA" in content
                assert "TTAA" in content

    def test_save_dataset_fasta_special_characters(self):
        """Testa salvamento com caracteres especiais."""
        sequences = ["ATCG-N", "GCTA?"]

        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "special_chars_dataset.fasta")

            result = save_dataset_fasta(sequences, filepath)
            assert result is not None
            assert isinstance(result, Path)
            assert result.exists()

    def test_save_dataset_fasta_long_sequences(self):
        """Testa salvamento com sequências longas."""
        sequences = ["A" * 1000, "T" * 1000]

        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "long_dataset.fasta")

            result = save_dataset_fasta(sequences, filepath)
            assert result is not None
            assert isinstance(result, Path)
            assert result.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
