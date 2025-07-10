"""
CSPBench - Testes para dataset de arquivo.

Este módulo contém testes para as funcionalidades relacionadas ao
carregamento de datasets de arquivos FASTA.
"""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from src.datasets.dataset_file import load_dataset, load_dataset_with_params
from src.datasets.dataset_utils import ensure_datasets_folder


class TestLoadDataset:
    """Testes para load_dataset."""

    def test_load_dataset_with_existing_file(self):
        """Testa carregamento de dataset com arquivo existente."""
        # Cria um arquivo FASTA na pasta saved_datasets
        datasets_dir = ensure_datasets_folder()
        test_file = datasets_dir / "test_sequences.fasta"

        with open(test_file, "w") as f:
            f.write(">seq1\nATCG\n>seq2\nGCTA\n")

        try:
            # Mock do input para retornar o arquivo criado
            with patch("builtins.input", return_value=str(test_file)):
                sequences, metadata = load_dataset(silent=False)
                assert sequences is not None
                assert len(sequences) == 2
                assert "ATCG" in sequences
                assert "GCTA" in sequences
        finally:
            # Limpa o arquivo de teste
            if test_file.exists():
                test_file.unlink()

    def test_load_dataset_silent(self):
        """Testa carregamento silencioso."""
        # Cria um arquivo com o nome padrão
        datasets_dir = ensure_datasets_folder()
        default_file = datasets_dir / "sequences.fasta"

        with open(default_file, "w") as f:
            f.write(">seq1\nATCG\n")

        try:
            sequences, metadata = load_dataset(silent=True)
            assert sequences is not None
            assert len(sequences) == 1
            assert "ATCG" in sequences
        finally:
            # Limpa o arquivo de teste
            if default_file.exists():
                default_file.unlink()

    def test_load_dataset_nonexistent_file(self):
        """Testa carregamento com arquivo inexistente."""
        with patch("builtins.input", return_value="/nonexistent/file.fasta"):
            with pytest.raises(FileNotFoundError):
                load_dataset(silent=False)


class TestLoadDatasetWithParams:
    """Testes para load_dataset_with_params."""

    def test_load_dataset_with_params_basic(self):
        """Testa carregamento básico com parâmetros."""
        # Cria um arquivo FASTA temporário
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTAA\n")
            temp_file = f.name

        try:
            params = {
                "filepath": temp_file,
                "n_sequences": 2,
            }

            sequences, metadata = load_dataset_with_params(params)
            assert sequences is not None
            assert len(sequences) == 2
            assert metadata is not None
            assert metadata["n_sequences"] == 2
        finally:
            os.unlink(temp_file)

    def test_load_dataset_with_params_all_sequences(self):
        """Testa carregamento de todas as sequências."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTAA\n")
            temp_file = f.name

        try:
            params = {
                "filepath": temp_file,
                "n_sequences": -1,  # Todas as sequências
            }

            sequences, metadata = load_dataset_with_params(params)
            assert sequences is not None
            assert len(sequences) == 3
            assert metadata is not None
            assert metadata["n_sequences"] == 3
        finally:
            os.unlink(temp_file)

    def test_load_dataset_with_params_invalid_file(self):
        """Testa parâmetros com arquivo inválido."""
        params = {
            "filepath": "/nonexistent/file.fasta",
            "n_sequences": 10,
        }

        with pytest.raises(FileNotFoundError):
            load_dataset_with_params(params)

    def test_load_dataset_with_params_empty_file(self):
        """Testa carregamento de arquivo vazio."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write("")  # Arquivo vazio
            temp_file = f.name

        try:
            params = {
                "filepath": temp_file,
                "n_sequences": 5,
            }

            # Arquivo vazio deve gerar exceção
            with pytest.raises(ValueError):
                load_dataset_with_params(params)
        finally:
            os.unlink(temp_file)

    def test_load_dataset_with_params_malformed_fasta(self):
        """Testa carregamento de arquivo FASTA malformado."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write("ATCG\nGCTA\n")  # Sem cabeçalhos
            temp_file = f.name

        try:
            params = {
                "filepath": temp_file,
                "n_sequences": 2,
            }

            # Deve processar mesmo com formato incorreto
            sequences, metadata = load_dataset_with_params(params)
            assert sequences is not None
            assert metadata is not None
        finally:
            os.unlink(temp_file)

    def test_load_dataset_with_params_different_lengths(self):
        """Testa carregamento com sequências de comprimentos diferentes."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\nATCG\n>seq2\nGCTAAAA\n")  # Comprimentos diferentes
            temp_file = f.name

        try:
            params = {
                "filepath": temp_file,
                "n_sequences": 2,
            }

            sequences, metadata = load_dataset_with_params(params)
            assert sequences is not None
            assert (
                len(sequences) == 1
            )  # Uma sequência foi filtrada por comprimento diferente
            assert metadata is not None
        finally:
            os.unlink(temp_file)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
