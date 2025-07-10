"""
Testes Unitários para Módulo de Dataset Sintético - CSPBench

Testa as funcionalidades do módulo src.datasets.dataset_synthetic
"""

import random
from unittest.mock import MagicMock, patch

import pytest

from src.datasets.dataset_synthetic import (
    generate_dataset,
    generate_dataset_from_params,
    generate_dataset_with_params,
    generate_perturbations,
)


class TestGenerateDataset:
    """Testes para função generate_dataset."""

    @patch("src.datasets.dataset_synthetic.safe_input")
    def test_generate_dataset_default(self, mock_input):
        """Testa geração padrão de dataset."""
        # Mock para entrada do usuário
        mock_input.side_effect = ["10", "5", "0.1", "synthetic_test"]

        sequences, metadata = generate_dataset(silent=True)

        assert isinstance(sequences, list)
        assert isinstance(metadata, dict)
        assert len(sequences) > 0
        assert "n" in metadata
        assert "L" in metadata
        assert "noise" in metadata

    def test_generate_dataset_silent(self):
        """Testa geração silenciosa de dataset."""
        sequences, metadata = generate_dataset(silent=True)

        assert isinstance(sequences, list)
        assert isinstance(metadata, dict)
        assert len(sequences) > 0
        assert "n" in metadata
        assert "L" in metadata
        assert "noise" in metadata


class TestGenerateDatasetWithParams:
    """Testes para função generate_dataset_with_params."""

    def test_generate_dataset_with_params_basic(self):
        """Testa geração com parâmetros básicos."""
        params = {"n": 5, "L": 10, "noise": 0.1, "filename": "test_dataset"}

        sequences, metadata = generate_dataset_with_params(params)

        assert isinstance(sequences, list)
        assert isinstance(metadata, dict)
        assert len(sequences) == 5
        assert metadata["n"] == 5
        assert metadata["L"] == 10
        assert metadata["noise"] == 0.1

    def test_generate_dataset_with_params_invalid(self):
        """Testa geração com parâmetros inválidos."""
        params = {"n": 0, "L": 10, "noise": 0.1, "filename": "test_dataset"}  # Inválido

        # n=0 deve gerar 0 sequências (comportamento válido)
        sequences, metadata = generate_dataset_with_params(params)
        assert len(sequences) == 0  # n=0 resulta em 0 sequências
        assert metadata["n"] == 0

    def test_generate_dataset_with_params_missing(self):
        """Testa geração com parâmetros faltando."""
        params = {
            "n": 5,
            # L faltando
            "noise": 0.1,
            "filename": "test_dataset",
        }

        # Deve usar valores padrão ou gerar erro
        try:
            sequences, metadata = generate_dataset_with_params(params)
            assert isinstance(sequences, list)
            assert isinstance(metadata, dict)
        except (ValueError, KeyError):
            # Aceita erro se parâmetros obrigatórios faltam
            pass


class TestGeneratePerturbations:
    """Testes para função generate_perturbations."""

    def test_generate_perturbations_basic(self):
        """Testa geração básica de perturbações."""
        original_sequence = "ATCGATCG"
        rng = random.Random(42)
        perturbations = generate_perturbations(original_sequence, 3, 2, "ATCG", rng)

        assert isinstance(perturbations, list)
        assert len(perturbations) == 3
        for seq in perturbations:
            assert isinstance(seq, str)
            assert len(seq) == len(original_sequence)

    def test_generate_perturbations_zero_distance(self):
        """Testa geração com distância zero."""
        original_sequence = "ATCGATCG"
        rng = random.Random(42)
        perturbations = generate_perturbations(original_sequence, 3, 0, "ATCG", rng)

        assert isinstance(perturbations, list)
        assert len(perturbations) == 3
        # Com distância zero, todas devem ser idênticas ao original
        for seq in perturbations:
            assert seq == original_sequence

    def test_generate_perturbations_max_distance(self):
        """Testa geração com distância máxima."""
        original_sequence = "ATCGATCG"
        rng = random.Random(42)
        max_distance = len(original_sequence)
        perturbations = generate_perturbations(
            original_sequence, 3, max_distance, "ATCG", rng
        )

        assert isinstance(perturbations, list)
        assert len(perturbations) == 3
        for seq in perturbations:
            assert isinstance(seq, str)
            assert len(seq) == len(original_sequence)

    def test_generate_perturbations_empty_sequence(self):
        """Testa geração com sequência vazia."""
        original_sequence = ""
        rng = random.Random(42)
        perturbations = generate_perturbations(original_sequence, 3, 0, "ATCG", rng)

        assert isinstance(perturbations, list)
        assert len(perturbations) == 3
        for seq in perturbations:
            assert seq == ""


class TestGenerateDatasetFromParams:
    """Testes para função generate_dataset_from_params."""

    def test_generate_dataset_from_params_basic(self):
        """Testa geração com parâmetros básicos."""
        sequences, metadata = generate_dataset_from_params(
            n=5, L=10, alphabet="ATCG", noise=0.1
        )

        assert isinstance(sequences, list)
        assert isinstance(metadata, dict)
        assert len(sequences) == 5
        for seq in sequences:
            assert isinstance(seq, str)
            assert len(seq) == 10
            # Verifica se usa apenas caracteres do alfabeto
            assert all(c in "ATCG" for c in seq)

    def test_generate_dataset_from_params_custom_alphabet(self):
        """Testa geração com alfabeto personalizado."""
        sequences, metadata = generate_dataset_from_params(
            n=3, L=8, alphabet="AB", noise=0.2
        )

        assert isinstance(sequences, list)
        assert isinstance(metadata, dict)
        assert len(sequences) == 3
        for seq in sequences:
            assert isinstance(seq, str)
            assert len(seq) == 8
            # Verifica se usa apenas caracteres do alfabeto personalizado
            assert all(c in "AB" for c in seq)

    def test_generate_dataset_from_params_fully_random(self):
        """Testa geração com modo totalmente aleatório."""
        sequences, metadata = generate_dataset_from_params(
            n=5, L=10, alphabet="ATCG", noise=0.1, fully_random=True
        )

        assert isinstance(sequences, list)
        assert isinstance(metadata, dict)
        assert len(sequences) == 5
        assert metadata["fully_random"] is True

    def test_generate_dataset_from_params_with_seed(self):
        """Testa geração com semente específica."""
        sequences1, metadata1 = generate_dataset_from_params(
            n=5, L=10, alphabet="ATCG", noise=0.1, seed=42
        )

        sequences2, metadata2 = generate_dataset_from_params(
            n=5, L=10, alphabet="ATCG", noise=0.1, seed=42
        )

        # Com a mesma semente, resultados devem ser idênticos
        assert sequences1 == sequences2
        assert metadata1["seed"] == metadata2["seed"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
