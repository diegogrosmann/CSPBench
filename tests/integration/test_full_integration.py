"""
Testes de Integração - CSPBench

Testa a integração completa entre módulos principais do sistema
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from src.datasets.dataset_file import load_dataset_with_params
from src.datasets.dataset_synthetic import (
    generate_dataset,
    generate_dataset_from_params,
)
from src.datasets.dataset_utils import ensure_datasets_folder, save_dataset_fasta
from src.utils.config import merge_configs, validate_config
from src.utils.distance import hamming_distance, max_distance


class TestDatasetIntegration:
    """Testa integração entre módulos de dataset."""

    def test_synthetic_to_file_roundtrip(self):
        """Testa geração sintética -> salvamento -> carregamento."""
        # Gera dataset sintético
        sequences, metadata = generate_dataset_from_params(
            n=5, L=10, alphabet="ATCG", noise=0.1, seed=42
        )

        assert len(sequences) == 5
        assert all(len(seq) == 10 for seq in sequences)

        # Salva em arquivo
        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "test_dataset.fasta")

            success = save_dataset_fasta(sequences, filepath)
            assert success is not None
            assert isinstance(success, Path)
            assert success.exists()

            # Carrega de volta
            params = {"filepath": filepath}
            loaded_sequences, loaded_metadata = load_dataset_with_params(params)

            assert len(loaded_sequences) == len(sequences)
            assert loaded_sequences == sequences
            assert loaded_metadata["filepath"] == filepath

    def test_dataset_folder_creation(self):
        """Testa criação e uso da pasta de datasets."""
        # Garante que pasta existe
        datasets_folder = ensure_datasets_folder()

        assert isinstance(datasets_folder, Path)
        assert datasets_folder.exists()
        assert datasets_folder.is_dir()

        # Testa salvamento na pasta
        sequences = ["ATCG", "GCTA", "TTAA"]

        success = save_dataset_fasta(sequences, "integration_test.fasta")
        assert success is not None
        assert isinstance(success, Path)
        assert success.exists()

        # Limpa arquivo de teste
        success.unlink()

    def test_dataset_generation_with_analysis(self):
        """Testa geração de dataset com análise de distâncias."""
        # Gera dataset pequeno
        sequences, metadata = generate_dataset_from_params(
            n=3, L=8, alphabet="ATCG", noise=0.2, seed=123
        )

        assert len(sequences) == 3

        # Calcula distâncias entre todas as sequências
        distances = []
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                dist = hamming_distance(sequences[i], sequences[j])
                distances.append(dist)

        assert len(distances) == 3  # C(3,2) = 3 pares
        assert all(isinstance(d, int) for d in distances)
        assert all(d >= 0 for d in distances)

        # Calcula distância máxima
        max_dist = max_distance(sequences[0], sequences)
        assert isinstance(max_dist, int)
        assert max_dist >= 0
        assert max_dist <= 8  # Não pode ser maior que o comprimento


class TestConfigIntegration:
    """Testa integração com sistema de configuração."""

    def test_config_validation_with_dataset_params(self):
        """Testa validação de configuração com parâmetros de dataset."""
        # Configuração válida
        config = {
            "algorithm": "BLF-GA",
            "dataset": {"type": "synthetic", "n": 10, "L": 20, "noise": 0.1},
            "optimization": {"n_trials": 50},
        }

        assert validate_config(config) is True

        # Configuração inválida (chaves faltando)
        invalid_config = {"dataset": {}}

        assert (
            validate_config(invalid_config, required_keys=["algorithm", "optimization"])
            is False
        )

    def test_config_merge_with_dataset_override(self):
        """Testa merge de configurações com override de dataset."""
        base_config = {
            "algorithm": "BLF-GA",
            "dataset": {"type": "synthetic", "n": 10, "L": 20},
        }

        override_config = {
            "dataset": {"n": 20, "noise": 0.2},
            "optimization": {"n_trials": 100},
        }

        merged = merge_configs(base_config, override_config)

        assert merged["algorithm"] == "BLF-GA"
        assert merged["dataset"]["n"] == 20  # Sobrescrito
        assert merged["dataset"]["noise"] == 0.2  # Novo
        assert merged["optimization"]["n_trials"] == 100


class TestDistanceAnalysisIntegration:
    """Testa integração com análise de distâncias."""

    def test_distance_analysis_workflow(self):
        """Testa fluxo completo de análise de distâncias."""
        # Gera dataset controlado
        sequences = ["ATCG", "ATCG", "GCTA", "TTTT"]

        # Calcula distâncias ponto a ponto
        distances = []
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                dist = hamming_distance(sequences[i], sequences[j])
                distances.append((i, j, dist))

        # Verifica resultados esperados
        assert len(distances) == 6  # C(4,2) = 6 pares

        # Distância entre sequências idênticas deve ser 0
        identical_pairs = [d for d in distances if d[2] == 0]
        assert len(identical_pairs) == 1  # Apenas as duas primeiras são idênticas

        # Calcula distância máxima
        max_dist = max_distance(sequences[0], sequences)
        assert max_dist == 4  # Distância máxima entre "ATCG" e "TTTT"

    def test_distance_with_different_lengths(self):
        """Testa cálculo de distância com sequências de comprimentos diferentes."""
        # Sequências de comprimentos diferentes
        seq1 = "ATCG"
        seq2 = "ATCGATCG"

        # Deve tratar sequências de comprimentos diferentes
        try:
            dist = hamming_distance(seq1, seq2)
            assert isinstance(dist, int)
            assert dist >= 0
        except (ValueError, TypeError):
            # Aceita erro se função não suporta comprimentos diferentes
            pass

    def test_distance_edge_cases(self):
        """Testa casos extremos de cálculo de distância."""
        # Sequências vazias
        assert hamming_distance("", "") == 0

        # Sequências idênticas
        assert hamming_distance("ATCG", "ATCG") == 0

        # Sequências completamente diferentes
        dist = hamming_distance("AAAA", "TTTT")
        assert dist == 4

        # Teste com lista vazia
        assert max_distance("", []) == 0

        # Teste com uma sequência
        assert max_distance("ATCG", ["ATCG"]) == 0


class TestFullWorkflowIntegration:
    """Testa workflows completos do sistema."""

    def test_complete_synthetic_workflow(self):
        """Testa workflow completo com dataset sintético."""
        # 1. Gera dataset sintético
        sequences, metadata = generate_dataset_from_params(
            n=5, L=8, alphabet="ATCG", noise=0.1, seed=42
        )

        # 2. Valida dataset gerado
        assert len(sequences) == 5
        assert all(len(seq) == 8 for seq in sequences)
        assert all(c in "ATCG" for seq in sequences for c in seq)

        # 3. Salva dataset
        with tempfile.TemporaryDirectory() as temp_dir:
            filepath = os.path.join(temp_dir, "workflow_test.fasta")
            success = save_dataset_fasta(sequences, filepath)
            assert success is not None
            assert isinstance(success, Path)
            assert success.exists()

            # 4. Carrega dataset
            params = {"filepath": filepath}
            loaded_sequences, loaded_metadata = load_dataset_with_params(params)

            # 5. Verifica integridade
            assert loaded_sequences == sequences

            # 6. Analisa distâncias
            max_dist = max_distance(loaded_sequences[0], loaded_sequences)
            assert isinstance(max_dist, int)
            assert max_dist >= 0
            assert max_dist <= 8


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
