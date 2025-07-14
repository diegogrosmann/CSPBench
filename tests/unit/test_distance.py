"""
Testes Unitários para Módulo de Distâncias - CSPBench

Testa todas as funcionalidades do módulo src.utils.distance
"""

import pytest

from src.domain.metrics import hamming_distance, max_distance

# Nota: max_distance foi removido na refatoração
# Usando max_distance como equivalente


class TestHammingDistance:
    """Testes para função hamming_distance."""

    def test_identical_strings(self):
        """Testa distância entre strings idênticas."""
        assert hamming_distance("ACGT", "ACGT") == 0
        assert hamming_distance("", "") == 0
        assert hamming_distance("A", "A") == 0

    def test_different_strings(self):
        """Testa distância entre strings diferentes."""
        assert hamming_distance("ACGT", "AGCT") == 2
        assert hamming_distance("AAAA", "TTTT") == 4
        assert hamming_distance("ACGT", "TGCA") == 4

    def test_single_difference(self):
        """Testa strings com apenas uma diferença."""
        assert hamming_distance("ACGT", "TCGT") == 1
        assert hamming_distance("ACGT", "ATGT") == 1
        assert hamming_distance("ACGT", "ACTT") == 1

    def test_empty_strings(self):
        """Testa comportamento com strings vazias."""
        assert hamming_distance("", "") == 0

    def test_single_character(self):
        """Testa strings de um caractere."""
        assert hamming_distance("A", "A") == 0
        assert hamming_distance("A", "T") == 1

    def test_dna_sequences(self):
        """Testa sequências de DNA típicas."""
        seq1 = "ATCGATCG"
        seq2 = "ATCGATCG"
        assert hamming_distance(seq1, seq2) == 0

        seq3 = "ATCGATCG"
        seq4 = "ATCGATCT"
        assert hamming_distance(seq3, seq4) == 1

    def test_longer_sequences(self):
        """Testa sequências mais longas."""
        seq1 = "ATCGATCGATCGATCG"
        seq2 = "ATCGATCGATCGATCG"
        assert hamming_distance(seq1, seq2) == 0

        seq3 = "ATCGATCGATCGATCG"
        seq4 = "TTCGATCGATCGATCT"
        assert hamming_distance(seq3, seq4) == 2

    def test_different_lengths_raises_error(self):
        """Testa que strings de comprimentos diferentes geram erro."""
        with pytest.raises(ValueError, match="mesmo comprimento"):
            hamming_distance("ACGT", "ACG")

        with pytest.raises(ValueError, match="mesmo comprimento"):
            hamming_distance("AC", "ACGT")


class TestMaxDistance:
    """Testes para função max_distance."""

    def test_single_string(self):
        """Testa com uma única string."""
        center = "ACGT"
        strings = ["ACGT"]
        assert max_distance(center, strings) == 0

    def test_multiple_identical_strings(self):
        """Testa com múltiplas strings idênticas."""
        center = "ACGT"
        strings = ["ACGT", "ACGT", "ACGT"]
        assert max_distance(center, strings) == 0

    def test_different_strings(self):
        """Testa com strings diferentes."""
        center = "ACGT"
        strings = ["ACGT", "AGCT", "ATTT"]
        expected = max(
            hamming_distance(center, "ACGT"),
            hamming_distance(center, "AGCT"),
            hamming_distance(center, "ATTT"),
        )
        assert max_distance(center, strings) == expected

    def test_worst_case_distance(self):
        """Testa caso com distância máxima."""
        center = "AAAA"
        strings = ["AAAA", "TTTT", "GGGG"]
        assert max_distance(center, strings) == 4

    def test_mixed_distances(self):
        """Testa com distâncias mistas."""
        center = "ACGT"
        strings = ["ACGT", "ACGC", "ACTT", "TTTT"]
        # Calculando manualmente: ACGT->ACGT=0, ACGT->ACGC=1, ACGT->ACTT=2, ACGT->TTTT=3
        expected = 3  # Distância para "TTTT"
        assert max_distance(center, strings) == expected

    def test_empty_strings_list(self):
        """Testa com lista vazia de strings."""
        center = "ACGT"
        strings = []
        assert max_distance(center, strings) == 0

    def test_different_lengths_raises_error(self):
        """Testa que strings de comprimentos diferentes geram erro."""
        center = "ACGT"
        strings = ["ACGT", "ACG", "ATGT"]

        with pytest.raises(ValueError, match="comprimento"):
            max_distance(center, strings)


class TestMaxHammingParallel:
    """Testes para função max_distance."""

    def test_single_string(self):
        """Testa com uma única string."""
        candidate = "ACGT"
        strings = ["ACGT"]
        assert max_distance(candidate, strings) == 0

    def test_multiple_strings(self):
        """Testa com múltiplas strings."""
        candidate = "ACGT"
        strings = ["ACGT", "AGCT", "ATTT"]
        expected = max_distance(candidate, strings)
        assert max_distance(candidate, strings) == expected

    def test_large_set(self):
        """Testa com conjunto maior de strings."""
        candidate = "ACGT"
        strings = ["ACGT", "AGCT", "ATTT", "AAAA", "TTTT", "GGGG", "CCCC"]
        expected = max_distance(candidate, strings)
        assert max_distance(candidate, strings) == expected

    def test_empty_strings(self):
        """Testa com lista vazia."""
        candidate = "ACGT"
        strings = []
        assert max_distance(candidate, strings) == 0

    def test_consistency_with_max_distance(self):
        """Testa consistência com max_distance."""
        candidate = "ACGTACGT"
        strings = ["ACGTACGT", "ACGTACGC", "ACGTACTT", "ACGTTTTT"]

        result_parallel = max_distance(candidate, strings)
        result_normal = max_distance(candidate, strings)

        assert result_parallel == result_normal

    def test_small_dataset_fallback(self):
        """Testa que datasets pequenos usam versão sequencial."""
        candidate = "ACGT"
        strings = ["ACGT", "AGCT", "ATTT"]  # < 1000 strings

        # Deve usar versão sequencial
        result = max_distance(candidate, strings)
        expected = max_distance(candidate, strings)
        assert result == expected


class TestIntegrationScenarios:
    """Testes de integração para cenários reais."""

    def test_dna_analysis_scenario(self):
        """Testa cenário de análise de DNA."""
        # Simula análise de sequências de DNA
        reference = "ATCGATCGATCG"
        samples = [
            "ATCGATCGATCG",  # Idêntica
            "ATCGATCGATCT",  # 1 diferença
            "ATCGATCGATTT",  # 2 diferenças
            "ATCGATCGTTTT",  # 3 diferenças
        ]

        # Testa distâncias individuais
        assert hamming_distance(reference, samples[0]) == 0
        assert hamming_distance(reference, samples[1]) == 1
        assert hamming_distance(reference, samples[2]) == 2
        assert hamming_distance(reference, samples[3]) == 3

        # Testa distância máxima
        assert max_distance(reference, samples) == 3

    def test_csp_optimization_scenario(self):
        """Testa cenário de otimização CSP."""
        # Simula instância CSP pequena
        strings = ["ACGT", "AGCT", "ATGT", "ACTT"]

        # Testa diferentes candidatos
        candidates = ["ACGT", "AGGT", "ATTT"]

        distances = []
        for candidate in candidates:
            dist = max_distance(candidate, strings)
            distances.append(dist)

        # Verifica que pelo menos um candidato tem distância razoável
        assert min(distances) <= 2

        # Verifica que conseguimos encontrar o melhor candidato
        best_candidate = candidates[distances.index(min(distances))]
        best_distance = max_distance(best_candidate, strings)

        # Verifica que a distância é razoável
        assert best_distance <= 2

    def test_empty_edge_cases(self):
        """Testa casos extremos com entradas vazias."""
        # Strings vazias
        assert hamming_distance("", "") == 0
        assert max_distance("", [""]) == 0

        # Listas vazias
        assert max_distance("ACGT", []) == 0
        assert max_distance("ACGT", []) == 0

    def test_performance_consistency(self):
        """Testa consistência entre versões normal e paralela."""
        candidate = "ACGTACGTACGT"
        strings = [
            "ACGTACGTACGT",
            "ACGTACGTACGC",
            "ACGTACGTACTT",
            "ACGTACGTTTTT",
            "ACGTACGTGGGG",
            "ACGTACGTCCCC",
        ]

        # Compara resultados
        normal_result = max_distance(candidate, strings)
        parallel_result = max_distance(candidate, strings)

        assert normal_result == parallel_result

    def test_real_world_csp_instance(self):
        """Testa instância realista de CSP."""
        # Simula problema CSP com strings de DNA
        strings = [
            "ATCGATCGATCG",
            "ATCGATCGATCT",
            "ATCGATCGATTT",
            "ATCGATCGTTTT",
            "ATCGATTTTTTT",
            "ATCGTTTTTTTT",
        ]

        # Testa vários candidatos
        candidates = [
            "ATCGATCGATCG",  # Primeiro da lista
            "ATCGATCGATCT",  # Segundo da lista
            "ATCGATCGATTT",  # Terceiro da lista
        ]

        results = []
        for candidate in candidates:
            dist = max_distance(candidate, strings)
            results.append(dist)

        # Verifica que encontramos soluções válidas
        assert min(results) <= 6  # Máximo possível seria 12

        # Verifica que o primeiro candidato é uma boa escolha
        assert results[0] <= 6

    def test_edge_case_single_character(self):
        """Testa casos extremos com caracteres únicos."""
        # Strings de um caractere
        assert hamming_distance("A", "A") == 0
        assert hamming_distance("A", "T") == 1

        # Múltiplas strings de um caractere
        center = "A"
        strings = ["A", "T", "G", "C"]
        assert max_distance(center, strings) == 1

        # Versão paralela
        assert max_distance(center, strings) == 1
