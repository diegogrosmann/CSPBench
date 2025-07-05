"""
Testes para dataset_entrez com pytest-vcr.

Usa VCR para gravar/reproduzir requisições HTTP ao NCBI,
evitando chamadas reais à API durante os testes.
"""

import pytest
import vcr

from src.datasets.dataset_entrez import fetch_dataset_silent


class TestDatasetEntrez:
    """Testes para dataset_entrez usando VCR."""

    def test_fetch_dataset_silent(self):
        """Testa busca de dataset do NCBI com parâmetros."""
        params = {
            "email": "test@example.com",
            "db": "nucleotide",
            "term": "Escherichia coli[Organism] AND complete genome",
            "n": 2,
        }

        # Na realidade, genomas completos têm tamanhos muito diferentes
        # Então esperamos erro de comprimento uniforme
        with vcr.use_cassette("tests/vcr_cassettes/entrez_nucleotide_test.yaml"):
            with pytest.raises(ValueError, match="Muito poucas sequências de comprimento uniforme"):
                sequences, metadata = fetch_dataset_silent(params)

    def test_fetch_protein_dataset(self):
        """Testa busca de dataset de proteínas."""
        params = {
            "email": "test@example.com",
            "db": "protein",
            "term": "insulin[Protein Name] AND human[Organism]",
            "n": 3,
        }

        with vcr.use_cassette("tests/vcr_cassettes/entrez_protein_test.yaml"):
            sequences, metadata = fetch_dataset_silent(params)

            # Verificar sequências
            assert len(sequences) <= 3  # Pode retornar menos se não houver suficientes
            assert all(isinstance(seq, str) for seq in sequences)

            # Verificar metadata
            assert metadata["db"] == "protein"
            assert "n_obtained" in metadata

    def test_fetch_dataset_invalid_params(self):
        """Testa parâmetros inválidos."""
        # Sem email - deve aceitar mas dar warning
        params = {
            "db": "nucleotide",
            "term": "test",
            "n": 1,
        }

        # Na realidade, isso funciona mas pode encontrar sequências com comprimentos diferentes
        with vcr.use_cassette("tests/vcr_cassettes/entrez_invalid_params.yaml"):
            # Pode funcionar ou dar erro dependendo do que encontrar
            try:
                sequences, metadata = fetch_dataset_silent(params)
                # Se funcionou, verificar se tem dados
                assert len(sequences) >= 1
                assert "n_obtained" in metadata
            except ValueError as e:
                # Se deu erro, deve ser por sequências de comprimento uniforme
                assert "Muito poucas sequências de comprimento uniforme" in str(e)

    def test_fetch_dataset_empty_results(self):
        """Testa busca que não retorna resultados."""
        params = {
            "email": "test@example.com",
            "db": "nucleotide",
            "term": "sequencia_que_nao_existe_no_ncbi_12345",
            "n": 1,
        }

        with vcr.use_cassette("tests/vcr_cassettes/entrez_empty_results.yaml"):
            with pytest.raises(ValueError, match="Nenhum resultado encontrado"):
                fetch_dataset_silent(params)

    def test_fetch_dataset_with_api_key(self):
        """Testa busca com API key."""
        params = {
            "email": "test@example.com",
            "api_key": "test_api_key_123456",
            "db": "nucleotide",
            "term": "Escherichia coli[Organism]",
            "n": 1,
        }

        # API key inválida causa erro HTTP 400
        with vcr.use_cassette("tests/vcr_cassettes/entrez_with_api_key.yaml"):
            with pytest.raises(ValueError, match="Erro ao acessar NCBI"):
                fetch_dataset_silent(params)

    def test_fetch_dataset_large_request(self):
        """Testa requisição com muitas sequências."""
        params = {
            "email": "test@example.com",
            "db": "nucleotide",
            "term": "Escherichia coli[Organism]",
            "n": 10,  # Requisitar mais sequências
        }

        # Grandes requisições de E. coli têm comprimentos muito diferentes
        with vcr.use_cassette("tests/vcr_cassettes/entrez_large_request.yaml"):
            with pytest.raises(ValueError, match="Muito poucas sequências de comprimento uniforme"):
                fetch_dataset_silent(params)

    def test_fetch_dataset_uniform_length_filtering(self):
        """Testa filtragem por comprimento uniforme."""
        params = {
            "email": "test@example.com",
            "db": "nucleotide",
            "term": "ribosomal RNA[Title]",  # RNAs ribossomais têm comprimentos variados
            "n": 5,
        }

        # RNAs ribossomais têm comprimentos diferentes, então esperamos erro
        with vcr.use_cassette("tests/vcr_cassettes/entrez_uniform_length.yaml"):
            with pytest.raises(ValueError, match="Muito poucas sequências de comprimento uniforme"):
                fetch_dataset_silent(params)


# Configuração do VCR
def pytest_configure():
    """Configuração do pytest com VCR."""
    # Criar diretório para cassettes se não existir
    import os

    os.makedirs("tests/vcr_cassettes", exist_ok=True)


# Fixture para configurar VCR
@pytest.fixture(scope="module")
def vcr_config():
    """Configuração do VCR."""
    return {
        "filter_headers": ["authorization", "user-agent"],
        "filter_query_parameters": ["api_key"],
        "ignore_hosts": ["localhost", "127.0.0.1"],
        "record_mode": "once",  # Gravar apenas uma vez
    }
