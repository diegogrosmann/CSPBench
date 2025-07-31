"""
Configuração de Testes para CSPBench

Este módulo configura pytest e fornece fixtures comuns para todos os testes.
"""

import shutil
import tempfile
from pathlib import Path
from typing import Any, Dict, Generator, List

import pytest


@pytest.fixture
def sample_sequences() -> List[str]:
    """Fixture com sequências de exemplo para testes."""
    return ["ACGTACGT", "AGGTACGT", "ACGTACCT", "ACTTACGT"]


@pytest.fixture
def small_sequences() -> List[str]:
    """Fixture com sequências pequenas para testes rápidos."""
    return ["ACGT", "AGCT", "ATCT"]


@pytest.fixture
def large_sequences() -> List[str]:
    """Fixture com sequências maiores para testes de performance."""
    return [
        "ACGTACGTACGTACGTACGTACGT",
        "AGGTACGTACGTACGTACGTACGT",
        "ACGTACCTACGTACGTACGTACGT",
        "ACTTACGTACGTACGTACGTACGT",
        "ACGTACGTACCTACGTACGTACGT",
    ]


@pytest.fixture
def protein_sequences() -> List[str]:
    """Fixture com sequências de proteínas para testes."""
    return [
        "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHG",
        "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKGHPETLEKFDRVKHLKTEAEMKASEDLKKHG",
        "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHD",
    ]


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Fixture que cria um diretório temporário para testes."""
    temp_dir = Path(tempfile.mkdtemp())
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def sample_fasta_content() -> str:
    """Fixture com conteúdo FASTA para testes."""
    return """>seq_1 Test sequence 1
ACGTACGT
>seq_2 Test sequence 2
AGGTACGT
>seq_3 Test sequence 3
ACGTACCT
"""


@pytest.fixture
def sample_txt_content() -> str:
    """Fixture com conteúdo texto para testes."""
    return """ACGTACGT
AGGTACGT
ACGTACCT
ACTTACGT
"""


@pytest.fixture
def synthetic_params() -> Dict[str, Any]:
    """Fixture com parâmetros para geração sintética."""
    return {"n": 10, "L": 20, "alphabet": "ACGT", "noise": 0.15, "seed": 42}


@pytest.fixture
def algorithm_params() -> Dict[str, Any]:
    """Fixture com parâmetros de algoritmos."""
    return {
        "BLF-GA": {
            "population_size": 50,
            "max_generations": 100,
            "mutation_rate": 0.1,
            "crossover_rate": 0.8,
        },
        "CSC": {"n_clusters": 3, "max_iterations": 50},
        "H3-CSP": {"max_iterations": 100, "early_stopping": True},
    }


@pytest.fixture
def mock_entrez_params() -> Dict[str, Any]:
    """Fixture com parâmetros mock para Entrez."""
    return {
        "email": "test@example.com",
        "db": "nucleotide",
        "term": "test query",
        "n": 5,
        "api_key": None,
    }


# Configuração do pytest
def pytest_configure(config):
    """Configuração personalizada do pytest."""
    # Adicionar marcadores customizados
    config.addinivalue_line("markers", "slow: marca testes que demoram para executar")
    config.addinivalue_line("markers", "integration: marca testes de integração")
    config.addinivalue_line("markers", "unit: marca testes unitários")
    config.addinivalue_line("markers", "network: marca testes que requerem internet")


def pytest_collection_modifyitems(config, items):
    """Modifica itens de teste coletados."""
    # Adicionar marker 'slow' para testes que demoram
    for item in items:
        if "integration" in item.nodeid:
            item.add_marker(pytest.mark.integration)
        if "unit" in item.nodeid:
            item.add_marker(pytest.mark.unit)
