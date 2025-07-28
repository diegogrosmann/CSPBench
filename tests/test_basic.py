"""
Teste básico para verificar se o ambiente está funcionando.
"""

import pytest
from src.domain import Dataset, SyntheticDatasetGenerator, hamming_distance


def test_basic_domain_functionality():
    """Teste básico do domínio."""
    # Teste dataset
    sequences = ["ACGT", "ATGT", "GCGT"]
    dataset = Dataset(sequences=sequences, metadata={})
    
    assert dataset.size == 3
    assert dataset.length == 4
    
    # Teste hamming distance
    assert hamming_distance("ACGT", "ACGT") == 0
    assert hamming_distance("ACGT", "ATGT") == 1
    
    # Teste synthetic generator
    synthetic = SyntheticDatasetGenerator.generate_random(
        n=3, length=5, alphabet="ACGT", seed=42
    )
    
    assert synthetic.size == 3
    assert synthetic.length == 5


def test_config_parser():
    """Teste básico do parser de configuração."""
    from src.application.services.config_parser import ConfigurationParser
    
    parser = ConfigurationParser()
    assert parser is not None


def test_experiment_service():
    """Teste básico do serviço de experimentos."""
    from src.application.services.experiment_service import ExperimentService
    from tests.unit.application.fakes import (
        FakeDatasetRepository, 
        FakeExportPort, 
        FakeExecutorPort,
        FakeAlgorithmRegistry
    )
    
    dataset_repo = FakeDatasetRepository()
    exporter = FakeExportPort()
    executor = FakeExecutorPort()
    algo_registry = FakeAlgorithmRegistry()
    
    service = ExperimentService(
        dataset_repo=dataset_repo,
        exporter=exporter,
        executor=executor,
        algo_registry=algo_registry
    )
    
    assert service is not None
