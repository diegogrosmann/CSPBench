"""
Comprehensive tests for src/domain/config.py to improve code coverage.
These tests focus on areas not covered by existing tests.
"""

import os
import tempfile
from dataclasses import asdict
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from src.domain.config import (
    AlgParamsConfig,
    CSPBenchConfig,
    ExperimentTaskConfig,
    FileDatasetConfig,
    MetadataConfig,
    OptimizationTaskConfig,
    OutputConfig,
    ResultsConfig,
    ResultsFormats,
    SyntheticDatasetConfig,
    TaskConfig,
    TasksGroupConfig,
)
from src.domain.errors import ConfigurationError


class TestMetadataConfig:
    """Test MetadataConfig functionality."""

    def test_metadata_creation_with_all_fields(self):
        """Test metadata creation with all fields."""
        metadata = MetadataConfig(
            name="Test Config",
            description="A test configuration",
            author="Test Author",
            version="1.0.0",
            creation_date="2024-01-01",
            tags=["test", "config"],
        )
        assert metadata.name == "Test Config"
        assert metadata.description == "A test configuration"
        assert metadata.author == "Test Author"
        assert metadata.version == "1.0.0"
        assert metadata.creation_date == "2024-01-01"
        assert metadata.tags == ["test", "config"]

    def test_metadata_creation_minimal_fields(self):
        """Test metadata creation with minimal required fields."""
        metadata = MetadataConfig(
            name="Minimal Config",
            description="Minimal configuration",
            author="Min Author",
            version="1.0.0",
            creation_date="2024-01-01",
        )
        assert metadata.name == "Minimal Config"
        assert metadata.author == "Min Author"
        assert metadata.version == "1.0.0"
        assert metadata.creation_date == "2024-01-01"

    def test_metadata_serialization(self):
        """Test metadata serialization."""
        metadata = MetadataConfig(
            name="Serialize Config",
            description="For serialization testing",
            author="Serialize Author",
            version="2.0.0",
            creation_date="2024-02-01",
        )
        result = asdict(metadata)
        assert result["name"] == "Serialize Config"
        assert result["author"] == "Serialize Author"
        assert result["version"] == "2.0.0"


class TestOutputConfig:
    """Test OutputConfig functionality."""

    def test_output_config_creation(self):
        """Test output config creation."""
        formats = ResultsFormats(csv=True, json=False, parquet=True, pickle=False)
        results_config = ResultsConfig(formats=formats, partial_results=True)
        output_config = OutputConfig(logging=True, results=results_config)
        assert output_config.logging is True
        assert output_config.results.partial_results is True
        assert output_config.results.formats.csv is True
        assert output_config.results.formats.json is False

    def test_output_config_serialization(self):
        """Test output config serialization."""
        formats = ResultsFormats(csv=False, json=True, parquet=False, pickle=True)
        results_config = ResultsConfig(formats=formats, partial_results=False)
        output_config = OutputConfig(logging=True, results=results_config)
        result = asdict(output_config)
        assert result["logging"] is True
        assert result["results"]["formats"]["json"] is True


class TestTaskConfig:
    """Test TaskConfig functionality."""

    def test_experiment_task_config_creation(self):
        """Test experiment task config creation."""
        task_config = ExperimentTaskConfig(
            id="exp1",
            name="Experiment 1",
            datasets=["dataset1", "dataset2"],
            algorithms=["algo1", "algo2"],
            repetitions=10,
        )
        assert task_config.type == "experiment"
        assert task_config.datasets == ["dataset1", "dataset2"]
        assert task_config.algorithms == ["algo1", "algo2"]
        assert task_config.repetitions == 10

    def test_optimization_task_config_creation(self):
        """Test optimization task config creation."""
        task_config = OptimizationTaskConfig(
            id="opt1",
            name="Optimization 1",
            datasets=["opt_dataset"],
            algorithms=["genetic"],
            parameters={"population_size": 50, "generations": 100},
        )
        assert task_config.type == "optimization"
        assert task_config.datasets == ["opt_dataset"]
        assert task_config.algorithms == ["genetic"]
        assert task_config.parameters["population_size"] == 50
        assert task_config.parameters["generations"] == 100

    def test_task_config_serialization(self):
        """Test task config serialization."""
        task_config = ExperimentTaskConfig(
            id="bench1",
            name="Benchmark Task",
            datasets=["bench1"],
            algorithms=["baseline"],
            repetitions=5,
        )
        result = asdict(task_config)
        assert result["type"] == "experiment"
        assert result["datasets"] == ["bench1"]
        assert result["algorithms"] == ["baseline"]
        assert result["repetitions"] == 5


class TestDatasetConfig:
    """Test DatasetConfig functionality."""

    def test_synthetic_dataset_config_creation(self):
        """Test synthetic dataset configuration."""
        dataset_config = SyntheticDatasetConfig(
            id="synthetic_test",
            name="Synthetic Test Dataset",
            mode="random",
            n=100,
            L=50,
            alphabet="ACGT",
        )
        assert dataset_config.id == "synthetic_test"
        assert dataset_config.name == "Synthetic Test Dataset"
        assert dataset_config.type == "synthetic"
        assert dataset_config.mode == "random"
        assert dataset_config.n == 100
        assert dataset_config.L == 50
        assert dataset_config.alphabet == "ACGT"

    def test_file_dataset_config_creation(self):
        """Test file-based dataset configuration."""
        dataset_config = FileDatasetConfig(
            id="file_test", name="File Test Dataset", filename="/path/to/dataset.fasta"
        )
        assert dataset_config.id == "file_test"
        assert dataset_config.name == "File Test Dataset"
        assert dataset_config.type == "file"
        assert dataset_config.filename == "/path/to/dataset.fasta"

    def test_dataset_config_with_noise(self):
        """Test dataset config with noise parameters."""
        dataset_config = SyntheticDatasetConfig(
            id="noisy_dataset",
            name="Noisy Dataset",
            mode="noise",
            n=50,
            L=30,
            parameters_mode={"noise_rate": 0.1, "base_sequences": ["ACGT", "TGCA"]},
        )
        assert dataset_config.id == "noisy_dataset"
        assert dataset_config.mode == "noise"
        assert dataset_config.parameters_mode["noise_rate"] == 0.1
        assert dataset_config.parameters_mode["base_sequences"] == ["ACGT", "TGCA"]

    def test_dataset_config_serialization(self):
        """Test dataset config serialization."""
        dataset_config = SyntheticDatasetConfig(
            id="serialize_test",
            name="Serialization Test",
            mode="clustered",
            n=25,
            L=40,
            parameters_mode={"clusters": 3, "cluster_distance": 0.3},
        )
        result = asdict(dataset_config)
        assert result["id"] == "serialize_test"
        assert result["name"] == "Serialization Test"
        assert result["type"] == "synthetic"
        assert result["mode"] == "clustered"
        assert result["n"] == 25
        assert result["L"] == 40
        assert result["parameters_mode"]["clusters"] == 3
        assert result["parameters_mode"]["cluster_distance"] == 0.3


class TestCSPBenchConfig:
    """Test CSPBenchConfig main configuration class."""

    def create_sample_config(self):
        """Create a sample configuration for testing."""
        metadata = MetadataConfig(
            name="Sample Config",
            description="A sample configuration for testing",
            author="Test Suite",
            version="1.0.0",
            creation_date="2024-01-01",
        )

        formats = ResultsFormats(csv=True, json=False, parquet=True, pickle=False)
        results_config = ResultsConfig(formats=formats, partial_results=False)

        output = OutputConfig(logging=True, results=results_config)

        task = ExperimentTaskConfig(
            id="test_task",
            name="Test Task",
            datasets=["test_dataset"],
            algorithms=["test_algorithm"],
            repetitions=5,
        )

        tasks_group = TasksGroupConfig(type="experiment", items=[task])

        datasets = {
            "test_dataset": SyntheticDatasetConfig(
                id="test_dataset", name="Test Dataset", mode="random", n=10, L=20
            )
        }

        algorithms = {
            "test_algorithm": AlgParamsConfig(name="test_algorithm", params={})
        }

        return CSPBenchConfig(
            metadata=metadata,
            output=output,
            tasks=tasks_group,
            datasets=datasets,
            algorithms=algorithms,
        )

    def test_cspbench_config_creation(self):
        """Test CSPBenchConfig creation."""
        config = self.create_sample_config()

        assert config.metadata.name == "Sample Config"
        assert config.output.logging is True
        assert config.tasks.type == "experiment"
        assert len(config.datasets) == 1
        assert len(config.algorithms) == 1

    def test_cspbench_config_serialization(self):
        """Test CSPBenchConfig serialization."""
        config = self.create_sample_config()
        result = asdict(config)

        assert result["metadata"]["name"] == "Sample Config"
        assert result["output"]["logging"] is True
        assert result["tasks"]["type"] == "experiment"
        assert "test_dataset" in result["datasets"]
        assert "test_algorithm" in result["algorithms"]

    def test_cspbench_config_validation(self):
        """Test configuration validation."""
        config = self.create_sample_config()

        # Should not raise any exceptions for valid config
        try:
            # Just accessing properties should work
            _ = config.metadata.name
            _ = config.output.logging
            _ = config.tasks.type
        except Exception as e:
            pytest.fail(f"Valid config should not raise exceptions: {e}")

    def test_cspbench_config_empty_collections(self):
        """Test configuration with empty datasets and algorithms."""
        formats = ResultsFormats(csv=True, json=False, parquet=False, pickle=False)
        results_config = ResultsConfig(formats=formats, partial_results=False)
        output = OutputConfig(logging=True, results=results_config)

        task = ExperimentTaskConfig(
            id="empty", name="Empty Task", datasets=[], algorithms=[]
        )
        tasks_group = TasksGroupConfig(type="experiment", items=[task])

        config = CSPBenchConfig(
            metadata=MetadataConfig(
                name="Empty",
                description="Empty config",
                author="Test",
                version="1.0.0",
                creation_date="2024-01-01",
            ),
            output=output,
            tasks=tasks_group,
            datasets={},
            algorithms={},
        )

        assert len(config.datasets) == 0
        assert len(config.algorithms) == 0
        assert config.tasks.items[0].datasets == []
        assert config.tasks.items[0].algorithms == []
