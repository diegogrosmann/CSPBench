"""
Simple tests for src/domain/config.py to improve code coverage.
Focus on testing the actual structure of the classes as they are.
"""

import pytest
from src.domain.config import (
    MetadataConfig,
    AlgParamsConfig,
    AlgorithmsPresetConfig,
    SyntheticDatasetConfig,
    FileDatasetConfig,
    TaskConfig,
    CSPBenchConfig,
    ExperimentTasksConfig
)


class TestMetadataConfig:
    """Test MetadataConfig functionality."""
    
    def test_metadata_creation_complete(self):
        """Test metadata creation with all required fields."""
        metadata = MetadataConfig(
            name="Test Config",
            description="A test configuration", 
            author="Test Author",
            version="1.0.0",
            creation_date="2025-01-01"
        )
        assert metadata.name == "Test Config"
        assert metadata.description == "A test configuration"
        assert metadata.author == "Test Author"
        assert metadata.version == "1.0.0"
        assert metadata.creation_date == "2025-01-01"
        
    def test_metadata_with_tags(self):
        """Test metadata creation with tags."""
        metadata = MetadataConfig(
            name="Tagged Config",
            description="Configuration with tags",
            author="Tag Author", 
            version="2.0.0",
            creation_date="2025-01-02",
            tags=["test", "benchmark", "experimental"]
        )
        assert metadata.tags == ["test", "benchmark", "experimental"]
        
    def test_metadata_default_tags(self):
        """Test metadata creation with default empty tags."""
        metadata = MetadataConfig(
            name="Default Tags",
            description="Configuration with default tags",
            author="Default Author",
            version="1.0.0", 
            creation_date="2025-01-01"
        )
        assert metadata.tags == []


class TestAlgParamsConfig:
    """Test AlgParamsConfig functionality."""
    
    def test_alg_params_creation(self):
        """Test algorithm parameters creation."""
        alg_params = AlgParamsConfig(
            name="test_algorithm",
            params={"param1": 10, "param2": "test"}
        )
        assert alg_params.name == "test_algorithm"
        assert alg_params.params["param1"] == 10
        assert alg_params.params["param2"] == "test"
        
    def test_alg_params_empty_params(self):
        """Test algorithm parameters with empty parameters."""
        alg_params = AlgParamsConfig(name="simple_algorithm")
        assert alg_params.name == "simple_algorithm"
        assert alg_params.params == {}
        
    def test_alg_params_complex_params(self):
        """Test algorithm parameters with complex parameters."""
        complex_params = {
            "population_size": 100,
            "mutation_rate": 0.01,
            "selection": "tournament", 
            "nested": {"subparam": "value", "numbers": [1, 2, 3]}
        }
        alg_params = AlgParamsConfig(name="genetic_algorithm", params=complex_params)
        assert alg_params.params["population_size"] == 100
        assert alg_params.params["nested"]["subparam"] == "value"


class TestAlgorithmsPresetConfig:
    """Test AlgorithmsPresetConfig functionality."""
    
    def test_algorithms_preset_creation(self):
        """Test algorithms preset creation."""
        preset = AlgorithmsPresetConfig(
            id="test_preset",
            name="Test Preset",
            description="A test preset"
        )
        assert preset.id == "test_preset"
        assert preset.name == "Test Preset"
        assert preset.description == "A test preset"
        assert preset.items == []
        
    def test_algorithms_preset_with_items(self):
        """Test algorithms preset with algorithm items."""
        alg1 = AlgParamsConfig(name="algo1", params={"param1": 1})
        alg2 = AlgParamsConfig(name="algo2", params={"param2": 2})
        
        preset = AlgorithmsPresetConfig(
            id="multi_preset", 
            name="Multi Algorithm Preset",
            description="Multiple algorithms",
            items=[alg1, alg2]
        )
        assert len(preset.items) == 2
        assert preset.items[0].name == "algo1"
        assert preset.items[1].name == "algo2"


class TestDatasetConfigs:
    """Test dataset configuration classes."""
    
    def test_synthetic_dataset_creation(self):
        """Test synthetic dataset configuration."""
        dataset = SyntheticDatasetConfig(
            id="syn_test",
            name="Synthetic Test Dataset"
        )
        assert dataset.id == "syn_test"
        assert dataset.name == "Synthetic Test Dataset"
        assert dataset.type == "synthetic"
        assert dataset.mode == "random"  # default value
        
    def test_synthetic_dataset_with_params(self):
        """Test synthetic dataset with specific parameters."""
        dataset = SyntheticDatasetConfig(
            id="custom_syn",
            name="Custom Synthetic Dataset",
            mode="clustered",
            n=100,
            L=50,
            alphabet="ACGT"
        )
        assert dataset.mode == "clustered"
        assert dataset.n == 100
        assert dataset.L == 50
        assert dataset.alphabet == "ACGT"
        
    def test_file_dataset_creation(self):
        """Test file dataset configuration."""
        dataset = FileDatasetConfig(
            id="file_test",
            name="File Test Dataset",
            filename="data.fasta"
        )
        assert dataset.id == "file_test"
        assert dataset.name == "File Test Dataset"
        assert dataset.filename == "data.fasta"
        assert dataset.type == "file"


class TestTaskConfig:
    """Test TaskConfig functionality."""
    
    def test_task_creation_basic(self):
        """Test basic task configuration creation."""
        task = TaskConfig(
            id="test_task",
            name="Test Task",
            type="experiment"
        )
        assert task.id == "test_task"
        assert task.name == "Test Task"
        assert task.type == "experiment"
        assert task.datasets == []
        assert task.algorithms == []
        
    def test_task_with_datasets_algorithms(self):
        """Test task configuration with datasets and algorithms."""
        task = TaskConfig(
            id="full_task",
            name="Full Task",
            type="optimization",
            datasets=["dataset1", "dataset2"],
            algorithms=["algo1", "algo2"]
        )
        assert task.datasets == ["dataset1", "dataset2"]
        assert task.algorithms == ["algo1", "algo2"]


class TestCSPBenchConfig:
    """Test CSPBenchConfig main configuration class."""
    
    def test_cspbench_config_minimal(self):
        """Test minimal CSPBenchConfig creation."""
        metadata = MetadataConfig(
            name="Minimal Config",
            description="Minimal test config",
            author="Test Author",
            version="1.0.0",
            creation_date="2025-01-01"
        )
        
        tasks = ExperimentTasksConfig(type="experiment", items=[])
        
        config = CSPBenchConfig(
            metadata=metadata,
            datasets={},
            algorithms={},
            tasks=tasks
        )
        assert config.metadata.name == "Minimal Config"
        assert config.output is None
        assert config.datasets == {}
        assert config.algorithms == {}
        
    def test_cspbench_config_to_dict(self):
        """Test CSPBenchConfig to_dict method."""
        metadata = MetadataConfig(
            name="Dict Config",
            description="Test to_dict",
            author="Dict Author",
            version="1.0.0",
            creation_date="2025-01-01"
        )
        
        tasks = ExperimentTasksConfig(type="experiment", items=[])
        
        config = CSPBenchConfig(
            metadata=metadata,
            datasets={},
            algorithms={},
            tasks=tasks
        )
        
        result = config.to_dict()
        assert result["metadata"]["name"] == "Dict Config"
        assert result["datasets"] == {}
        assert result["algorithms"] == {}


class TestAdditionalConfigs:
    """Test additional configuration classes for better coverage."""
    
    def test_experiment_tasks_config(self):
        """Test ExperimentTasksConfig creation."""
        task = TaskConfig(
            id="exp_task",
            name="Experiment Task",
            type="experiment"
        )
        
        exp_tasks = ExperimentTasksConfig(
            type="experiment",
            items=[task]
        )
        assert exp_tasks.type == "experiment"
        assert len(exp_tasks.items) == 1
        assert exp_tasks.items[0].id == "exp_task"
        
    def test_synthetic_dataset_modes(self):
        """Test different modes of synthetic datasets."""
        # Test random mode (default)
        random_dataset = SyntheticDatasetConfig(
            id="random_test",
            name="Random Dataset"
        )
        assert random_dataset.mode == "random"
        
        # Test clustered mode
        clustered_dataset = SyntheticDatasetConfig(
            id="clustered_test", 
            name="Clustered Dataset",
            mode="clustered"
        )
        assert clustered_dataset.mode == "clustered"
        
        # Test noise mode
        noise_dataset = SyntheticDatasetConfig(
            id="noise_test",
            name="Noise Dataset", 
            mode="noise"
        )
        assert noise_dataset.mode == "noise"
        
        # Test mutations mode
        mutations_dataset = SyntheticDatasetConfig(
            id="mutations_test",
            name="Mutations Dataset",
            mode="mutations"
        )
        assert mutations_dataset.mode == "mutations"