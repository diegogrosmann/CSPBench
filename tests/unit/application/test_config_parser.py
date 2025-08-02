"""
Testes unitários para ConfigurationParser.

Testa o parsing e validação de configurações de batch.
"""

import json
import tempfile
from pathlib import Path

import pytest
import yaml

from src.application.services.config_parser import (
    BatchMetadata,
    ConfigurationParser,
    ConfigurationValidator,
    ExportConfig,
    InfrastructureConfig,
    LoggingConfig,
    MonitoringConfig,
    OptimizationConfig,
    PlotsConfig,
    SensitivityConfig,
    SystemConfig,
)
from src.domain.errors import (
    BatchConfigurationError,
    OptimizationConfigurationError,
    SensitivityConfigurationError,
)


class TestConfigurationParser:
    """Tests for ConfigurationParser."""

    @pytest.fixture
    def valid_batch_config(self):
        """Valid batch configuration for testing."""
        return {
            "metadados": {
                "nome": "Test Batch",
                "descricao": "Test description",
                "autor": "Test Author",
                "versao": "1.0",
                "data_criacao": "2025-01-28",
                "tags": ["test", "example"],
                "timeout_global": 3600,
            },
            "infrastructure": {
                "history": {
                    "save_history": True,
                    "plot_history": True,
                    "history_frequency": 1,
                },
                "result": {
                    "save_partial_results": True,
                    "partial_file": "partial.json",
                },
            },
            "datasets": [
                {
                    "id": "test_dataset",
                    "nome": "Test Dataset",
                    "tipo": "synthetic",
                    "parametros": {
                        "n": 10,
                        "L": 20,
                        "alphabet": "ACGT",
                        "noise": 0.1,
                    },
                }
            ],
            "algorithms": [
                {
                    "id": "test_config",
                    "nome": "Test Config",
                    "algorithms": ["Baseline"],
                    "algorithm_params": {"Baseline": {"tie_break": "lex"}},
                }
            ],
            "task": {"type": "execution"},
            "execution": {
                "executions": [
                    {
                        "nome": "Test Execution",
                        "datasets": ["test_dataset"],
                        "algorithms": ["test_config"],
                        "repetitions": 1,
                    }
                ]
            },
            "export": {
                "enabled": True,
                "destination": "outputs/test",
                "formats": {"json": True, "csv": False},
            },
            "plots": {
                "enabled": True,
                "plot_convergence": True,
                "style": "seaborn-v0_8",
            },
            "monitoring": {
                "enabled": True,
                "interface": "simple",
                "update_interval": 5,
            },
            "logging": {
                "level": "INFO",
                "output": {"console": True, "file": True},
            },
            "system": {
                "reproducibility": {"global_seed": 42, "strict_mode": True},
                "checkpointing": {"enabled": True, "interval": 5},
            },
            "resources": {
                "parallel": {
                    "enabled": True,
                    "max_workers": 4,
                    "internal_jobs": 2,
                }
            },
        }

    def test_load_file_yaml(self, valid_batch_config):
        """Test loading YAML file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            loaded_config = ConfigurationParser.load_file(temp_path)
            assert loaded_config == valid_batch_config
        finally:
            Path(temp_path).unlink()

    def test_load_file_json(self, valid_batch_config):
        """Test loading JSON file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            loaded_config = ConfigurationParser.load_file(temp_path)
            assert loaded_config == valid_batch_config
        finally:
            Path(temp_path).unlink()

    def test_load_file_not_found(self):
        """Test loading non-existent file."""
        with pytest.raises(BatchConfigurationError, match="File not found"):
            ConfigurationParser.load_file("nonexistent.yaml")

    def test_load_file_unsupported_format(self):
        """Test loading unsupported file format."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("test content")
            temp_path = f.name

        try:
            with pytest.raises(BatchConfigurationError, match="Unsupported format"):
                ConfigurationParser.load_file(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_parse_metadata_success(self, valid_batch_config):
        """Test successful metadata parsing."""
        metadata = ConfigurationParser.parse_metadata(valid_batch_config)

        assert isinstance(metadata, BatchMetadata)
        assert metadata.nome == "Test Batch"
        assert metadata.descricao == "Test description"
        assert metadata.autor == "Test Author"
        assert metadata.versao == "1.0"
        assert metadata.timeout_global == 3600
        assert "test" in metadata.tags

    def test_parse_metadata_legacy_format(self):
        """Test parsing metadata with legacy batch_info format."""
        config = {
            "batch_info": {
                "nome": "Legacy Batch",
                "descricao": "Legacy description",
                "autor": "Legacy Author",
                "versao": "0.1",
                "data_criacao": "2024-12-01",
                "tags": ["legacy"],
                "timeout_global": 1800,
            }
        }

        metadata = ConfigurationParser.parse_metadata(config)
        assert metadata.nome == "Legacy Batch"
        assert metadata.timeout_global == 1800

    def test_parse_metadata_missing(self):
        """Test parsing metadata when section is missing."""
        config = {"other_section": {}}

        with pytest.raises(
            BatchConfigurationError, match="metadados/batch_info section not found"
        ):
            ConfigurationParser.parse_metadata(config)

    def test_parse_infrastructure_config(self, valid_batch_config):
        """Test parsing infrastructure configuration."""
        config = ConfigurationParser.parse_infrastructure_config(valid_batch_config)

        assert isinstance(config, InfrastructureConfig)
        assert config.history is not None
        assert config.history["save_history"] is True
        assert config.result is not None
        assert config.result["save_partial_results"] is True

    def test_parse_export_config(self, valid_batch_config):
        """Test parsing export configuration."""
        config = ConfigurationParser.parse_export_config(valid_batch_config)

        assert isinstance(config, ExportConfig)
        assert config.enabled is True
        assert config.destination == "outputs/test"
        assert config.formats["json"] is True

    def test_parse_plots_config(self, valid_batch_config):
        """Test parsing plots configuration."""
        config = ConfigurationParser.parse_plots_config(valid_batch_config)

        assert isinstance(config, PlotsConfig)
        assert config.enabled is True
        assert config.plot_convergence is True
        assert config.style == "seaborn-v0_8"

    def test_parse_monitoring_config(self, valid_batch_config):
        """Test parsing monitoring configuration."""
        config = ConfigurationParser.parse_monitoring_config(valid_batch_config)

        assert isinstance(config, MonitoringConfig)
        assert config.enabled is True
        assert config.interface == "simple"
        assert config.update_interval == 5

    def test_parse_logging_config(self, valid_batch_config):
        """Test parsing logging configuration."""
        config = ConfigurationParser.parse_logging_config(valid_batch_config)

        assert isinstance(config, LoggingConfig)
        assert config.level == "INFO"
        assert config.output["console"] is True

    def test_parse_system_config(self, valid_batch_config):
        """Test parsing system configuration."""
        config = ConfigurationParser.parse_system_config(valid_batch_config)

        assert isinstance(config, SystemConfig)
        assert config.reproducibility["global_seed"] == 42
        assert config.checkpointing["enabled"] is True

    def test_parse_resources_config(self, valid_batch_config):
        """Test parsing resources configuration."""
        config = ConfigurationParser.parse_resources_config(valid_batch_config)

        assert config["enabled"] is True
        assert config["max_workers"] == 4
        assert config["internal_jobs"] == 2

    def test_parse_default_configs(self):
        """Test parsing with missing optional sections returns defaults."""
        minimal_config = {
            "metadados": {
                "nome": "Minimal",
                "descricao": "Minimal config",
                "autor": "Test",
                "versao": "1.0",
                "data_criacao": "2025-01-28",
                "tags": [],
                "timeout_global": 3600,
            }
        }

        # Test that missing sections return default configurations
        infrastructure = ConfigurationParser.parse_infrastructure_config(minimal_config)
        export_config = ConfigurationParser.parse_export_config(minimal_config)
        plots_config = ConfigurationParser.parse_plots_config(minimal_config)
        monitoring_config = ConfigurationParser.parse_monitoring_config(minimal_config)
        logging_config = ConfigurationParser.parse_logging_config(minimal_config)
        system_config = ConfigurationParser.parse_system_config(minimal_config)

        # Check defaults
        assert isinstance(infrastructure, InfrastructureConfig)
        assert export_config.enabled is True
        assert plots_config.enabled is True
        assert monitoring_config.enabled is True
        assert logging_config.level == "INFO"
        assert isinstance(system_config, SystemConfig)

    def test_parse_optimization_configs_new_structure(self):
        """Test parsing optimization configurations with new structure."""
        config = {
            "optimization": {
                "method": "optuna",
                "optimizations": [
                    {
                        "nome": "Test Optimization",
                        "study_name": "test_study",
                        "direction": "minimize",
                        "n_trials": 50,
                        "timeout_per_trial": 300,
                        "target_datasets": ["dataset1"],
                        "target_algorithm": "BLF-GA",
                        "parameters": {
                            "BLF-GA": {
                                "pop_size": {"type": "int", "low": 50, "high": 200}
                            }
                        },
                    }
                ],
            }
        }

        configs = ConfigurationParser.parse_optimization_configs(config)
        assert len(configs) == 1
        assert isinstance(configs[0], OptimizationConfig)
        assert configs[0].nome == "Test Optimization"
        assert configs[0].n_trials == 50

    def test_parse_sensitivity_configs_new_structure(self):
        """Test parsing sensitivity configurations with new structure."""
        config = {
            "sensitivity": {
                "method": "SALib",
                "global_salib_config": {"n_samples": 1000, "repetitions_per_sample": 3},
                "analyses": [
                    {
                        "nome": "Test Analysis",
                        "analysis_method": "morris",
                        "target_datasets": ["dataset1"],
                        "target_algorithm": "BLF-GA",
                        "parameters": {
                            "pop_size": {
                                "type": "integer",
                                "bounds": [50, 200],
                                "default": 100,
                            }
                        },
                        "output_metrics": ["distance", "execution_time"],
                    }
                ],
            }
        }

        configs = ConfigurationParser.parse_sensitivity_configs(config)
        assert len(configs) == 1
        assert isinstance(configs[0], SensitivityConfig)
        assert configs[0].nome == "Test Analysis"
        assert configs[0].analysis_method == "morris"


class TestConfigurationValidator:
    """Tests for ConfigurationValidator."""

    def test_validate_batch_structure_execution(self):
        """Test validation of execution batch structure."""
        config = {
            "metadados": {"nome": "test"},
            "datasets": [],
            "algorithms": [],
            "task": {"type": "execution"},
            "execution": {},
        }

        task_type = ConfigurationValidator.validate_batch_structure(config)
        assert task_type == "execution"

    def test_validate_batch_structure_optimization(self):
        """Test validation of optimization batch structure."""
        config = {
            "metadados": {"nome": "test"},
            "datasets": [],
            "algorithms": [],
            "task": {"type": "optimization"},
            "optimization": {},
        }

        task_type = ConfigurationValidator.validate_batch_structure(config)
        assert task_type == "optimization"

    def test_validate_batch_structure_sensitivity(self):
        """Test validation of sensitivity batch structure."""
        config = {
            "metadados": {"nome": "test"},
            "datasets": [],
            "algorithms": [],
            "task": {"type": "sensitivity"},
            "sensitivity": {},
        }

        task_type = ConfigurationValidator.validate_batch_structure(config)
        assert task_type == "sensitivity"

    def test_validate_batch_structure_legacy(self):
        """Test validation of legacy batch structure."""
        config = {"experiments": [{"algorithm": "baseline", "dataset": "test"}]}

        task_type = ConfigurationValidator.validate_batch_structure(config)
        assert task_type == "execution"

    def test_validate_batch_structure_missing_sections(self):
        """Test validation with missing required sections."""
        config = {
            "metadados": {"nome": "test"},
            "task": {"type": "execution"},
            # Missing datasets and algorithms
        }

        with pytest.raises(BatchConfigurationError, match="Required sections missing"):
            ConfigurationValidator.validate_batch_structure(config)

    def test_validate_batch_structure_backward_compatibility(self):
        """Test validation with batch_info instead of metadados."""
        config = {
            "batch_info": {"nome": "test"},  # legacy name
            "datasets": [],
            "algorithms": [],
            "task": {"type": "execution"},
        }

        task_type = ConfigurationValidator.validate_batch_structure(config)
        assert task_type == "execution"

    def test_validate_batch_structure_unrecognized(self):
        """Test validation with unrecognized structure."""
        config = {"unknown_section": {}}

        with pytest.raises(
            BatchConfigurationError, match="Configuration structure not recognized"
        ):
            ConfigurationValidator.validate_batch_structure(config)

    def test_validate_optimization_config(self):
        """Test optimization configuration validation."""
        config = OptimizationConfig(
            nome="Test",
            study_name="test_study",
            direction="minimize",
            n_trials=50,
            timeout_per_trial=300,
            target_datasets=["dataset1"],
            target_algorithm="BLF-GA",
            parameters={"param1": {"type": "int", "low": 1, "high": 10}},
        )

        # Should not raise exception
        ConfigurationValidator.validate_optimization_config(config)

    def test_validate_optimization_config_empty_datasets(self):
        """Test optimization configuration validation with empty datasets."""
        config = OptimizationConfig(
            nome="Test",
            study_name="test_study",
            direction="minimize",
            n_trials=50,
            timeout_per_trial=300,
            target_datasets=[],  # Empty
            target_algorithm="BLF-GA",
            parameters={"param1": {"type": "int"}},
        )

        with pytest.raises(
            OptimizationConfigurationError, match="target_datasets list cannot be empty"
        ):
            ConfigurationValidator.validate_optimization_config(config)

    def test_validate_sensitivity_config(self):
        """Test sensitivity configuration validation."""
        config = SensitivityConfig(
            nome="Test",
            analysis_method="morris",
            target_datasets=["dataset1"],
            target_algorithm="BLF-GA",
            n_samples=1000,
            repetitions_per_sample=3,
            parameters={"param1": {"type": "integer", "bounds": [1, 10]}},
            output_metrics=["distance"],
        )

        # Should not raise exception
        ConfigurationValidator.validate_sensitivity_config(config)

    def test_validate_sensitivity_config_invalid_method(self):
        """Test sensitivity configuration validation with invalid method."""
        config = SensitivityConfig(
            nome="Test",
            analysis_method="invalid_method",
            target_datasets=["dataset1"],
            target_algorithm="BLF-GA",
            n_samples=1000,
            repetitions_per_sample=3,
            parameters={"param1": {"type": "integer"}},
            output_metrics=["distance"],
        )

        with pytest.raises(
            SensitivityConfigurationError, match="Invalid analysis method"
        ):
            ConfigurationValidator.validate_sensitivity_config(config)
