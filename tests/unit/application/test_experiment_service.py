"""
Testes unitários para ExperimentService.

Testa a lógica de negócio do serviço usando fakes para as portas.
"""

import tempfile
from pathlib import Path

import pytest
import yaml

from src.application.services import ExperimentService
from src.domain import (
    AlgorithmNotFoundError,
    BatchConfigurationError,
    DatasetNotFoundError,
    OptimizationConfigurationError,
    SensitivityConfigurationError,
)

from .fakes import (
    FakeAlgorithmRegistry,
    FakeDatasetRepository,
    FakeExecutorPort,
    FakeExportPort,
    FakeMonitoringService,
)


class TestExperimentService:
    """Testes para o ExperimentService."""

    @pytest.fixture
    def service(self):
        """Fixture que cria uma instância do serviço com fakes."""
        dataset_repo = FakeDatasetRepository()
        algo_registry = FakeAlgorithmRegistry()
        exporter = FakeExportPort()
        executor = FakeExecutorPort()
        monitoring_service = FakeMonitoringService()

        return ExperimentService(
            dataset_repo=dataset_repo,
            exporter=exporter,
            executor=executor,
            algo_registry=algo_registry,
            monitoring_service=monitoring_service,
        )

    @pytest.fixture
    def valid_batch_config(self):
        """Fixture com configuração válida de batch (novo formato)."""
        return {
            "metadados": {
                "nome": "Test Batch",
                "descricao": "Test description",
                "autor": "Test Author",
                "versao": "1.0",
                "data_criacao": "2025-01-28",
                "tags": ["test"],
                "timeout_global": 3600,
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
                        "seed": 42,
                    },
                }
            ],
            "algorithms": [
                {
                    "id": "test_config",
                    "nome": "Test Config",
                    "algorithms": ["baseline"],
                    "algorithm_params": {"baseline": {"tie_break": "lex"}},
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
                "formats": {"json": True},
                "include": ["summary", "detailed_results"],
            },
            "plots": {
                "enabled": True,
                "plot_convergence": True,
            },
            "monitoring": {
                "enabled": True,
                "interface": "simple",
                "update_interval": 3,
            },
            "logging": {
                "level": "INFO",
            },
            "system": {
                "reproducibility": {"global_seed": 42},
            },
            "resources": {
                "parallel": {
                    "enabled": True,
                    "max_workers": 2,
                    "internal_jobs": 2,
                }
            },
        }

    @pytest.fixture
    def legacy_batch_config(self):
        """Fixture com configuração legada de batch."""
        return {
            "experiments": [
                {
                    "algorithm": "baseline",
                    "dataset": "test_dataset",
                    "params": {"tie_break": "lex"},
                }
            ]
        }

    @pytest.fixture
    def valid_optimization_config(self):
        """Fixture com configuração válida de otimização."""
        return {
            "algorithm": "baseline",
            "dataset": "test_dataset",
            "optimization": {
                "method": "optuna",
                "trials": 50,
                "parameters": {"param1": {"type": "int", "low": 1, "high": 100}},
            },
            "export": {"enabled": True, "destination": "test_optimization"},
        }

    @pytest.fixture
    def valid_sensitivity_config(self):
        """Fixture com configuração válida de sensibilidade."""
        return {
            "algorithm": "baseline",
            "dataset": "test_dataset",
            "parameters": {
                "param1": {"values": [1, 5, 10, 20, 50]},
                "param2": {"values": [0.1, 0.3, 0.5, 0.7, 0.9]},
            },
            "export": {
                "enabled": True,
                "format": "csv",
                "destination": "test_sensitivity",
            },
        }

    def test_run_single_experiment_success(self, service):
        """Testa execução bem-sucedida de experimento único."""
        result = service.run_single_experiment(
            algorithm_name="baseline",
            dataset_id="test_dataset",
            params={"param1": 15},
            timeout=60,
        )

        assert result["algorithm"] == "baseline"
        assert result["status"] == "success"
        assert result["metadata"]["params"]["param1"] == 15
        assert result["metadata"]["timeout"] == 60

    def test_run_single_experiment_algorithm_not_found(self, service):
        """Testa erro quando algoritmo não existe."""
        with pytest.raises(AlgorithmNotFoundError):
            service.run_single_experiment(
                algorithm_name="nonexistent_algo", dataset_id="test_dataset"
            )

    def test_run_single_experiment_dataset_not_found(self, service):
        """Testa erro quando dataset não existe."""
        with pytest.raises(DatasetNotFoundError):
            service.run_single_experiment(
                algorithm_name="baseline", dataset_id="nonexistent_dataset"
            )

    def test_run_batch_success_new_format(self, service, valid_batch_config):
        """Testa execução bem-sucedida de batch no novo formato."""
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            result = service.run_batch(temp_path)

            assert "total_experiments" in result
            assert result["total_experiments"] >= 1
            assert "summary" in result
            assert result["summary"]["successful"] >= 0
            assert result["summary"]["failed"] >= 0
        finally:
            Path(temp_path).unlink()

    def test_run_batch_success_legacy_format(self, service, legacy_batch_config):
        """Testa execução bem-sucedida de batch no formato legado."""
        # Add required sections for legacy format
        legacy_batch_config.update(
            {
                "metadados": {
                    "nome": "Teste Legacy",
                    "descricao": "Teste de formato legado",
                },
                "datasets": ["test_dataset"],
                "algorithms": ["baseline"],
            }
        )

        # Create temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(legacy_batch_config, f)
            temp_path = f.name

        try:
            result = service.run_batch(temp_path)

            assert "total_experiments" in result
            assert result["total_experiments"] >= 1
        finally:
            Path(temp_path).unlink()

    def test_run_batch_with_export(self, service, valid_batch_config):
        """Testa execução de batch com exportação habilitada."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            result = service.run_batch(temp_path)

            # Verify export was called
            exported_files = service._exporter.get_exported_files()
            assert len(exported_files) > 0
            assert any("test" in f["path"] for f in exported_files)
        finally:
            Path(temp_path).unlink()

    def test_run_batch_with_monitoring(self, service, valid_batch_config):
        """Testa execução de batch com monitoramento."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            result = service.run_batch(temp_path)

            # Verify monitoring was used
            monitoring = service._monitoring_service
            assert monitoring.get_interface() == "simple"
            assert monitoring.get_update_interval() == 3
        finally:
            Path(temp_path).unlink()

    def test_run_batch_invalid_config(self, service):
        """Testa erro com configuração inválida de batch."""
        invalid_config = {"invalid": "config"}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(invalid_config, f)
            temp_path = f.name

        try:
            with pytest.raises(BatchConfigurationError):
                service.run_batch(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_run_batch_file_not_found(self, service):
        """Testa erro quando arquivo de batch não existe."""
        with pytest.raises(
            BatchConfigurationError, match="Erro ao parsear configuração"
        ):
            service.run_batch("nonexistent.yaml")

    def test_run_single_experiment_success(self, service):
        """Testa execução bem-sucedida de experimento único."""
        result = service.run_single_experiment(
            algorithm_name="baseline",
            dataset_id="test_dataset",
            params={"param1": 15},
            timeout=60,
        )

        assert isinstance(result["result"], str)
        assert len(result["result"]) > 0
        assert isinstance(result["distance"], int)
        assert result["distance"] >= 0
        assert result["execution_time"] > 0
        assert "metadata" in result

    def test_run_single_experiment_algorithm_not_found(self, service):
        """Testa erro quando algoritmo não existe."""
        with pytest.raises(AlgorithmNotFoundError):
            service.run_single_experiment(
                algorithm_name="nonexistent_algo", dataset_id="test_dataset"
            )

    def test_run_single_experiment_dataset_not_found(self, service):
        """Testa erro quando dataset não existe."""
        with pytest.raises(DatasetNotFoundError):
            service.run_single_experiment(
                algorithm_name="baseline", dataset_id="nonexistent_dataset"
            )

    def test_get_algorithm_info_success(self, service):
        """Testa obtenção de informações de algoritmo."""
        info = service.get_algorithm_info("baseline")

        assert info["name"] == "baseline"
        assert "class" in info

    def test_get_algorithm_info_not_found(self, service):
        """Testa erro ao obter informações de algoritmo inexistente."""
        with pytest.raises(AlgorithmNotFoundError):
            service.get_algorithm_info("nonexistent_algorithm")

    def test_list_datasets(self, service):
        """Testa listagem de datasets disponíveis."""
        datasets = service.list_datasets()

        assert isinstance(datasets, list)
        assert "test_dataset" in datasets
        assert "small_dataset" in datasets

    def test_list_algorithms(self, service):
        """Testa listagem de algoritmos disponíveis."""
        algorithms = service.list_algorithms()

        assert isinstance(algorithms, list)
        assert "baseline" in algorithms
        assert "test_algo" in algorithms

    def test_load_dataset_success(self, service):
        """Testa carregamento bem-sucedido de dataset."""
        dataset = service.load_dataset("test_dataset")

        assert dataset is not None
        assert dataset.size > 0
        assert len(dataset.sequences) > 0

    def test_load_dataset_not_found(self, service):
        """Testa erro ao carregar dataset inexistente."""
        with pytest.raises(DatasetNotFoundError):
            service.load_dataset("nonexistent_dataset")

    def test_configuration_parsing_integration(self, service, valid_batch_config):
        """Testa integração do parsing de configuração."""
        # Test that all configuration sections are parsed correctly
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            # This should not raise any configuration errors
            result = service.run_batch(temp_path)
            assert result is not None
        finally:
            Path(temp_path).unlink()

    def test_export_configuration_applied(self, service, valid_batch_config):
        """Test that export configuration is properly applied."""
        # Modify export config for testing
        valid_batch_config["export"]["formats"]["csv"] = True
        valid_batch_config["export"]["include"] = ["summary"]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            service.run_batch(temp_path)

            # Check that export was called with correct format
            exported_files = service._exporter.get_exported_files()
            assert len(exported_files) > 0
        finally:
            Path(temp_path).unlink()

    def test_plots_configuration_applied(self, service, valid_batch_config):
        """Test that plots configuration is properly applied."""
        valid_batch_config["plots"]["enabled"] = True
        valid_batch_config["plots"]["plot_convergence"] = True
        valid_batch_config["plots"]["style"] = "ggplot"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            # This should process plots configuration without error
            result = service.run_batch(temp_path)
            assert result is not None
        finally:
            Path(temp_path).unlink()

    def test_system_configuration_applied(self, service, valid_batch_config):
        """Test that system configuration is properly applied."""
        # Add system config that exists
        valid_batch_config["system"] = {
            "reproducibility": {"global_seed": 123},
            "performance": {"parallel_execution": True},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            # This should process system configuration without error
            result = service.run_batch(temp_path)
            assert result is not None
        finally:
            Path(temp_path).unlink()

    def test_resources_configuration_applied(self, service, valid_batch_config):
        """Test that resources configuration is properly applied."""
        valid_batch_config["resources"]["parallel"]["max_workers"] = 8
        valid_batch_config["resources"]["parallel"]["internal_jobs"] = 4

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            # This should process resources configuration without error
            result = service.run_batch(temp_path)
            assert result is not None
        finally:
            Path(temp_path).unlink()

    def test_logging_configuration_applied(self, service, valid_batch_config):
        """Test that logging configuration is properly applied."""
        valid_batch_config["logging"]["level"] = "DEBUG"
        valid_batch_config["logging"]["output"] = {"console": True, "file": False}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(valid_batch_config, f)
            temp_path = f.name

        try:
            # This should process logging configuration without error
            result = service.run_batch(temp_path)
            assert result is not None
        finally:
            Path(temp_path).unlink()


class TestExperimentServiceOptimization:
    """Tests for optimization functionality."""

    @pytest.fixture
    def service(self):
        """Service with optimization support."""
        dataset_repo = FakeDatasetRepository()
        algo_registry = FakeAlgorithmRegistry()
        exporter = FakeExportPort()
        executor = FakeExecutorPort()
        monitoring_service = FakeMonitoringService()

        return ExperimentService(
            dataset_repo=dataset_repo,
            exporter=exporter,
            executor=executor,
            algo_registry=algo_registry,
            monitoring_service=monitoring_service,
        )

    @pytest.fixture
    def optimization_batch_config(self):
        """Optimization batch configuration."""
        return {
            "metadados": {
                "nome": "Optimization Test",
                "descricao": "Test optimization",
                "autor": "Test",
                "versao": "1.0",
                "data_criacao": "2025-01-28",
                "tags": ["optimization"],
                "timeout_global": 3600,
            },
            "datasets": [
                {
                    "id": "test_dataset",
                    "nome": "Test Dataset",
                    "tipo": "synthetic",
                    "parametros": {"n": 10, "L": 20, "alphabet": "ACGT", "seed": 42},
                }
            ],
            "algorithms": [
                {
                    "id": "test_config",
                    "nome": "Test Config",
                    "algorithms": ["BLF-GA"],
                    "algorithm_params": {"BLF-GA": {"pop_size": 100, "max_gens": 200}},
                }
            ],
            "task": {"type": "optimization"},
            "optimization": {
                "method": "optuna",
                "optimizations": [
                    {
                        "nome": "Test Optimization",
                        "study_name": "test_study",
                        "direction": "minimize",
                        "n_trials": 10,
                        "timeout_per_trial": 60,
                        "target_datasets": ["test_dataset"],
                        "target_algorithm": "test_config",
                        "parameters": {
                            "BLF-GA": {
                                "pop_size": {"type": "int", "low": 50, "high": 200}
                            }
                        },
                    }
                ],
            },
            "export": {"enabled": True, "destination": "outputs/opt_test"},
        }

    def test_optimization_batch_execution(self, service, optimization_batch_config):
        """Test optimization batch execution."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(optimization_batch_config, f)
            temp_path = f.name

        try:
            result = service.run_batch(temp_path)
            assert result is not None
            assert "total_experiments" in result
        finally:
            Path(temp_path).unlink()


class TestExperimentServiceSensitivity:
    """Tests for sensitivity analysis functionality."""

    @pytest.fixture
    def service(self):
        """Service with sensitivity support."""
        dataset_repo = FakeDatasetRepository()
        algo_registry = FakeAlgorithmRegistry()
        exporter = FakeExportPort()
        executor = FakeExecutorPort()
        monitoring_service = FakeMonitoringService()

        return ExperimentService(
            dataset_repo=dataset_repo,
            exporter=exporter,
            executor=executor,
            algo_registry=algo_registry,
            monitoring_service=monitoring_service,
        )

    @pytest.fixture
    def sensitivity_batch_config(self):
        """Sensitivity batch configuration."""
        return {
            "metadados": {
                "nome": "Sensitivity Test",
                "descricao": "Test sensitivity analysis",
                "autor": "Test",
                "versao": "1.0",
                "data_criacao": "2025-01-28",
                "tags": ["sensitivity"],
                "timeout_global": 3600,
            },
            "datasets": [
                {
                    "id": "test_dataset",
                    "nome": "Test Dataset",
                    "tipo": "synthetic",
                    "parametros": {"n": 10, "L": 20, "alphabet": "ACGT", "seed": 42},
                }
            ],
            "algorithms": [
                {
                    "id": "test_config",
                    "nome": "Test Config",
                    "algorithms": ["BLF-GA"],
                    "algorithm_params": {"BLF-GA": {"pop_size": 100, "max_gens": 200}},
                }
            ],
            "task": {"type": "sensitivity"},
            "sensitivity": {
                "method": "SALib",
                "analyses": [
                    {
                        "nome": "Test Analysis",
                        "analysis_method": "morris",
                        "target_datasets": ["test_dataset"],
                        "target_algorithm": "BLF-GA",
                        "n_samples": 100,
                        "repetitions_per_sample": 2,
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
            },
            "export": {"enabled": True, "destination": "outputs/sens_test"},
        }

    def test_sensitivity_batch_execution(self, service, sensitivity_batch_config):
        """Test sensitivity batch execution."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(sensitivity_batch_config, f)
            temp_path = f.name

        try:
            result = service.run_batch(temp_path)
            assert result is not None
            assert "total_experiments" in result
        finally:
            Path(temp_path).unlink()
