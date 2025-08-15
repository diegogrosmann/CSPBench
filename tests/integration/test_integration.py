"""
Testes de integração para o sistema completo.

Testa a integração entre diferentes componentes do sistema.
"""

import json
import os
import tempfile
from pathlib import Path

import pytest
import yaml

from src.application.services.experiment_service import ExperimentService
from src.domain import Dataset, SyntheticDatasetGenerator

# Importa LocalExecutor via pacote raiz (alias para Executor)
from src.infrastructure.orchestrators import LocalExecutor


class TestEndToEndExecution:
    """Tests for end-to-end experiment execution."""

    @pytest.fixture
    def test_dataset_file(self):
        """Create a temporary test dataset file."""
        sequences = ["ACGTACGTACGT", "ATGTACGTACGT", "ACGTCCGTACGT", "ACGTACGTCCGT"]

        fasta_content = ""
        for i, seq in enumerate(sequences):
            fasta_content += f">seq{i+1}\n{seq}\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(fasta_content)
            f.flush()
            yield f.name

        os.unlink(f.name)

    @pytest.fixture
    def simple_batch_config(self, test_dataset_file):
        """Create a simple batch configuration."""
        config = {
            "metadata": {
                "name": "Integration Test Batch",
                "description": "Testing end-to-end execution",
            },
            "infrastructure": {"executor": {"type": "local", "max_workers": 2}},
            "export": {
                "enabled": True,
                "formats": ["json"],
                "output_dir": "test_output",
            },
            "experiments": [
                {
                    "name": "baseline_test",
                    "dataset": {"path": test_dataset_file},
                    "algorithm": {
                        "name": "baseline",
                        "parameters": {"tie_break": "lex"},
                    },
                },
                {
                    "name": "blf_ga_test",
                    "dataset": {"path": test_dataset_file},
                    "algorithm": {
                        "name": "blf_ga",
                        "parameters": {
                            "pop_size": 10,
                            "max_gens": 5,
                            "cross_prob": 0.8,
                            "mut_prob": 0.1,
                        },
                    },
                },
            ],
        }
        return config

    @pytest.fixture
    def batch_config_file(self, simple_batch_config):
        """Create a temporary batch configuration file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(simple_batch_config, f)
            f.flush()
            yield f.name

        os.unlink(f.name)

    def test_complete_experiment_execution(self, batch_config_file):
        """Test complete experiment execution from config file."""
        # Use fake services from unit tests
        from tests.unit.application.fakes import (
            FakeAlgorithmRegistry,
            FakeDatasetRepository,
            FakeExecutorPort,
            FakeExportPort,
        )

        dataset_repo = FakeDatasetRepository()
        exporter = FakeExportPort()
        executor = FakeExecutorPort()
        algo_registry = FakeAlgorithmRegistry()

        service = ExperimentService(
            dataset_repo=dataset_repo,
            exporter=exporter,
            executor=executor,
            algo_registry=algo_registry,
        )

        # Execute batch
        results = service.run_batch(batch_config_file)

        # Verify results
        assert isinstance(results, dict)
        assert "total_experiments" in results

    def test_experiment_with_monitoring(self, simple_batch_config):
        """Test experiment execution with monitoring enabled."""
        # Add monitoring to config
        simple_batch_config["monitoring"] = {
            "enabled": True,
            "progress_bar": True,
            "log_level": "INFO",
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(simple_batch_config, f)
            f.flush()

            try:
                executor = LocalExecutor()
                service = ExperimentService(executor)

                results = service.execute_batch_from_file(f.name)

                assert len(results) == 2
                for result in results:
                    assert "monitoring" in result["metadata"]

            finally:
                os.unlink(f.name)

    def test_experiment_with_export(self, simple_batch_config, test_dataset_file):
        """Test experiment execution with result export."""
        # Configure export
        simple_batch_config["export"]["formats"] = ["json", "csv"]

        with tempfile.TemporaryDirectory() as temp_dir:
            simple_batch_config["export"]["output_dir"] = temp_dir

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".yaml", delete=False
            ) as f:
                yaml.dump(simple_batch_config, f)
                f.flush()

                try:
                    executor = LocalExecutor()
                    service = ExperimentService(executor)

                    results = service.execute_batch_from_file(f.name)

                    # Verify export files were created
                    output_files = list(Path(temp_dir).glob("*"))
                    assert len(output_files) > 0

                    # Check for expected file types
                    json_files = list(Path(temp_dir).glob("*.json"))
                    csv_files = list(Path(temp_dir).glob("*.csv"))

                    assert len(json_files) > 0
                    assert len(csv_files) > 0

                finally:
                    os.unlink(f.name)

    def test_experiment_error_handling(self, test_dataset_file):
        """Test experiment error handling."""
        config = {
            "metadata": {"name": "Error Test"},
            "infrastructure": {"executor": {"type": "local"}},
            "experiments": [
                {
                    "name": "invalid_algorithm",
                    "dataset": {"path": test_dataset_file},
                    "algorithm": {"name": "nonexistent_algorithm", "parameters": {}},
                }
            ],
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config, f)
            f.flush()

            try:
                executor = LocalExecutor()
                service = ExperimentService(executor)

                # Should handle errors gracefully
                with pytest.raises(Exception):  # Expect some kind of error
                    service.execute_batch_from_file(f.name)

            finally:
                os.unlink(f.name)

    def test_parallel_execution(self, simple_batch_config):
        """Test parallel execution of experiments."""
        # Add more experiments for parallel testing
        simple_batch_config["infrastructure"]["executor"]["max_workers"] = 2

        # Duplicate experiments to test parallelism
        original_experiments = simple_batch_config["experiments"].copy()
        for i in range(3):
            for exp in original_experiments:
                new_exp = exp.copy()
                new_exp["name"] = f"{exp['name']}_copy_{i}"
                simple_batch_config["experiments"].append(new_exp)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(simple_batch_config, f)
            f.flush()

            try:
                executor = LocalExecutor()
                service = ExperimentService(executor)

                results = service.execute_batch_from_file(f.name)

                # Should have results for all experiments
                assert len(results) == 8  # 2 original + 6 copies

                # All experiments should complete successfully
                for result in results:
                    assert "result" in result
                    assert "distance" in result
                    assert result["distance"] >= 0

            finally:
                os.unlink(f.name)


class TestConfigurationIntegration:
    """Tests for configuration integration."""

    def test_legacy_format_compatibility(self):
        """Test compatibility with legacy configuration format."""
        legacy_config = {
            "algorithms": ["baseline", "blf_ga"],
            "datasets": ["test_dataset.fasta"],
            "parameters": {
                "baseline": {"tie_break": "lex"},
                "blf_ga": {"pop_size": 10, "max_gens": 5},
            },
            "output_dir": "legacy_output",
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(legacy_config, f)
            f.flush()

            try:
                from src.application.services.config_parser import ConfigurationParser

                parser = ConfigurationParser()

                # Should handle legacy format
                config = parser.parse_file(f.name)

                assert config is not None
                assert hasattr(config, "experiments")

            finally:
                os.unlink(f.name)

    def test_json_configuration(self):
        """Test JSON configuration format support."""
        config = {
            "metadata": {"name": "JSON Test"},
            "infrastructure": {"executor": {"type": "local"}},
            "export": {"enabled": True, "formats": ["json"]},
            "experiments": [
                {
                    "name": "json_test",
                    "dataset": {"path": "test.fasta"},
                    "algorithm": {"name": "baseline", "parameters": {}},
                }
            ],
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config, f)
            f.flush()

            try:
                from src.application.services.config_parser import ConfigurationParser

                parser = ConfigurationParser()

                parsed_config = parser.parse_file(f.name)

                assert parsed_config is not None
                assert parsed_config.metadata.name == "JSON Test"
                assert len(parsed_config.experiments) == 1

            finally:
                os.unlink(f.name)

    def test_configuration_validation(self):
        """Test configuration validation."""
        invalid_configs = [
            # Missing required fields
            {"metadata": {"name": "Test"}},
            # Invalid algorithm
            {
                "metadata": {"name": "Test"},
                "experiments": [
                    {
                        "name": "test",
                        "dataset": {"path": "test.fasta"},
                        "algorithm": {"name": "invalid_algo"},
                    }
                ],
            },
            # Invalid parameters
            {
                "metadata": {"name": "Test"},
                "experiments": [
                    {
                        "name": "test",
                        "dataset": {"path": "test.fasta"},
                        "algorithm": {
                            "name": "blf_ga",
                            "parameters": {"pop_size": -1},  # Invalid negative value
                        },
                    }
                ],
            },
        ]

        from src.application.services.config_parser import ConfigurationParser

        parser = ConfigurationParser()

        for invalid_config in invalid_configs:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".yaml", delete=False
            ) as f:
                yaml.dump(invalid_config, f)
                f.flush()

                try:
                    with pytest.raises(Exception):
                        parser.parse_file(f.name)
                finally:
                    os.unlink(f.name)


class TestDatasetIntegration:
    """Tests for dataset integration."""

    def test_dataset_loading_with_algorithms(self):
        """Test dataset loading and algorithm execution integration."""
        # Create test dataset
        dataset = SyntheticDatasetGenerator.generate_random(
            n=5, length=8, alphabet="ACGT", seed=42
        )

        # Test with different algorithms
        from algorithms.baseline.algorithm import BaselineAlgorithm
        from algorithms.blf_ga.algorithm import BLFGAAlgorithm

        algorithms = [
            BaselineAlgorithm(dataset.sequences, "ACGT"),
            BLFGAAlgorithm(dataset.sequences, "ACGT", pop_size=10, max_gens=5),
        ]

        for algorithm in algorithms:
            result, distance, metadata = algorithm.run()

            # Verify integration
            assert len(result) == dataset.length
            assert all(c in dataset.alphabet for c in result)
            assert distance >= 0
            assert isinstance(metadata, dict)

    def test_fasta_file_integration(self):
        """Test FASTA file loading and processing integration."""
        sequences = ["ACGTACGT", "ATGTACGT", "ACGTCCGT"]

        # Create FASTA file
        fasta_content = ""
        for i, seq in enumerate(sequences):
            fasta_content += f">sequence_{i+1}\n{seq}\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(fasta_content)
            f.flush()

            try:
                # Load dataset
                dataset = Dataset.from_file(f.name)

                # Verify loading
                assert dataset.sequences == sequences
                assert dataset.size == 3
                assert dataset.length == 8

                # Test with algorithm
                from algorithms.baseline.algorithm import BaselineAlgorithm

                algorithm = BaselineAlgorithm(dataset.sequences, "ACGT")
                result, distance, metadata = algorithm.run()

                assert isinstance(result, str)
                assert len(result) == 8

            finally:
                os.unlink(f.name)

    def test_synthetic_dataset_generation_integration(self):
        """Test synthetic dataset generation with full workflow."""
        # Generate datasets with different properties
        datasets = [
            SyntheticDatasetGenerator.generate_random(
                n=4, length=6, alphabet="ACGT", seed=42
            ),
            SyntheticDatasetGenerator.generate_with_noise(
                base_sequence="ACGTACGT", n=4, noise_rate=0.2, alphabet="ACGT", seed=42
            ),
            SyntheticDatasetGenerator.generate_clustered(
                n=6, length=8, num_clusters=2, alphabet="ACGT", seed=42
            ),
        ]

        from algorithms.baseline.algorithm import BaselineAlgorithm

        for dataset in datasets:
            # Test algorithm execution on each dataset type
            algorithm = BaselineAlgorithm(dataset.sequences, "ACGT")
            result, distance, metadata = algorithm.run()

            # Verify results are consistent with dataset
            assert len(result) == dataset.length
            assert all(c in dataset.alphabet for c in result)
            assert distance >= 0

            # Test metrics calculation
            from src.domain.metrics import calculate_coverage, calculate_entropy

            coverage = calculate_coverage(dataset.sequences, result)
            entropy = calculate_entropy(dataset.sequences)

            assert coverage["max_distance"] == distance
            assert len(entropy["position_entropies"]) == dataset.length


class TestPerformanceIntegration:
    """Tests for performance and scalability."""

    def test_small_scale_performance(self):
        """Test performance with small datasets."""
        import time

        dataset = SyntheticDatasetGenerator.generate_random(
            n=10, length=20, alphabet="ACGT", seed=42
        )

        from algorithms.baseline.algorithm import BaselineAlgorithm
        from algorithms.blf_ga.algorithm import BLFGAAlgorithm

        algorithms = [
            BaselineAlgorithm(dataset.sequences, "ACGT"),
            BLFGAAlgorithm(dataset.sequences, "ACGT", pop_size=20, max_gens=10),
        ]

        for algorithm in algorithms:
            start_time = time.time()
            result, distance, metadata = algorithm.run()
            execution_time = time.time() - start_time

            # Should complete reasonably quickly
            assert execution_time < 10.0  # 10 seconds max
            assert "execution_time" in metadata
            assert metadata["execution_time"] > 0

    def test_memory_usage_integration(self):
        """Test memory usage with different dataset sizes."""
        import gc

        # Test with progressively larger datasets
        sizes = [5, 10, 20]

        for n in sizes:
            dataset = SyntheticDatasetGenerator.generate_random(
                n=n, length=10, alphabet="ACGT", seed=42
            )

            from algorithms.baseline.algorithm import BaselineAlgorithm

            algorithm = BaselineAlgorithm(dataset.sequences, "ACGT")

            # Force garbage collection before test
            gc.collect()

            result, distance, metadata = algorithm.run()

            # Verify results regardless of size
            assert len(result) == 10
            assert distance >= 0

            # Clean up
            del algorithm, result, distance, metadata
            gc.collect()

    def test_concurrent_execution(self):
        """Test concurrent execution of multiple experiments."""
        import concurrent.futures
        import time

        def run_experiment(exp_id):
            dataset = SyntheticDatasetGenerator.generate_random(
                n=5, length=8, alphabet="ACGT", seed=exp_id
            )

            from algorithms.baseline.algorithm import BaselineAlgorithm

            algorithm = BaselineAlgorithm(dataset.sequences, "ACGT")

            result, distance, metadata = algorithm.run()
            return {
                "experiment_id": exp_id,
                "result": result,
                "distance": distance,
                "execution_time": metadata.get("execution_time", 0),
            }

        # Run multiple experiments concurrently
        start_time = time.time()

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            futures = [executor.submit(run_experiment, i) for i in range(4)]

            results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]

        total_time = time.time() - start_time

        # Verify all experiments completed
        assert len(results) == 4

        for result in results:
            assert "experiment_id" in result
            assert "result" in result
            assert "distance" in result
            assert result["distance"] >= 0

        # Concurrent execution should be faster than sequential
        # (This is a rough check, actual speedup depends on system)
        assert total_time < 20.0  # Should complete within reasonable time


class TestRegressionIntegration:
    """Regression tests for integration scenarios."""

    def test_algorithm_consistency(self):
        """Test that algorithms produce consistent results."""
        # Use fixed seed for reproducibility
        dataset = SyntheticDatasetGenerator.generate_random(
            n=5, length=10, alphabet="ACGT", seed=12345
        )

        from algorithms.baseline.algorithm import BaselineAlgorithm

        # Run same algorithm multiple times
        results = []
        for _ in range(3):
            algorithm = BaselineAlgorithm(dataset.sequences, "ACGT", tie_break="lex")
            result, distance, metadata = algorithm.run()
            results.append((result, distance))

        # All results should be identical for deterministic algorithm
        first_result = results[0]
        for result in results[1:]:
            assert result == first_result

    def test_configuration_format_stability(self):
        """Test that configuration format parsing is stable."""
        config_variations = [
            # Minimal config
            {
                "metadata": {"name": "Minimal"},
                "experiments": [
                    {
                        "name": "test",
                        "dataset": {"path": "test.fasta"},
                        "algorithm": {"name": "baseline"},
                    }
                ],
            },
            # Full config
            {
                "metadata": {
                    "name": "Full Test",
                    "description": "Complete configuration",
                },
                "infrastructure": {"executor": {"type": "local", "max_workers": 1}},
                "export": {
                    "enabled": True,
                    "formats": ["json"],
                    "output_dir": "test_output",
                },
                "monitoring": {"enabled": True, "progress_bar": True},
                "experiments": [
                    {
                        "name": "full_test",
                        "dataset": {"path": "test.fasta"},
                        "algorithm": {
                            "name": "blf_ga",
                            "parameters": {"pop_size": 20, "max_gens": 10},
                        },
                    }
                ],
            },
        ]

        from src.application.services.config_parser import ConfigurationParser

        parser = ConfigurationParser()

        for config in config_variations:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".yaml", delete=False
            ) as f:
                yaml.dump(config, f)
                f.flush()

                try:
                    parsed_config = parser.parse_file(f.name)

                    # Basic validation
                    assert parsed_config is not None
                    assert hasattr(parsed_config, "metadata")
                    assert hasattr(parsed_config, "experiments")
                    assert len(parsed_config.experiments) > 0

                finally:
                    os.unlink(f.name)

    def test_error_recovery_integration(self):
        """Test error recovery in integrated scenarios."""
        # Test scenarios that should handle errors gracefully
        problematic_scenarios = [
            # Missing dataset file
            {
                "metadata": {"name": "Missing Dataset"},
                "experiments": [
                    {
                        "name": "test",
                        "dataset": {"path": "nonexistent.fasta"},
                        "algorithm": {"name": "baseline"},
                    }
                ],
            }
        ]

        for scenario in problematic_scenarios:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".yaml", delete=False
            ) as f:
                yaml.dump(scenario, f)
                f.flush()

                try:
                    executor = LocalExecutor()
                    service = ExperimentService(executor)

                    # Should handle errors without crashing
                    with pytest.raises(Exception):
                        service.execute_batch_from_file(f.name)

                finally:
                    os.unlink(f.name)
