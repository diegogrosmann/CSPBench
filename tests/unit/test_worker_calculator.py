"""
Testes Unitários para Módulo de Calculadora de Workers - CSPBench

Testa as funcionalidades do módulo src.utils.worker_calculator
"""

from unittest.mock import MagicMock, patch

import pytest

from src.utils.worker_calculator import (
    calculate_internal_workers,
    calculate_optuna_workers,
    calculate_salib_workers,
    calculate_workers,
    get_cpu_count,
    get_worker_config,
)


class TestGetCpuCount:
    """Testes para função get_cpu_count."""

    def test_get_cpu_count_positive(self):
        """Testa se retorna número positivo de CPUs."""
        cpu_count = get_cpu_count()
        assert cpu_count > 0
        assert isinstance(cpu_count, int)

    def test_get_cpu_count_fallback(self):
        """Testa que get_cpu_count sempre retorna um valor válido."""
        cpu_count = get_cpu_count()
        assert cpu_count >= 1  # Pelo menos 1 CPU sempre


class TestCalculateOptunaWorkers:
    """Testes para função calculate_optuna_workers."""

    def test_calculate_optuna_workers_default(self):
        """Testa cálculo padrão de workers para Optuna."""
        optuna_workers, internal_workers = calculate_optuna_workers()
        assert optuna_workers > 0
        assert internal_workers > 0
        assert isinstance(optuna_workers, int)
        assert isinstance(internal_workers, int)

    def test_calculate_optuna_workers_with_algorithm(self):
        """Testa cálculo com algoritmo específico."""
        optuna_workers, internal_workers = calculate_optuna_workers(
            algorithm_name="BLF-GA"
        )
        assert optuna_workers > 0
        assert internal_workers > 0
        assert isinstance(optuna_workers, int)
        assert isinstance(internal_workers, int)

    def test_calculate_optuna_workers_with_cpu_count(self):
        """Testa cálculo com número específico de CPUs."""
        optuna_workers, internal_workers = calculate_optuna_workers(total_cpus=4)
        assert optuna_workers > 0
        assert internal_workers > 0
        assert isinstance(optuna_workers, int)
        assert isinstance(internal_workers, int)


class TestCalculateSalibWorkers:
    """Testes para função calculate_salib_workers."""

    def test_calculate_salib_workers_default(self):
        """Testa cálculo padrão de workers para SALib."""
        workers = calculate_salib_workers()
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_salib_workers_with_samples(self):
        """Testa cálculo com número específico de amostras."""
        workers = calculate_salib_workers(n_samples=1000)
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_salib_workers_with_cpu_count(self):
        """Testa cálculo com número específico de CPUs."""
        workers = calculate_salib_workers(total_cpus=8)
        assert workers > 0
        assert isinstance(workers, int)


class TestCalculateInternalWorkers:
    """Testes para função calculate_internal_workers."""

    def test_calculate_internal_workers_default(self):
        """Testa cálculo padrão de workers internos."""
        workers = calculate_internal_workers("BLF-GA", 4)
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_internal_workers_with_algorithm(self):
        """Testa cálculo com algoritmo específico."""
        workers = calculate_internal_workers("BLF-GA", 8)
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_internal_workers_with_external_parallelism(self):
        """Testa cálculo com paralelismo externo."""
        workers = calculate_internal_workers("BLF-GA", 8, external_parallelism=True)
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_internal_workers_sequential_algorithm(self):
        """Testa cálculo para algoritmo sequencial."""
        workers = calculate_internal_workers("CSC", 8)
        assert workers > 0
        assert isinstance(workers, int)


class TestGetWorkerConfig:
    """Testes para função get_worker_config."""

    def test_get_worker_config_default(self):
        """Testa configuração padrão de workers."""
        config = get_worker_config()

        assert isinstance(config, dict)
        assert "optuna_workers" in config
        assert "internal_workers" in config
        assert "total_cpus" in config

    def test_get_worker_config_with_algorithm(self):
        """Testa configuração com algoritmo específico."""
        config = get_worker_config(algorithm_name="BLF-GA")

        assert isinstance(config, dict)
        assert config["optuna_workers"] > 0
        assert config["internal_workers"] > 0

    def test_get_worker_config_optuna_context(self):
        """Testa configuração para contexto Optuna."""
        config = get_worker_config(context="optuna", algorithm_name="BLF-GA")

        assert isinstance(config, dict)
        assert config["optuna_workers"] > 0
        assert config["internal_workers"] > 0

    def test_get_worker_config_salib_context(self):
        """Testa configuração para contexto SALib."""
        config = get_worker_config(
            context="salib", algorithm_name="BLF-GA", n_samples=1000
        )

        assert isinstance(config, dict)
        assert config["salib_workers"] > 0


class TestCalculateWorkers:
    """Testes para função calculate_workers."""

    def test_calculate_workers_default(self):
        """Testa cálculo padrão de workers."""
        workers = calculate_workers()
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_workers_blf_ga(self):
        """Testa cálculo de workers para BLF-GA."""
        workers = calculate_workers(algorithm_name="BLF-GA")
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_workers_csc(self):
        """Testa cálculo de workers para CSC."""
        workers = calculate_workers(algorithm_name="CSC")
        assert workers > 0
        assert isinstance(workers, int)

    def test_calculate_workers_unknown_algorithm(self):
        """Testa cálculo de workers para algoritmo desconhecido."""
        workers = calculate_workers(algorithm_name="UNKNOWN")
        assert workers > 0
        assert isinstance(workers, int)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
