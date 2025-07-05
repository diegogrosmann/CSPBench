"""
Teste para validar centralização de limites de recursos (T15).
"""

import os
import unittest
from unittest.mock import patch

from src.utils.resource_limits_config import get_merged_resource_limits
from src.utils.resource_monitor import ResourceLimits


class TestResourceLimitsCentralization(unittest.TestCase):
    """Testa a centralização de limites de recursos."""

    def test_resource_limits_from_config(self):
        """Testa que ResourceLimits usa configurações centralizadas."""
        limits = ResourceLimits.from_config()

        # Verificar que os valores padrão são usados
        merged_limits = get_merged_resource_limits()
        self.assertEqual(limits.max_memory_mb, merged_limits["max_memory_mb"])
        self.assertEqual(limits.max_iterations, merged_limits["max_iterations"])
        self.assertEqual(limits.check_interval, merged_limits["check_interval"])

    def test_env_override_propagation(self):
        """Testa que mudanças no ambiente se propagam corretamente."""
        # Simular variável de ambiente
        with patch.dict(os.environ, {"CSP_MAX_MEMORY_MB": "4096"}):
            merged_limits = get_merged_resource_limits()
            self.assertEqual(merged_limits["max_memory_mb"], 4096)

            # Verificar que ResourceLimits usa o valor modificado
            limits = ResourceLimits.from_config()
            self.assertEqual(limits.max_memory_mb, 4096)

    def test_algorithm_executor_uses_centralized_limits(self):
        """Testa que AlgorithmExecutor usa limites centralizados."""
        from src.core.exec.algorithm_executor import AlgorithmExecutor

        # Criar executor
        executor = AlgorithmExecutor(timeout_seconds=30)

        # Verificar que os limites são baseados na configuração central
        base_limits = get_merged_resource_limits()
        expected_max_iterations = min(base_limits["max_iterations"], 30 * 1000)

        self.assertEqual(executor.resource_monitor.limits.max_iterations, expected_max_iterations)


if __name__ == "__main__":
    unittest.main()
