"""
Domain-based Algorithm Registry

Implements AlgorithmRegistry using the global domain registry.
"""

from typing import Any, Dict, List, Type

from src.domain import CSPAlgorithm, global_registry
from src.domain.errors import AlgorithmNotFoundError, AlgorithmRegistrationError


class DomainAlgorithmRegistry:
    """Algorithm registry based on the global domain registry."""

    def register_algorithm(self, algorithm_class: Type[CSPAlgorithm]) -> None:
        """Register an algorithm."""
        if not issubclass(algorithm_class, CSPAlgorithm):
            raise AlgorithmRegistrationError(
                f"Algorithm must inherit from CSPAlgorithm: {algorithm_class}"
            )

        # Use class name as key
        name = algorithm_class.__name__
        global_registry[name] = algorithm_class

    def get_algorithm(self, name: str) -> Type[CSPAlgorithm]:
        """Get algorithm class by name."""
        if name not in global_registry:
            raise AlgorithmNotFoundError(f"Algorithm not found: {name}")

        return global_registry[name]

    def list_algorithms(self) -> List[str]:
        """List available algorithms."""
        return list(global_registry.keys())

    def algorithm_exists(self, name: str) -> bool:
        """Check if algorithm exists."""
        return name in global_registry

    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """Get algorithm metadata."""
        if name not in global_registry:
            raise AlgorithmNotFoundError(f"Algorithm not found: {name}")

        algo_class = global_registry[name]
        return {
            "name": name,
            "class": algo_class.__name__,
            "module": algo_class.__module__,
            "description": algo_class.__doc__ or "No description",
        }

    # Compatibility methods for existing code
    def register(self, name: str, algorithm_class: Type[CSPAlgorithm]) -> None:
        """Register an algorithm (legacy method)."""
        if not issubclass(algorithm_class, CSPAlgorithm):
            raise AlgorithmRegistrationError(
                f"Algorithm must inherit from CSPAlgorithm: {algorithm_class}"
            )

        global_registry[name] = algorithm_class

    def get(self, name: str) -> Type[CSPAlgorithm]:
        """Get algorithm class by name (legacy method)."""
        return self.get_algorithm(name)

    def list_available(self) -> List[str]:
        """List available algorithms."""
        return list(global_registry.keys())

    def exists(self, name: str) -> bool:
        """Check if algorithm exists."""
        return name in global_registry

    def get_info(self, name: str) -> Dict[str, Any]:
        """Get information about algorithm."""
        if name not in global_registry:
            raise AlgorithmNotFoundError(f"Algorithm not found: {name}")

        algorithm_class = global_registry[name]

        return {
            "name": getattr(algorithm_class, "name", name),
            "default_params": getattr(algorithm_class, "default_params", {}),
            "is_deterministic": getattr(algorithm_class, "is_deterministic", False),
            "supports_internal_parallel": getattr(
                algorithm_class, "supports_internal_parallel", False
            ),
            "description": algorithm_class.__doc__ or "No description available",
        }

    def create_instance(
        self, name: str, strings: List[str], alphabet: str, **params
    ) -> CSPAlgorithm:
        """Create algorithm instance."""
        algorithm_class = self.get(name)
        return algorithm_class(strings, alphabet, **params)
