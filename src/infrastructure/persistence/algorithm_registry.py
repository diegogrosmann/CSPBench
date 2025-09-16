"""
Domain-based Algorithm Registry

Implements AlgorithmRegistry using the global domain registry.
"""

from typing import Any, Dict, List, Type

from src.domain import CSPAlgorithm, global_registry
from src.domain.errors import AlgorithmNotFoundError, AlgorithmRegistrationError


class DomainAlgorithmRegistry:
    """
    Algorithm registry based on the global domain registry.

    This class provides a registry for CSP algorithms, allowing registration,
    retrieval, and management of algorithm implementations. It uses the global
    domain registry as the underlying storage mechanism.
    """

    def register_algorithm(self, algorithm_class: Type[CSPAlgorithm]) -> None:
        """
        Register an algorithm class in the registry.

        Args:
            algorithm_class: CSPAlgorithm subclass to register

        Raises:
            AlgorithmRegistrationError: If the class doesn't inherit from CSPAlgorithm
        """
        if not issubclass(algorithm_class, CSPAlgorithm):
            raise AlgorithmRegistrationError(
                f"Algorithm must inherit from CSPAlgorithm: {algorithm_class}"
            )

        # Use class name as key
        name = algorithm_class.__name__
        global_registry[name] = algorithm_class

    def get_algorithm(self, name: str) -> Type[CSPAlgorithm]:
        """
        Get algorithm class by name.

        Args:
            name: Name of the algorithm to retrieve

        Returns:
            Type[CSPAlgorithm]: The algorithm class

        Raises:
            AlgorithmNotFoundError: If algorithm is not found
        """
        if name not in global_registry:
            raise AlgorithmNotFoundError(f"Algorithm not found: {name}")

        return global_registry[name]

    def list_algorithms(self) -> List[str]:
        """
        List all available algorithm names.

        Returns:
            List[str]: List of registered algorithm names
        """
        return list(global_registry.keys())

    def algorithm_exists(self, name: str) -> bool:
        """
        Check if an algorithm is registered.

        Args:
            name: Name of the algorithm to check

        Returns:
            bool: True if algorithm exists, False otherwise
        """
        return name in global_registry

    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """
        Get metadata information about an algorithm.

        Args:
            name: Name of the algorithm

        Returns:
            Dict[str, Any]: Dictionary containing algorithm metadata

        Raises:
            AlgorithmNotFoundError: If algorithm is not found
        """
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
        """
        Register an algorithm with a custom name (legacy method).

        Args:
            name: Custom name for the algorithm
            algorithm_class: CSPAlgorithm subclass to register

        Raises:
            AlgorithmRegistrationError: If the class doesn't inherit from CSPAlgorithm
        """
        if not issubclass(algorithm_class, CSPAlgorithm):
            raise AlgorithmRegistrationError(
                f"Algorithm must inherit from CSPAlgorithm: {algorithm_class}"
            )

        global_registry[name] = algorithm_class

    def get(self, name: str) -> Type[CSPAlgorithm]:
        """
        Get algorithm class by name (legacy method).

        Args:
            name: Name of the algorithm to retrieve

        Returns:
            Type[CSPAlgorithm]: The algorithm class
        """
        return self.get_algorithm(name)

    def list_available(self) -> List[str]:
        """
        List available algorithms.

        Returns:
            List[str]: List of registered algorithm names
        """
        return list(global_registry.keys())

    def exists(self, name: str) -> bool:
        """
        Check if algorithm exists.

        Args:
            name: Name of the algorithm to check

        Returns:
            bool: True if algorithm exists, False otherwise
        """
        return name in global_registry

    def get_info(self, name: str) -> Dict[str, Any]:
        """
        Get detailed information about an algorithm.

        Args:
            name: Name of the algorithm

        Returns:
            Dict[str, Any]: Dictionary containing algorithm information

        Raises:
            AlgorithmNotFoundError: If algorithm is not found
        """
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
        """
        Create an instance of the specified algorithm.

        Args:
            name: Name of the algorithm to instantiate
            strings: List of strings for the algorithm
            alphabet: Alphabet to use
            **params: Additional parameters for the algorithm

        Returns:
            CSPAlgorithm: Configured algorithm instance
        """
        algorithm_class = self.get(name)
        return algorithm_class(strings, alphabet, **params)
