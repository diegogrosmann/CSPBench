"""
Registry de Algoritmos baseado no domínio

Implementa AlgorithmRegistry usando o registry global do domínio.
"""

from typing import Any, Dict, List, Type

from src.domain import CSPAlgorithm, global_registry
from src.domain.errors import AlgorithmNotFoundError, AlgorithmRegistrationError


class DomainAlgorithmRegistry:
    """Registry de algoritmos baseado no registry global do domínio."""

    def register_algorithm(self, algorithm_class: Type[CSPAlgorithm]) -> None:
        """Registra um algoritmo."""
        if not issubclass(algorithm_class, CSPAlgorithm):
            raise AlgorithmRegistrationError(
                f"Algoritmo deve herdar de CSPAlgorithm: {algorithm_class}"
            )

        # Usa o nome da classe como chave
        name = algorithm_class.__name__
        global_registry[name] = algorithm_class

    def get_algorithm(self, name: str) -> Type[CSPAlgorithm]:
        """Obtém classe de algoritmo por nome."""
        if name not in global_registry:
            raise AlgorithmNotFoundError(f"Algoritmo não encontrado: {name}")

        return global_registry[name]

    def list_algorithms(self) -> List[str]:
        """Lista algoritmos disponíveis."""
        return list(global_registry.keys())

    def algorithm_exists(self, name: str) -> bool:
        """Verifica se algoritmo existe."""
        return name in global_registry

    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """Obtém metadados de algoritmo."""
        if name not in global_registry:
            raise AlgorithmNotFoundError(f"Algoritmo não encontrado: {name}")

        algo_class = global_registry[name]
        return {
            "name": name,
            "class": algo_class.__name__,
            "module": algo_class.__module__,
            "description": algo_class.__doc__ or "Sem descrição",
        }

    # Métodos de compatibilidade para código existente
    def register(self, name: str, algorithm_class: Type[CSPAlgorithm]) -> None:
        """Registra um algoritmo (método legacy)."""
        if not issubclass(algorithm_class, CSPAlgorithm):
            raise AlgorithmRegistrationError(
                f"Algoritmo deve herdar de CSPAlgorithm: {algorithm_class}"
            )

        global_registry[name] = algorithm_class

    def get(self, name: str) -> Type[CSPAlgorithm]:
        """Obtém classe de algoritmo por nome (método legacy)."""
        return self.get_algorithm(name)

    def list_available(self) -> List[str]:
        """Lista algoritmos disponíveis."""
        return list(global_registry.keys())

    def exists(self, name: str) -> bool:
        """Verifica se algoritmo existe."""
        return name in global_registry

    def get_info(self, name: str) -> Dict[str, Any]:
        """Obtém informações sobre algoritmo."""
        if name not in global_registry:
            raise AlgorithmNotFoundError(f"Algoritmo não encontrado: {name}")

        algorithm_class = global_registry[name]

        return {
            "name": getattr(algorithm_class, "name", name),
            "default_params": getattr(algorithm_class, "default_params", {}),
            "is_deterministic": getattr(algorithm_class, "is_deterministic", False),
            "supports_internal_parallel": getattr(
                algorithm_class, "supports_internal_parallel", False
            ),
            "description": algorithm_class.__doc__ or "Sem descrição disponível",
        }

    def create_instance(
        self, name: str, strings: List[str], alphabet: str, **params
    ) -> CSPAlgorithm:
        """Cria instância de algoritmo."""
        algorithm_class = self.get(name)
        return algorithm_class(strings, alphabet, **params)
