"""H2-CSP algorithm package.

Importa a classe ``H2CSPAlgorithm`` que, ao ser carregada, registra-se no
``global_registry`` via o decorator ``@register_algorithm``.

Mantido simples para evitar dependências circulares.
"""

from .algorithm import H2CSPAlgorithm  # noqa: F401  (re-export)
