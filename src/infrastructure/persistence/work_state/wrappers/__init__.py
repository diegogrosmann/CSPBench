"""Persistence wrapper classes.

Este módulo expõe as classes wrapper de persistência com suporte a um alias
retrocompatível: ``UnitScopedPersistence`` que aponta para
``ExecutionScopedPersistence``. O nome antigo ainda é importado em diversos
locais do código e testes; manter o alias evita quebrar estes pontos enquanto
o restante do refactor avança.
"""

from .combination_scoped import CombinationScopedPersistence
from .execution_scoped import ExecutionScopedPersistence
from .work_scoped import WorkScopedPersistence

# Alias backward-compatible
UnitScopedPersistence = ExecutionScopedPersistence  # type: ignore

__all__ = [
    "WorkScopedPersistence",
    "ExecutionScopedPersistence",
    "UnitScopedPersistence",  # alias legado
    "CombinationScopedPersistence",
]
