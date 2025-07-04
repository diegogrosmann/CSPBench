"""
Pacote de algoritmos CSP.

Este módulo inicializa o pacote algorithms e expõe o registry global.
"""

# algorithms package

import importlib
import pkgutil

from .base import global_registry, register_algorithm

# Auto-import subpackages to register algorithms
for _, name, ispkg in pkgutil.iter_modules(__path__):
    if ispkg:
        importlib.import_module(f"{__name__}.{name}")
