"""
Pacote de algoritmos CSP.

Este módulo inicializa o pacote algorithms e expõe o registry global.
"""

# algorithms package

from .base import global_registry, register_algorithm
import pkgutil, importlib

# Auto-import subpackages to register algorithms
for _, name, ispkg in pkgutil.iter_modules(__path__):
    if ispkg:
        importlib.import_module(f"{__name__}.{name}")
