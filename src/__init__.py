"""
CSPBench - Framework para Closest String Problem

Framework modular para teste e comparação de algoritmos do Closest String Problem.
Implementação baseada em arquitetura hexagonal com separação clara de responsabilidades.
"""

__version__ = "0.1.0"
__author__ = "CSPBench Development Team"

from . import application, domain, infrastructure, presentation

__all__ = ["domain", "application", "infrastructure", "presentation"]
