"""
CSPBench - Framework for Closest String Problem

Modular framework for testing and comparing Closest String Problem algorithms.
Implementation based on hexagonal architecture with clear separation of responsibilities.
"""

__version__ = "0.1.0"
__author__ = "CSPBench Development Team"

from . import application, domain, infrastructure, presentation

__all__ = ["domain", "application", "infrastructure", "presentation"]
