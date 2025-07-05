"""
Módulo de configurações para o sistema CSP-BLFGA.

Este módulo contém utilitários para carregar, validar e gerenciar
configurações de execução em lote.
"""

from .config_loader import ConfigError, ConfigLoader

__all__ = ["ConfigLoader", "ConfigError"]
