"""
Configurações padrão de limites de recursos para o monitoramento.

Atributos:
    RESOURCE_LIMITS_DEFAULTS (dict): Limites padrão de memória, iterações, etc.
    MEMORY_SAFETY_DEFAULTS (dict): Configurações de segurança de memória.
    GC_DEFAULTS (dict): Configurações de garbage collection.
"""

import logging
import os

logger = logging.getLogger(__name__)

# Configurações básicas de recursos
RESOURCE_LIMITS_DEFAULTS = {
    "max_memory_mb": 2048,
    "max_iterations": 100000,
    "check_interval": 2.0,
    "gc_threshold": 1000,
}

# Configurações de segurança de memória
MEMORY_SAFETY_DEFAULTS = {
    "memory_usage_ratio": 0.5,  # Usar no máximo 50% da memória disponível
    "max_safe_memory_mb": 4096,  # Máximo absoluto
    "min_free_memory_mb": 512,  # Mínimo a deixar livre
    "memory_check_aggressive": True,  # Verificação agressiva
}

# Configurações de garbage collection
GC_DEFAULTS = {
    "gc_auto_collect": True,  # Coleta automática
    "gc_frequency": 500,  # A cada 500 checks
    "gc_force_on_limit": True,  # Forçar GC ao aproximar do limite
    "gc_threshold_ratio": 0.8,  # Forçar GC a 80% do limite
}


def get_resource_limits_from_env() -> dict:
    """
    Obtém limites de recursos de variáveis de ambiente.

    Variáveis suportadas:
    - CSP_MAX_MEMORY_MB: Limite de memória em MB
    - CSP_MAX_ITERATIONS: Limite de iterações
    - CSP_CHECK_INTERVAL: Intervalo de verificação
    - CSP_GC_THRESHOLD: Limite para GC
    """
    env_limits = {}

    # Mapear variáveis de ambiente para configurações
    env_mappings = {
        "CSP_MAX_MEMORY_MB": ("max_memory_mb", int),
        "CSP_MAX_ITERATIONS": ("max_iterations", int),
        "CSP_CHECK_INTERVAL": ("check_interval", float),
        "CSP_GC_THRESHOLD": ("gc_threshold", int),
    }

    for env_var, (config_key, type_func) in env_mappings.items():
        if env_var in os.environ:
            try:
                env_limits[config_key] = type_func(os.environ[env_var])
                logger.info(
                    f"Configuração de ambiente: {config_key} = {env_limits[config_key]}"
                )
            except ValueError as e:
                logger.warning(f"Valor inválido para {env_var}: {e}")

    return env_limits


def get_merged_resource_limits() -> dict:
    """
    Obtém configurações de recursos mergeadas (defaults + ambiente).
    """
    limits = RESOURCE_LIMITS_DEFAULTS.copy()
    env_limits = get_resource_limits_from_env()
    limits.update(env_limits)
    return limits


def get_merged_memory_safety() -> dict:
    """
    Obtém configurações de segurança de memória mergeadas.
    """
    return MEMORY_SAFETY_DEFAULTS.copy()


def get_merged_gc_config() -> dict:
    """
    Obtém configurações de garbage collection mergeadas.
    """
    return GC_DEFAULTS.copy()
