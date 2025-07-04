"""
Configurações padrão de limites de recursos para o monitoramento.

Atributos:
    RESOURCE_LIMITS_DEFAULTS (dict): Limites padrão de memória, iterações, etc.
"""

RESOURCE_LIMITS_DEFAULTS = {
    "max_memory_mb": 2048,
    "max_iterations": 100000,
    "check_interval": 2.0,
    "gc_threshold": 1000,
}
