"""
Configurações padrão para o algoritmo DP-CSP.

Atributos:
    DP_CSP_DEFAULTS (dict): Parâmetros padrão do DP-CSP.
"""

# DP-CSP Configuration
DP_CSP_DEFAULTS = {
    "max_d": None,  # padrão: distância baseline
    "warn_threshold": 9,  # alerta se (d+1)^n > 10^9
    "max_time": 300,  # timeout em segundos
}
