"""
Default configurations for the DP-CSP algorithm.

Attributes:
    DP_CSP_DEFAULTS (dict): Default parameters for DP-CSP.
"""

# DP-CSP Configuration
DP_CSP_DEFAULTS = {
    "max_d": None,  # default: baseline distance
    "warn_threshold": 9,  # warn if (d+1)^n > 10^9
    "max_time": 300,  # timeout in seconds
}
