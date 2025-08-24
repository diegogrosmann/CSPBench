"""
Default configurations for the DP-CSP algorithm.

Attributes:
    DP_CSP_DEFAULTS (dict): Default parameters for DP-CSP.
"""

# DP-CSP Configuration (refactored)
# Only complexity control retained. Execution stops if (d+1)^n exceeds threshold.
DP_CSP_DEFAULTS = {
    "max_d": None,  # int|null: upper bound distance to test (None -> baseline distance)
    "complexity_abort_threshold": 2000000000,  # int: abort if (d+1)^n exceeds this value
    "search_strategy": "linear",  # "linear" | "binary" (binary performs binary search over d)
}
