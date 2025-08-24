"""
Default configurations for the baseline algorithm.

CONFIGURATION HIERARCHY:
- This file defines ALGORITHM DEFAULT CONFIG (lowest priority)
- Can be overridden by:
  * config/settings.yaml (execution defaults)
  * batches/*.yaml (batch specific overrides)
- Cannot override:
  * .env (system level settings)

Attributes:
    BASELINE_DEFAULTS (dict): Default baseline parameters.
"""

BASELINE_DEFAULTS = {
    "tie_break": "lex",  # tie-breaking criterion: 'lex', 'random', 'first'
}
