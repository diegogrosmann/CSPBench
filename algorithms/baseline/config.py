"""
Default configurations for the baseline algorithm.

HIERARQUIA DE CONFIGURAÇÕES:
- Este arquivo define CONFIGURAÇÕES PADRÃO DE ALGORITMOS (prioridade mais baixa)
- Pode ser sobreposto por:
  * config/settings.yaml (configurações padrão de execução)
  * batches/*.yaml (configurações de execução específicas)
- Não pode sobrepor:
  * .env (configurações de sistema)

Attributes:
    BASELINE_DEFAULTS (dict): Default baseline parameters.
"""

BASELINE_DEFAULTS = {
    "tie_break": "lex",  # tie-breaking criterion: 'lex', 'random', 'first'
}
