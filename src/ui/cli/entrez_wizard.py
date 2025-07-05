"""
Wizard para entrada de parâmetros do NCBI Entrez.

Este módulo fornece uma interface para coletar parâmetros
necessários para buscar datasets no NCBI.
"""

import logging
from typing import Any

from src.utils.config import ENTREZ_DEFAULTS, safe_input

logger = logging.getLogger(__name__)


def collect_entrez_parameters() -> dict[str, Any]:
    """
    Coleta parâmetros para busca no NCBI via interface interativa.

    Returns:
        Dicionário com parâmetros para busca no NCBI
    """
    defaults = ENTREZ_DEFAULTS

    # Coletar email
    email_input = safe_input(f"Informe seu e-mail para o NCBI [{defaults['email']}]: ")
    email = email_input or defaults["email"]

    if not email:
        raise ValueError("É necessário fornecer um e-mail para acessar o NCBI")

    # Coletar API key (opcional)
    api_key_input = safe_input(f"Informe sua NCBI API key (opcional) [{defaults.get('api_key', '')}]: ")
    api_key = api_key_input or defaults.get("api_key", "")

    # Coletar base de dados
    db_input = safe_input(f"Base (nucleotide / protein) [{defaults['db']}]: ")
    db = db_input or defaults["db"]

    # Coletar termo de busca
    term_input = safe_input(f"Termo de busca Entrez [{defaults['term']}]: ")
    term = term_input or defaults["term"]

    # Coletar número de registros
    n_input = safe_input(f"Quantos registros deseja baixar? [{defaults['n']}]: ")
    n = int(n_input) if n_input else defaults["n"]

    params = {
        "email": email,
        "db": db,
        "term": term,
        "n": n,
    }

    if api_key:
        params["api_key"] = api_key

    logger.debug(f"Parâmetros coletados: {params}")
    return params
