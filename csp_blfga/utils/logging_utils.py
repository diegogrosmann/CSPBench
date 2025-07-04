# logging_utils.py
"""
Utilitários para configuração de logging da aplicação CSP.

Funções:
    setup_logging(debug_mode, log_file): Configura logging global.
"""
import logging


def setup_logging(base_name: str, silent: bool = False, debug_mode: bool = False):
    """
    Configura o logging global da aplicação.

    Args:
        base_name (str): Nome base para o arquivo de log.
        silent (bool): Se True, suprime mensagens de log.
        debug_mode (bool): Se True, ativa modo debug detalhado.
    """
    if not silent:
        log_file = f"logs/{base_name}.log"
        if debug_mode:
            logging.basicConfig(
                level=logging.DEBUG,
                format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                filename=log_file,
                filemode="w",
            )
            print(f"Debug mode enabled. Log saved to {log_file}")
        else:
            # Configurar logging básico para arquivo
            logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                filename=log_file,
                filemode="w",
            )
    else:
        logging.disable(logging.CRITICAL)
