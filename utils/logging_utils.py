# logging_utils.py
"""
Utilitários para configuração de logging da aplicação CSP.

Funções:
    setup_logging(debug_mode, log_file): Configura logging global.
"""
import logging


def setup_logging(debug_mode: bool = False, log_file: str = "debug.log"):
    """
    Configura o logging global da aplicação.

    Args:
        debug_mode (bool): Se True, ativa modo debug detalhado.
        log_file (str): Caminho do arquivo de log.
    """
    if debug_mode:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            filename=log_file,
            filemode="w",
        )
        print(f"Debug mode enabled. Log saved to {log_file}")
    else:
        logging.disable(logging.CRITICAL)
        print("Debug mode disabled.")
