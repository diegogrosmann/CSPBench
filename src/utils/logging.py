# logging.py
"""
Sistema de logging padronizado para CSP-BLFGA.

Este módulo fornece configuração centralizada de logging para toda a aplicação,
com suporte a diferentes níveis de verbosidade e saídas.

Funções:
    setup_logging(base_name, silent, debug): Configura o sistema de logging global.
"""
import logging
import os


def setup_logging(base_name: str, silent: bool = False, debug: bool = False) -> None:
    """
    Configura o sistema de logging global da aplicação.

    Args:
        base_name (str): Nome base para o arquivo de log (sem extensão).
        silent (bool, optional): Se True, desabilita logging. Defaults to False.
        debug (bool, optional): Se True, ativa logging detalhado. Defaults to False.

    Note:
        - Em modo silent, o logging é completamente desabilitado
        - Em modo debug, logs são salvos com nível DEBUG
        - Em modo normal, logs são salvos com nível INFO
        - Arquivos de log são salvos no diretório 'outputs/logs/'
    """
    if silent:
        logging.disable(logging.CRITICAL)
        return

    # Garantir que o diretório logs existe
    os.makedirs("outputs/logs", exist_ok=True)

    log_file = f"outputs/logs/{base_name}.log"
    log_level = logging.DEBUG if debug else logging.INFO

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        filename=log_file,
        filemode="w",
        force=True,  # Reconfigurar logging se já foi configurado
    )

    if debug:
        print(f"Debug logging enabled. Logs saved to: {log_file}")
    else:
        print(f"Logging enabled. Logs saved to: {log_file}")
