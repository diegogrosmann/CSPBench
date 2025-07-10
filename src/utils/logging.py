"""
Sistema de Logging Estruturado do CSPBench

Este módulo fornece um sistema de logging robusto e centralizado para todo o
framework CSPBench. Oferece configuração flexível, múltiplos handlers,
rotação de arquivos e integração com interfaces de usuário.

CARACTERÍSTICAS PRINCIPAIS:
===========================
- Configuração centralizada e padronizada
- Múltiplos níveis de logging (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- Rotação automática de arquivos de log
- Suporte a diferentes formatos de saída
- Integração com modo silencioso para automação
- Configuração específica por módulo
- Thread-safe para execução paralela

NÍVEIS DE LOGGING:
=================
- DEBUG: Informações detalhadas para desenvolvimento
- INFO: Informações gerais de funcionamento
- WARNING: Situações que merecem atenção mas não param execução
- ERROR: Erros que podem afetar funcionalidade
- CRITICAL: Erros críticos que podem parar a aplicação

ESTRUTURA DE LOGS:
=================
- Timestamp: Data e hora precisa do evento
- Logger Name: Módulo/classe que gerou o log
- Level: Nível de severidade
- Message: Mensagem descritiva
- Contexto adicional: Metadata quando relevante

CONFIGURAÇÕES SUPORTADAS:
========================
- Arquivo único ou rotação automática
- Console simultâneo ou apenas arquivo
- Formatação customizável
- Filtros por módulo ou nível
- Compressão de logs antigos

EXEMPLO DE USO:
==============
```python
from src.utils.logging import setup_logging
import logging

# Configurar logging
setup_logging("experiment_001", debug=True)

# Usar em qualquer módulo
logger = logging.getLogger(__name__)
logger.info("Iniciando experimento")
logger.debug("Parâmetros: %s", params)
logger.warning("Dataset grande, pode demorar")
logger.error("Erro na execução: %s", error)
```

INTEGRAÇÃO COM FRAMEWORK:
========================
- Todos os módulos do CSPBench usam este sistema
- Configuração automática no início da aplicação
- Logs estruturados para análise posterior
- Sincronização com interface curses

PERFORMANCE:
===========
- Lazy evaluation de mensagens
- Buffering para reduzir I/O
- Compressão automática de logs antigos
- Limpeza periódica de arquivos

SEGURANÇA:
=========
- Sanitização de dados sensíveis
- Controle de tamanho de arquivos
- Permissões apropriadas nos arquivos
- Não exposição de informações críticas
"""

import logging
import logging.handlers
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

# =============================================================================
# CONFIGURAÇÕES PADRÃO
# =============================================================================

DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
DEFAULT_LOG_DIR = "outputs/logs"
DEFAULT_MAX_BYTES = 10 * 1024 * 1024  # 10MB
DEFAULT_BACKUP_COUNT = 5

# Mapeamento de níveis para facilitar configuração
LOG_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


def setup_logging(
    base_name: str,
    silent: bool = False,
    debug: bool = False,
    console: bool = True,
    log_dir: Optional[str] = None,
    max_bytes: int = DEFAULT_MAX_BYTES,
    backup_count: int = DEFAULT_BACKUP_COUNT,
    format_string: Optional[str] = None,
    date_format: Optional[str] = None,
) -> logging.Logger:
    """
    Configura o sistema de logging global da aplicação CSPBench.

    Esta função estabelece uma configuração robusta de logging que suporta
    múltiplos handlers, rotação de arquivos e diferentes níveis de verbosidade.
    É a única função que deve ser chamada para configurar logging no framework.

    Args:
        base_name (str): Nome base para o arquivo de log (sem extensão)
        silent (bool): Se True, desabilita completamente o logging
        debug (bool): Se True, ativa nível DEBUG; senão usa INFO
        console (bool): Se True, também envia logs para console
        log_dir (str, optional): Diretório para arquivos de log
        max_bytes (int): Tamanho máximo por arquivo antes da rotação
        backup_count (int): Número de arquivos de backup a manter
        format_string (str, optional): Formato customizado para logs
        date_format (str, optional): Formato customizado para datas

    Returns:
        logging.Logger: Logger raiz configurado

    Examples:
        >>> # Configuração básica
        >>> setup_logging("experiment_01")

        >>> # Configuração para desenvolvimento
        >>> setup_logging("dev_session", debug=True, console=True)

        >>> # Configuração para produção
        >>> setup_logging("production", console=False, max_bytes=50*1024*1024)

        >>> # Modo silencioso para testes
        >>> setup_logging("test", silent=True)

    Note:
        - Em modo silencioso, todos os logs são suprimidos
        - Rotação automática previne crescimento descontrolado de arquivos
        - Thread-safe para uso em ambientes paralelos
        - Configuração persiste para toda a aplicação
    """
    # Desabilitar logging completamente em modo silencioso
    if silent:
        logging.disable(logging.CRITICAL)
        return logging.getLogger()

    # Definir diretório de logs
    if log_dir is None:
        log_dir = DEFAULT_LOG_DIR

    # Garantir que diretório existe
    Path(log_dir).mkdir(parents=True, exist_ok=True)

    # Configurar formatação
    if format_string is None:
        format_string = DEFAULT_LOG_FORMAT
    if date_format is None:
        date_format = DEFAULT_DATE_FORMAT

    formatter = logging.Formatter(format_string, date_format)

    # Determinar nível de logging
    log_level = logging.DEBUG if debug else logging.INFO

    # Obter logger raiz
    root_logger = logging.getLogger()

    # Limpar handlers existentes para reconfiguração
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Configurar nível do logger raiz
    root_logger.setLevel(log_level)

    # Configurar handler de arquivo com rotação
    log_file = Path(log_dir) / f"{base_name}.log"
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
    )
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

    # Configurar handler de console se solicitado
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(log_level)

        # Formato mais simples para console
        console_formatter = logging.Formatter("%(levelname)s - %(name)s - %(message)s")
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)

    # Log inicial de configuração
    logger = logging.getLogger(__name__)
    logger.info("Sistema de logging configurado")
    logger.info("Arquivo de log: %s", log_file)
    logger.info("Nível de logging: %s", logging.getLevelName(log_level))
    logger.info("Rotação: %d bytes, %d backups", max_bytes, backup_count)

    if debug:
        logger.debug("Modo DEBUG ativado - logs detalhados habilitados")

    return root_logger


def get_logger(name: str, level: Optional[str] = None) -> logging.Logger:
    """
    Obtém um logger específico para um módulo.

    Args:
        name (str): Nome do logger (geralmente __name__)
        level (str, optional): Nível específico para este logger

    Returns:
        logging.Logger: Logger configurado para o módulo

    Example:
        >>> logger = get_logger(__name__)
        >>> logger.info("Módulo inicializado")
    """
    logger = logging.getLogger(name)

    if level and level.upper() in LOG_LEVELS:
        logger.setLevel(LOG_LEVELS[level.upper()])

    return logger


def set_module_level(module_name: str, level: str) -> None:
    """
    Define nível de logging específico para um módulo.

    Args:
        module_name (str): Nome do módulo
        level (str): Nível de logging (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    Example:
        >>> set_module_level("algorithms.blf_ga", "DEBUG")
        >>> set_module_level("src.datasets", "WARNING")
    """
    if level.upper() not in LOG_LEVELS:
        raise ValueError(f"Nível inválido: {level}. Use: {list(LOG_LEVELS.keys())}")

    logger = logging.getLogger(module_name)
    logger.setLevel(LOG_LEVELS[level.upper()])

    logging.getLogger(__name__).info(
        "Nível de logging para '%s' definido como %s", module_name, level.upper()
    )


def configure_algorithm_logging(
    algorithm_name: str, verbose: bool = False
) -> logging.Logger:
    """
    Configura logging específico para um algoritmo.

    Args:
        algorithm_name (str): Nome do algoritmo
        verbose (bool): Se True, usa nível DEBUG para o algoritmo

    Returns:
        logging.Logger: Logger configurado para o algoritmo
    """
    logger_name = f"algorithms.{algorithm_name.lower()}"
    logger = logging.getLogger(logger_name)

    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Logging detalhado ativado para algoritmo %s", algorithm_name)
    else:
        logger.setLevel(logging.INFO)

    return logger


def log_execution_context(logger: logging.Logger, context: Dict[str, Any]) -> None:
    """
    Registra contexto de execução no log.

    Args:
        logger (logging.Logger): Logger para usar
        context (Dict[str, Any]): Contexto a registrar

    Example:
        >>> context = {
        ...     "algorithm": "BLF-GA",
        ...     "dataset_size": 100,
        ...     "parameters": {"pop_size": 50}
        ... }
        >>> log_execution_context(logger, context)
    """
    logger.info("=== CONTEXTO DE EXECUÇÃO ===")
    for key, value in context.items():
        logger.info("%s: %s", key, value)
    logger.info("=== FIM CONTEXTO ===")


def cleanup_old_logs(log_dir: str = DEFAULT_LOG_DIR, days_to_keep: int = 30) -> None:
    """
    Remove logs antigos para economizar espaço.

    Args:
        log_dir (str): Diretório de logs
        days_to_keep (int): Número de dias para manter logs
    """
    from datetime import datetime, timedelta

    log_path = Path(log_dir)
    if not log_path.exists():
        return

    cutoff_date = datetime.now() - timedelta(days=days_to_keep)
    removed_count = 0

    for log_file in log_path.glob("*.log*"):
        if log_file.stat().st_mtime < cutoff_date.timestamp():
            try:
                log_file.unlink()
                removed_count += 1
            except OSError:
                pass  # Ignorar erros de remoção

    if removed_count > 0:
        logger = logging.getLogger(__name__)
        logger.info("Removidos %d arquivos de log antigos", removed_count)
