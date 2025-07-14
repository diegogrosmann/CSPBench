"""
Configuração de Logging para CSPBench

Configurações centralizadas de logging sem dependências de console externas.
"""

import logging
import logging.handlers
from pathlib import Path
from typing import Optional


def setup_basic_logging(
    level: str = "INFO",
    log_dir: str = "outputs/logs",
    base_name: str = "cspbench",
    max_bytes: int = 10 * 1024 * 1024,  # 10MB
    backup_count: int = 5,
    log_file_path: Optional[str] = None,  # Caminho completo do arquivo se especificado
) -> None:
    """
    Configura logging básico do sistema.

    Args:
        level: Nível de log (DEBUG, INFO, WARNING, ERROR)
        log_dir: Diretório para arquivos de log
        base_name: Nome base para arquivos de log
        max_bytes: Tamanho máximo do arquivo antes da rotação
        backup_count: Número de backups a manter
        log_file_path: Caminho completo do arquivo (sobrepõe log_dir + base_name)
    """
    # Determinar caminho do arquivo
    if log_file_path:
        log_file = Path(log_file_path)
        log_file.parent.mkdir(parents=True, exist_ok=True)
    else:
        # Criar diretório se não existir
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        log_file = log_path / f"{base_name}.log"

    # Configurar formatação
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Configurar handler para arquivo com rotação
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
    )
    file_handler.setFormatter(formatter)

    # Configurar logger raiz
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Limpar handlers existentes
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Adicionar novo handler
    root_logger.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    """
    Retorna logger configurado para módulo específico.

    Args:
        name: Nome do módulo/logger

    Returns:
        logging.Logger: Logger configurado
    """
    return logging.getLogger(name)


class LoggerConfig:
    """Configuração centralizada de loggers."""

    _initialized = False
    _log_level = "INFO"

    @classmethod
    def initialize(cls, level: str = "INFO", **kwargs) -> None:
        """Inicializa configuração de logging."""
        if not cls._initialized:
            setup_basic_logging(level=level, **kwargs)
            cls._log_level = level
            cls._initialized = True

    @classmethod
    def get_level(cls) -> str:
        """Retorna nível atual de logging."""
        return cls._log_level

    @classmethod
    def set_level(cls, level: str) -> None:
        """Altera nível de logging."""
        cls._log_level = level
        logging.getLogger().setLevel(getattr(logging, level.upper(), logging.INFO))

    @classmethod
    def is_initialized(cls) -> bool:
        """Verifica se logging foi inicializado."""
        return cls._initialized
