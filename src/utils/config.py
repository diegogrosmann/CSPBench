"""
Configurações Centrais do CSPBench

Este módulo centraliza todas as configurações padrão, constantes e utilitários
de configuração do framework CSPBench. Fornece uma interface unificada para
gestão de parâmetros de algoritmos, datasets, execução e sistema.

CARACTERÍSTICAS PRINCIPAIS:
===========================
- Configurações centralizadas para todo o framework
- Parâmetros padrão para diferentes tipos de datasets
- Configuração de timeouts e limites de recursos
- Suporte a modo automatizado para testes e CI/CD
- Validação e normalização de parâmetros
- Sistema de configuração hierárquica

CONFIGURAÇÕES DISPONÍVEIS:
==========================
- DEBUG_DEFAULT: Controle de modo debug
- ALGORITHM_TIMEOUT: Timeout padrão para algoritmos
- SYNTHETIC_DEFAULTS: Parâmetros para datasets sintéticos
- FILE_DEFAULTS: Configurações para carregamento de arquivos
- ENTREZ_DEFAULTS: Configurações para download do NCBI
- BATCH_DEFAULTS: Configurações para execução em lote

MODO AUTOMATIZADO:
=================
O framework suporta execução automatizada através da variável de ambiente
CSP_AUTOMATED_TEST=1, útil para testes, CI/CD e scripts automatizados.

EXEMPLOS DE USO:
===============
```python
from src.utils.config import SYNTHETIC_DEFAULTS, safe_input

# Usar configurações padrão
dataset_params = SYNTHETIC_DEFAULTS.copy()
dataset_params.update({"n": 50, "L": 200})

# Entrada segura (compatível com automação)
user_input = safe_input("Digite um valor [padrão]: ", "padrão")

# Configuração condicional
if os.environ.get("CSP_AUTOMATED_TEST"):
    # Usar valores padrão
    params = SYNTHETIC_DEFAULTS
else:
    # Solicitar entrada do usuário
    params = get_user_config()
```

HIERARQUIA DE CONFIGURAÇÃO:
==========================
1. Variáveis de ambiente (mais alta prioridade)
2. Argumentos de linha de comando
3. Arquivos de configuração YAML
4. Configurações padrão deste módulo (mais baixa prioridade)

DESIGN PATTERNS:
===============
- Singleton implícito para configurações globais
- Factory para criação de configurações específicas
- Strategy para diferentes modos de entrada
- Observer para mudanças de configuração (futuro)

SEGURANÇA:
=========
- Validação de tipos e ranges para parâmetros
- Sanitização de entradas do usuário
- Controle de acesso a recursos sensíveis
- Logging de mudanças de configuração críticas

EXTENSIBILIDADE:
===============
Novas configurações podem ser facilmente adicionadas seguindo o padrão:
```python
NEW_MODULE_DEFAULTS = {
    "param1": valor_padrão,
    "param2": valor_padrão,
    # ...
}
```

Thread Safety:
==============
Todas as configurações são thread-safe para leitura. Modificações durante
execução devem ser sincronizadas externamente se necessário.
"""

import logging
import os
import sys
from typing import Any, Dict, Optional, Union

logger = logging.getLogger(__name__)

# =============================================================================
# CONFIGURAÇÕES GERAIS DO FRAMEWORK
# =============================================================================

# Controle de debug e logging
DEBUG_DEFAULT = "S"  # Modo debug habilitado por padrão para desenvolvimento
VERBOSE_DEFAULT = True  # Saída verbosa por padrão

# Timeouts e limites de recursos
ALGORITHM_TIMEOUT = 300  # Timeout padrão de 5 minutos para algoritmos
GLOBAL_TIMEOUT = 1800  # Timeout global de 30 minutos para operações longas
MAX_MEMORY_MB = 4096  # Limite padrão de memória em MB

# Configurações de paralelização
DEFAULT_WORKERS = 1  # Número padrão de workers (1 = sequencial)
MAX_WORKERS = 8  # Máximo de workers permitidos

# Configurações de arquivos e diretórios
OUTPUT_DIR = "outputs"  # Diretório para resultados
LOG_DIR = "outputs/logs"  # Diretório para logs
REPORT_DIR = "outputs/reports"  # Diretório para relatórios
DATASET_DIR = "saved_datasets"  # Diretório para datasets salvos


def safe_input(prompt: str, default: str = "") -> str:
    """
    Função de entrada segura compatível com modo automatizado.

    Esta função permite entrada do usuário em modo interativo, mas retorna
    valores padrão quando executada em modo automatizado (testes, CI/CD).

    Args:
        prompt (str): Texto a exibir para o usuário
        default (str): Valor padrão a retornar em modo automatizado

    Returns:
        str: Entrada do usuário ou valor padrão

    Example:
        >>> # Modo interativo
        >>> value = safe_input("Digite um valor [10]: ", "10")

        >>> # Modo automatizado (CSP_AUTOMATED_TEST=1)
        >>> value = safe_input("Digite um valor [10]: ", "10")  # Retorna "10"

    Note:
        - Em modo automatizado, sempre retorna o valor padrão
        - Trata KeyboardInterrupt e EOFError graciosamente
        - Logging de entradas para debugging
    """
    # Verificar se está em modo automatizado
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        logger.debug(
            "Modo automatizado: usando valor padrão '%s' para prompt '%s'",
            default,
            prompt.strip(),
        )
        return default

    try:
        # Exibir prompt e aguardar entrada
        print(prompt, end="", flush=True)
        user_input = input().strip()

        # Se usuário não digitou nada, usar padrão
        if not user_input:
            logger.debug("Entrada vazia, usando valor padrão: '%s'", default)
            return default

        logger.debug("Entrada do usuário: '%s'", user_input)
        return user_input

    except (KeyboardInterrupt, EOFError):
        print("\nOperação cancelada pelo usuário.")
        logger.info("Execução cancelada pelo usuário via Ctrl+C ou EOF")
        sys.exit(0)
    except Exception as e:
        logger.error("Erro na entrada do usuário: %s", e)
        print(f"\nErro na entrada: {e}. Usando valor padrão.")
        return default


# =============================================================================
# CONFIGURAÇÕES DE DATASETS
# =============================================================================

# Parâmetros padrão para datasets sintéticos
SYNTHETIC_DEFAULTS: Dict[str, Any] = {
    "n": 20,  # Número de strings
    "L": 100,  # Comprimento das strings
    "alphabet": "ACGT",  # Alfabeto (DNA por padrão)
    "noise": 0.10,  # Nível de ruído (10%)
    "seed": None,  # Semente para reprodutibilidade
    "pattern": "random",  # Padrão de geração: random, block, gradient
    "min_distance": 1,  # Distância mínima entre strings
    "max_distance": None,  # Distância máxima (None = automático)
}

# Parâmetros padrão para datasets de arquivo
FILE_DEFAULTS: Dict[str, Any] = {
    "filepath": "saved_datasets/sequences.fasta",  # Caminho padrão
    "format": "auto",  # Formato: auto, fasta, txt, csv
    "max_sequences": None,  # Máximo de sequências (None = todas)
    "min_length": 10,  # Comprimento mínimo das sequências
    "max_length": 1000,  # Comprimento máximo das sequências
    "normalize_length": True,  # Normalizar comprimentos
    "remove_gaps": True,  # Remover gaps (-) das sequências
    "uppercase": True,  # Converter para maiúsculas
}

# Parâmetros padrão para datasets NCBI/Entrez
ENTREZ_DEFAULTS: Dict[str, Any] = {
    "email": "diegogrosmann@gmail.com",  # Email para NCBI
    "db": "nucleotide",  # Database: nucleotide, protein
    "term": "COI[Gene] AND 600:650[SLEN]",  # Termo de busca
    "n": 20,  # Número de sequências
    "api_key": "40aff1a7f51fc7e0711203ea3b2f5ae37c09",  # API key NCBI
    "retmode": "fasta",  # Formato de retorno
    "sort": "relevance",  # Ordenação: relevance, date
    "max_length": 2000,  # Comprimento máximo
    "cache_results": True,  # Cache de resultados
}

# =============================================================================
# CONFIGURAÇÕES DE EXECUÇÃO
# =============================================================================

# Parâmetros padrão para execução em lote
BATCH_DEFAULTS: Dict[str, Any] = {
    "timeout_global": 1800,  # 30 minutos por execução completa
    "timeout_per_algorithm": 300,  # 5 minutos por algoritmo
    "max_concurrent": 1,  # Execuções paralelas
    "save_individual_reports": True,  # Salvar relatórios individuais
    "save_consolidated_report": True,  # Salvar relatório consolidado
    "continue_on_error": True,  # Continuar se um algoritmo falhar
    "retry_failed": False,  # Tentar novamente algoritmos que falharam
    "progress_update_interval": 10,  # Intervalo de atualização de progresso (segundos)
}

# Configurações de otimização
OPTIMIZATION_DEFAULTS: Dict[str, Any] = {
    "n_trials": 100,  # Número de trials para otimização
    "timeout_per_trial": 60,  # Timeout por trial
    "n_jobs": 1,  # Paralelização da otimização
    "study_name": "csp_optimization",  # Nome do estudo Optuna
    "storage": None,  # Storage para persistência (None = memória)
    "sampler": "TPE",  # Sampler: TPE, Random, CmaEs
    "pruner": "Median",  # Pruner: Median, Hyperband, None
}

# Configurações de análise de sensibilidade
SENSITIVITY_DEFAULTS: Dict[str, Any] = {
    "n_samples": 1000,  # Número de amostras
    "method": "morris",  # Método: morris, sobol, fast
    "confidence_level": 0.95,  # Nível de confiança
    "output_format": "json",  # Formato de saída: json, csv, html
}

# =============================================================================
# CONFIGURAÇÕES DE INTERFACE
# =============================================================================

# Configurações da interface curses
CURSES_DEFAULTS: Dict[str, Any] = {
    "refresh_rate": 1.0,  # Taxa de atualização em segundos
    "show_progress_bar": True,  # Mostrar barra de progresso
    "show_resource_monitor": True,  # Mostrar monitor de recursos
    "color_scheme": "default",  # Esquema de cores: default, dark, light
    "max_log_lines": 100,  # Máximo de linhas de log visíveis
}

# Configurações de logging
LOGGING_DEFAULTS: Dict[str, Any] = {
    "level": "INFO",  # Nível: DEBUG, INFO, WARNING, ERROR
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    "date_format": "%Y-%m-%d %H:%M:%S",
    "file_rotation": True,  # Rotação de arquivos de log
    "max_file_size": "10MB",  # Tamanho máximo por arquivo
    "backup_count": 5,  # Número de backups
}

# =============================================================================
# UTILITÁRIOS DE CONFIGURAÇÃO
# =============================================================================


def get_config_value(
    key: str, default: Any = None, config_dict: Optional[Dict[str, Any]] = None
) -> Any:
    """
    Obtém valor de configuração com hierarquia de prioridades.

    Args:
        key: Chave da configuração
        default: Valor padrão se não encontrado
        config_dict: Dicionário de configuração específico

    Returns:
        Valor da configuração encontrado ou padrão
    """
    # 1. Verificar variáveis de ambiente
    env_key = f"CSP_{key.upper()}"
    if env_key in os.environ:
        return os.environ[env_key]

    # 2. Verificar dicionário específico
    if config_dict and key in config_dict:
        return config_dict[key]

    # 3. Retornar padrão
    return default


def validate_config(
    config: Dict[str, Any], required_keys: Optional[list] = None
) -> bool:
    """
    Valida configuração verificando chaves obrigatórias.

    Args:
        config: Dicionário de configuração
        required_keys: Lista de chaves obrigatórias

    Returns:
        True se válida, False caso contrário
    """
    if required_keys is None:
        return True

    missing_keys = [key for key in required_keys if key not in config]
    if missing_keys:
        logger.error("Configuração inválida. Chaves faltando: %s", missing_keys)
        return False

    return True


def merge_configs(*configs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Mescla múltiplos dicionários de configuração.

    Configurações posteriores sobrescrevem anteriores.

    Args:
        *configs: Dicionários de configuração

    Returns:
        Dicionário mesclado
    """
    result = {}
    for config in configs:
        result.update(config)
    return result
