"""
Utilitário central para cálculo de workers para paralelização.

Este módulo fornece funções centralizadas para calcular o número apropriado
de workers para diferentes fases de execução (Optuna, SALib, etc.), evitando
oversubscription e garantindo uso eficiente dos recursos.

Funções:
    get_cpu_count(): Obtém número de CPUs disponíveis
    calculate_optuna_workers(): Calcula workers para Optuna
    calculate_salib_workers(): Calcula workers para SALib
    calculate_internal_workers(): Calcula workers internos para algoritmos
    get_worker_config(): Obtém configuração completa de workers
"""

import logging
import multiprocessing
import os
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)


def get_cpu_count() -> int:
    """
    Obtém o número de CPUs disponíveis no sistema.

    Returns:
        int: Número de CPUs disponíveis
    """
    try:
        # Tentar obter do limite de cgroup primeiro (containers)
        if os.path.exists("/sys/fs/cgroup/cpu/cpu.cfs_quota_us"):
            with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us", "r") as f:
                quota = int(f.read().strip())
            with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us", "r") as f:
                period = int(f.read().strip())
            if quota > 0 and period > 0:
                return max(1, quota // period)

        # Fallback para multiprocessing
        return multiprocessing.cpu_count()
    except (OSError, ValueError):
        logger.warning("Não foi possível determinar número de CPUs, usando 1")
        return 1


def calculate_optuna_workers(
    total_cpus: Optional[int] = None,
    yaml_config: Optional[Dict] = None,
    algorithm_name: str = "BLF-GA",
) -> Tuple[int, int]:
    """
    Calcula número de workers para Optuna e workers internos para algoritmos.

    Args:
        total_cpus: Número total de CPUs (auto-detect se None)
        yaml_config: Configuração do YAML com seção parallel
        algorithm_name: Nome do algoritmo para ajustar workers internos

    Returns:
        Tuple[int, int]: (optuna_workers, internal_workers)
    """
    if total_cpus is None:
        total_cpus = get_cpu_count()

    # Obter configuração do YAML se disponível
    yaml_optuna_workers = None
    if yaml_config and "optimization_config" in yaml_config:
        opt_config = yaml_config["optimization_config"]
        if "parallel" in opt_config:
            yaml_optuna_workers = opt_config["parallel"].get("n_jobs")

    # Usar configuração do YAML ou calcular dinamicamente
    if yaml_optuna_workers is not None:
        if yaml_optuna_workers == -1:
            # -1 significa usar todos os CPUs
            optuna_workers = total_cpus
        else:
            # Limitar ao número de CPUs disponíveis
            optuna_workers = min(yaml_optuna_workers, total_cpus)
    else:
        # Usar todos os CPUs disponíveis (máximo eficiente)
        optuna_workers = total_cpus

    # Calcular workers internos para algoritmos
    # Algoritmos com paralelismo interno (como BLF-GA) podem usar 2 workers por núcleo
    if algorithm_name in ["BLF-GA", "CSC"]:
        # Permitir até 2 workers internos por worker externo
        internal_workers = min(2, max(1, (total_cpus * 2) // optuna_workers))
    else:
        # Algoritmos determinísticos ou sem paralelismo interno
        internal_workers = 1

    logger.info(
        f"Calculado workers: Optuna={optuna_workers}, Internal={internal_workers}"
    )
    return optuna_workers, internal_workers


def calculate_salib_workers(
    total_cpus: Optional[int] = None,
    yaml_config: Optional[Dict] = None,
    n_samples: int = 1000,
) -> int:
    """
    Calcula número de workers para SALib.

    Args:
        total_cpus: Número total de CPUs (auto-detect se None)
        yaml_config: Configuração do YAML com seção parallel
        n_samples: Número de amostras para análise

    Returns:
        int: Número de workers para SALib
    """
    if total_cpus is None:
        total_cpus = get_cpu_count()

    # Obter configuração do YAML se disponível
    yaml_salib_workers = None
    if yaml_config and "sensitivity_config" in yaml_config:
        sens_config = yaml_config["sensitivity_config"]
        if "parallel" in sens_config:
            yaml_salib_workers = sens_config["parallel"].get("n_jobs")

    # Usar configuração do YAML ou calcular dinamicamente
    if yaml_salib_workers is not None:
        if yaml_salib_workers == -1:
            # -1 significa usar todos os CPUs
            salib_workers = total_cpus
        else:
            # Limitar ao número de CPUs disponíveis
            salib_workers = min(yaml_salib_workers, total_cpus)
    else:
        # Usar todos os CPUs disponíveis, mas limitar baseado no número de amostras
        # Mínimo 1 worker, máximo igual ao número de CPUs
        salib_workers = min(total_cpus, max(1, n_samples // 5))  # 5 amostras por worker

    logger.info(f"Calculado workers SALib: {salib_workers}")
    return salib_workers


def calculate_internal_workers(
    algorithm_name: str, available_cpus: int, external_parallelism: bool = False
) -> int:
    """
    Calcula workers internos para algoritmos específicos.

    Args:
        algorithm_name: Nome do algoritmo
        available_cpus: CPUs disponíveis para o algoritmo
        external_parallelism: Se há paralelismo externo (Optuna/SALib)

    Returns:
        int: Número de workers internos
    """
    # Algoritmos com paralelismo interno nativo
    internal_parallel_algorithms = ["BLF-GA", "CSC"]

    if algorithm_name not in internal_parallel_algorithms:
        return 1

    if external_parallelism:
        # Com paralelismo externo, permitir até 2 workers internos
        return min(2, max(1, available_cpus))
    else:
        # Sem paralelismo externo, usar até 2 workers internos por CPU
        return min(2, max(1, available_cpus))


def get_worker_config(
    yaml_config: Optional[Dict] = None,
    context: str = "optuna",  # "optuna", "salib", "algorithm"
    algorithm_name: str = "BLF-GA",
    n_samples: Optional[int] = None,
) -> Dict[str, int]:
    """
    Obtém configuração completa de workers baseada no contexto.

    Args:
        yaml_config: Configuração do YAML
        context: Contexto de execução
        algorithm_name: Nome do algoritmo
        n_samples: Número de amostras (para SALib)

    Returns:
        Dict[str, int]: Configuração de workers
    """
    total_cpus = get_cpu_count()

    if context == "optuna":
        optuna_workers, internal_workers = calculate_optuna_workers(
            total_cpus, yaml_config, algorithm_name
        )
        return {
            "optuna_workers": optuna_workers,
            "internal_workers": internal_workers,
            "total_cpus": total_cpus,
        }

    elif context == "salib":
        salib_workers = calculate_salib_workers(
            total_cpus, yaml_config, n_samples or 1000
        )
        return {"salib_workers": salib_workers, "total_cpus": total_cpus}

    elif context == "algorithm":
        internal_workers = calculate_internal_workers(
            algorithm_name, total_cpus, external_parallelism=False
        )
        return {"internal_workers": internal_workers, "total_cpus": total_cpus}

    else:
        raise ValueError(f"Contexto desconhecido: {context}")
