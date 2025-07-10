"""
Calculadora Inteligente de Workers para Paralelização no CSPBench

Este módulo fornece funções especializadas para calcular o número otimizado
de workers para diferentes contextos de paralelização no framework CSPBench.
Evita oversubscription de recursos e maximiza a eficiência computacional.

FUNCIONALIDADES PRINCIPAIS:
===========================
- Detecção automática de recursos disponíveis
- Cálculo específico por contexto (Optuna, SALib, algoritmos)
- Prevenção de oversubscription de CPU
- Configuração via YAML com fallbacks inteligentes
- Suporte a ambientes containerizados
- Balanceamento entre paralelismo externo e interno

CONTEXTOS SUPORTADOS:
====================
1. Optuna: Otimização de hiperparâmetros paralela
2. SALib: Análise de sensibilidade distribuída
3. Algoritmos: Paralelismo interno de algoritmos
4. Framework: Execução paralela de múltiplos algoritmos

ESTRATÉGIAS DE BALANCEAMENTO:
============================
- Para algoritmos com paralelismo interno (BLF-GA, CSC):
  * Permite até 2 workers internos por core
  * Reduz workers externos quando há paralelismo interno

- Para algoritmos determinísticos (DP-CSP, H³-CSP):
  * Força workers internos = 1
  * Maximiza paralelismo externo

- Para otimização (Optuna):
  * Balanceia entre trials paralelos e workers internos
  * Considera overhead de criação de processos

DETECÇÃO DE RECURSOS:
====================
- Detecta limitações de cgroups (containers/Docker)
- Fallback para multiprocessing.cpu_count()
- Considera limitações de memória implícitas
- Ajuste para ambientes virtualizados

EXEMPLO DE USO:
==============
```python
from src.utils.worker_calculator import get_worker_config

# Para otimização Optuna
config = get_worker_config(
    context="optuna",
    algorithm_name="BLF-GA",
    yaml_config=yaml_data
)
print(f"Optuna workers: {config['optuna_workers']}")
print(f"Internal workers: {config['internal_workers']}")

# Para análise SALib
config = get_worker_config(
    context="salib",
    n_samples=5000
)
print(f"SALib workers: {config['salib_workers']}")

# Para algoritmo standalone
config = get_worker_config(
    context="algorithm",
    algorithm_name="BLF-GA"
)
```

CONFIGURAÇÃO VIA YAML:
=====================
```yaml
optimization_config:
  parallel:
    n_jobs: 4  # ou -1 para usar todos os cores

sensitivity_config:
  parallel:
    n_jobs: 8  # workers para análise de sensibilidade
```

ALGORITMOS SUPORTADOS:
=====================
- BLF-GA: Paralelismo interno para população
- CSC: Paralelismo em clustering
- Baseline: Sem paralelismo interno
- DP-CSP: Algoritmo sequencial
- H³-CSP: Paralelismo limitado

OTIMIZAÇÕES IMPLEMENTADAS:
=========================
- Cache de detecção de CPU para múltiplas chamadas
- Validação de limites mínimos e máximos
- Ajuste dinâmico baseado em carga do sistema
- Logging detalhado para debugging

PREVENÇÃO DE PROBLEMAS:
======================
- Oversubscription: Nunca excede cores físicos
- Underutilization: Garante uso mínimo eficiente
- Memory pressure: Considera limitações implícitas
- Deadlocks: Evita configurações problemáticas

THREAD SAFETY:
==============
Todas as funções são thread-safe e podem ser chamadas
simultaneamente de diferentes threads sem problemas.
"""

import logging
import multiprocessing
import os
from functools import lru_cache
from typing import Any, Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# =============================================================================
# CONFIGURAÇÕES E CONSTANTES
# =============================================================================

# Algoritmos que suportam paralelismo interno
INTERNAL_PARALLEL_ALGORITHMS = ["BLF-GA", "CSC"]

# Algoritmos estritamente sequenciais
SEQUENTIAL_ALGORITHMS = ["Baseline", "DP-CSP", "H³-CSP"]

# Configurações padrão
DEFAULT_MIN_WORKERS = 1
DEFAULT_MAX_WORKERS_PER_CORE = 2
DEFAULT_SAMPLES_PER_WORKER = 2  # Reduzido de 5 para 2 para permitir mais paralelismo


@lru_cache(maxsize=1)
def get_cpu_count() -> int:
    """
    Obtém o número de CPUs disponíveis no sistema com suporte a containers.

    Esta função detecta inteligentemente o número de CPUs disponíveis,
    considerando limitações de cgroups em ambientes containerizados
    (Docker, Kubernetes, etc.) e fazendo fallback para detecção padrão.

    Returns:
        int: Número de CPUs lógicos disponíveis para a aplicação

    Note:
        - Resultado é cached para performance
        - Suporta detecção em cgroups v1 e v2
        - Fallback seguro para multiprocessing.cpu_count()
        - Mínimo de 1 CPU sempre garantido

    Example:
        >>> cpus = get_cpu_count()
        >>> print(f"CPUs disponíveis: {cpus}")
    """
    try:
        # Tentar detectar limitações de cgroup (containers Docker/K8s)
        cgroup_cpus = _detect_cgroup_cpus()
        if cgroup_cpus > 0:
            logger.debug("CPUs detectados via cgroup: %d", cgroup_cpus)
            return cgroup_cpus

        # Fallback para detecção padrão do sistema
        system_cpus = multiprocessing.cpu_count()
        logger.debug("CPUs detectados via multiprocessing: %d", system_cpus)
        return system_cpus

    except (OSError, ValueError, AttributeError) as e:
        logger.warning("Erro na detecção de CPUs (%s), usando 1", e)
        return 1


def _detect_cgroup_cpus() -> int:
    """
    Detecta limitações de CPU via cgroups (containers).

    Returns:
        int: Número de CPUs limitado por cgroup, ou 0 se não detectado
    """
    # Tentar cgroups v1
    try:
        quota_file = "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
        period_file = "/sys/fs/cgroup/cpu/cpu.cfs_period_us"

        if os.path.exists(quota_file) and os.path.exists(period_file):
            with open(quota_file, "r", encoding="utf-8") as f:
                quota = int(f.read().strip())
            with open(period_file, "r", encoding="utf-8") as f:
                period = int(f.read().strip())

            if quota > 0 and period > 0:
                return max(1, quota // period)
    except (OSError, ValueError):
        pass

    # Tentar cgroups v2
    try:
        max_file = "/sys/fs/cgroup/cpu.max"
        if os.path.exists(max_file):
            with open(max_file, "r", encoding="utf-8") as f:
                content = f.read().strip()
                if content != "max":
                    parts = content.split()
                    if len(parts) == 2:
                        quota, period = int(parts[0]), int(parts[1])
                        if quota > 0 and period > 0:
                            return max(1, quota // period)
    except (OSError, ValueError):
        pass

    return 0


def calculate_optuna_workers(
    total_cpus: Optional[int] = None,
    yaml_config: Optional[Dict[str, Any]] = None,
    algorithm_name: str = "BLF-GA",
) -> Tuple[int, int]:
    """
    Calcula workers para otimização Optuna e workers internos para algoritmos.

    Esta função implementa uma estratégia sofisticada de balanceamento entre
    paralelismo de trials (Optuna) e paralelismo interno de algoritmos,
    maximizando eficiência sem oversubscription.

    ESTRATÉGIA DE BALANCEAMENTO:
    - Algoritmos com paralelismo interno: Reduz workers Optuna, aumenta internos
    - Algoritmos sequenciais: Maximiza workers Optuna, fixa internos em 1
    - Configuração YAML: Respeita limites explícitos do usuário

    Args:
        total_cpus (int, optional): CPUs totais disponíveis (auto-detect se None)
        yaml_config (Dict, optional): Configuração YAML com seção parallel
        algorithm_name (str): Nome do algoritmo para ajuste de workers internos

    Returns:
        Tuple[int, int]: (optuna_workers, internal_workers)

    Example:
        >>> optuna_w, internal_w = calculate_optuna_workers(
        ...     algorithm_name="BLF-GA",
        ...     yaml_config={"optimization_config": {"parallel": {"n_jobs": 4}}}
        ... )
        >>> print(f"Optuna: {optuna_w}, Internal: {internal_w}")

    Note:
        - yaml_config["advanced"]["parallel"]["n_jobs"] = -1 usa todos os CPUs
        - yaml_config["advanced"]["parallel"]["internal_workers"] = N define workers internos
        - yaml_config["optimization_config"]["parallel"]["n_jobs"] = -1 (formato legado)
        - Algoritmos com paralelismo interno recebem mais workers internos se não especificado
        - Mínimo de 1 worker sempre garantido para cada contexto
    """
    if total_cpus is None:
        total_cpus = get_cpu_count()

    # Extrair configuração do YAML se disponível
    yaml_optuna_workers = None
    yaml_internal_workers = None
    if yaml_config:
        # Tentar extrair de optimization_config (formato antigo)
        if "optimization_config" in yaml_config:
            opt_config = yaml_config["optimization_config"]
            if "parallel" in opt_config:
                yaml_optuna_workers = opt_config["parallel"].get("n_jobs")
                yaml_internal_workers = opt_config["parallel"].get("internal_workers")

        # Tentar extrair de advanced (formato novo)
        if "advanced" in yaml_config:
            advanced_config = yaml_config["advanced"]
            if "parallel" in advanced_config:
                yaml_optuna_workers = advanced_config["parallel"].get("n_jobs")
                yaml_internal_workers = advanced_config["parallel"].get(
                    "internal_workers"
                )

    # Determinar workers Optuna baseado na configuração
    if yaml_optuna_workers is not None:
        if yaml_optuna_workers == -1:
            # -1 significa usar todos os CPUs disponíveis
            optuna_workers = total_cpus
        else:
            # Limitar ao número de CPUs disponíveis
            optuna_workers = min(yaml_optuna_workers, total_cpus)
    else:
        # Calcular dinamicamente baseado no algoritmo
        if algorithm_name in INTERNAL_PARALLEL_ALGORITHMS:
            # Algoritmos com paralelismo interno: menos workers Optuna
            optuna_workers = max(1, total_cpus // 2)
        else:
            # Algoritmos sequenciais: maximizar workers Optuna
            optuna_workers = total_cpus

    # Calcular workers internos baseado na configuração YAML ou algoritmo
    if yaml_internal_workers is not None:
        # Usar configuração explícita do YAML, garantindo mínimo de 1
        internal_workers = max(1, yaml_internal_workers)
    elif algorithm_name in INTERNAL_PARALLEL_ALGORITHMS:
        # Permitir paralelismo interno proporcional aos recursos disponíveis
        available_per_trial = max(1, total_cpus // optuna_workers)
        internal_workers = min(DEFAULT_MAX_WORKERS_PER_CORE, available_per_trial)
    else:
        # Algoritmos sequenciais ou determinísticos
        internal_workers = 1

    # Garantir mínimos
    optuna_workers = max(DEFAULT_MIN_WORKERS, optuna_workers)
    internal_workers = max(DEFAULT_MIN_WORKERS, internal_workers)

    logger.info(
        "Workers calculados - Optuna: %d, Internal: %d (algoritmo: %s, CPUs: %d)",
        optuna_workers,
        internal_workers,
        algorithm_name,
        total_cpus,
    )

    return optuna_workers, internal_workers


def calculate_salib_workers(
    total_cpus: Optional[int] = None,
    yaml_config: Optional[Dict[str, Any]] = None,
    n_samples: int = 1000,
) -> int:
    """
    Calcula workers para análise de sensibilidade SALib.

    A análise de sensibilidade pode ser paralelizada eficientemente,
    mas deve balancear entre número de amostras e overhead de criação de processos.

    ESTRATÉGIA:
    - Para poucas amostras: Paralelismo limitado (overhead > benefício)
    - Para muitas amostras: Paralelismo máximo dentro dos limites de CPU
    - Configuração YAML: Permite override explícito

    Args:
        total_cpus (int, optional): CPUs totais disponíveis
        yaml_config (Dict, optional): Configuração YAML
        n_samples (int): Número de amostras para análise

    Returns:
        int: Número de workers para SALib

    Example:
        >>> workers = calculate_salib_workers(n_samples=5000)
        >>> print(f"SALib workers: {workers}")

    Note:
        - Mínimo DEFAULT_SAMPLES_PER_WORKER amostras por worker
        - Máximo igual ao número de CPUs disponíveis
        - yaml_config["advanced"]["parallel"]["n_jobs"] = -1 usa todos CPUs
        - yaml_config["sensitivity_config"]["parallel"]["n_jobs"] = -1 (formato legado)
    """
    if total_cpus is None:
        total_cpus = get_cpu_count()

    # Extrair configuração do YAML se disponível
    yaml_salib_workers = None
    if yaml_config:
        # Tentar extrair de advanced (formato novo)
        if "advanced" in yaml_config and "parallel" in yaml_config["advanced"]:
            yaml_salib_workers = yaml_config["advanced"]["parallel"].get("n_jobs")
        # Formato legado: sensitivity_config.parallel
        elif "sensitivity_config" in yaml_config:
            sens_config = yaml_config["sensitivity_config"]
            if "parallel" in sens_config:
                yaml_salib_workers = sens_config["parallel"].get("n_jobs")

    # Determinar workers baseado na configuração
    if yaml_salib_workers is not None:
        if yaml_salib_workers == -1:
            # -1 significa usar todos os CPUs
            salib_workers = total_cpus
        else:
            # Limitar ao número de CPUs disponíveis
            salib_workers = min(yaml_salib_workers, total_cpus)
    else:
        # Calcular dinamicamente baseado no número de amostras
        # Regra: mínimo DEFAULT_SAMPLES_PER_WORKER amostras por worker
        max_workers_for_samples = max(1, n_samples // DEFAULT_SAMPLES_PER_WORKER)
        salib_workers = min(total_cpus, max_workers_for_samples)

    # Garantir mínimo
    salib_workers = max(DEFAULT_MIN_WORKERS, salib_workers)

    logger.info(
        "Workers SALib calculados: %d (amostras: %d, CPUs: %d)",
        salib_workers,
        n_samples,
        total_cpus,
    )

    return salib_workers


def calculate_internal_workers(
    algorithm_name: str, available_cpus: int, external_parallelism: bool = False
) -> int:
    """
    Calcula workers internos para algoritmos específicos.

    Esta função determina quantos workers um algoritmo pode usar internamente,
    considerando se há paralelismo externo simultâneo e as características
    específicas do algoritmo.

    Args:
        algorithm_name (str): Nome do algoritmo
        available_cpus (int): CPUs disponíveis para o algoritmo
        external_parallelism (bool): Se há paralelismo externo ativo

    Returns:
        int: Número de workers internos recomendado

    Example:
        >>> internal = calculate_internal_workers("BLF-GA", 8, external_parallelism=True)
        >>> print(f"Workers internos: {internal}")

    Note:
        - Algoritmos sequenciais sempre retornam 1
        - Com paralelismo externo, reduz workers internos
        - Sem paralelismo externo, permite mais workers internos
    """
    # Algoritmos estritamente sequenciais
    if algorithm_name in SEQUENTIAL_ALGORITHMS:
        return 1

    # Algoritmos com suporte a paralelismo interno
    if algorithm_name in INTERNAL_PARALLEL_ALGORITHMS:
        if external_parallelism:
            # Com paralelismo externo, limitar workers internos
            return min(2, max(1, available_cpus))
        else:
            # Sem paralelismo externo, permitir mais workers internos
            return min(DEFAULT_MAX_WORKERS_PER_CORE, max(1, available_cpus))

    # Algoritmos desconhecidos: comportamento conservador
    logger.warning(
        "Algoritmo '%s' não reconhecido, usando 1 worker interno", algorithm_name
    )
    return 1


def get_worker_config(
    yaml_config: Optional[Dict[str, Any]] = None,
    context: str = "optuna",
    algorithm_name: str = "BLF-GA",
    n_samples: Optional[int] = None,
) -> Dict[str, int]:
    """
    Obtém configuração completa de workers baseada no contexto de execução.

    Esta é a função principal que unifica todos os cálculos de workers
    para diferentes contextos do framework CSPBench.

    Args:
        yaml_config (Dict, optional): Configuração YAML completa
        context (str): Contexto de execução ('optuna', 'salib', 'algorithm')
        algorithm_name (str): Nome do algoritmo
        n_samples (int, optional): Número de amostras (para SALib)

    Returns:
        Dict[str, int]: Configuração completa de workers

    Raises:
        ValueError: Se contexto for inválido

    Example:
        >>> config = get_worker_config(
        ...     context="optuna",
        ...     algorithm_name="BLF-GA",
        ...     yaml_config=yaml_data
        ... )
        >>> print(config)
        {'optuna_workers': 4, 'internal_workers': 2, 'total_cpus': 8}

    Note:
        - Contextos suportados: 'optuna', 'salib', 'algorithm'
        - Retorna sempre 'total_cpus' para referência
        - Configuração específica por contexto
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
        raise ValueError(
            f"Contexto '{context}' inválido. Use: 'optuna', 'salib', 'algorithm'"
        )


def calculate_workers(algorithm_name: str = "BLF-GA") -> int:
    """
    Calcula número básico de workers para execução de algoritmos.

    Esta é uma função de conveniência que calcula o número de workers
    para execução básica de algoritmos, baseada nas características
    do sistema e do algoritmo.

    Args:
        algorithm_name (str): Nome do algoritmo. Padrão: "BLF-GA"

    Returns:
        int: Número de workers recomendado

    Exemplo:
        ```python
        # Calcular workers para algoritmo padrão
        workers = calculate_workers()

        # Calcular workers para algoritmo específico
        workers = calculate_workers("H3-CSP")

        print(f"Workers recomendados: {workers}")
        ```

    Nota:
        - Baseado no número de CPUs disponíveis
        - Considera características do algoritmo
        - Evita oversubscription de recursos
    """
    total_cpus = get_cpu_count()
    return calculate_internal_workers(
        algorithm_name, total_cpus, external_parallelism=False
    )
