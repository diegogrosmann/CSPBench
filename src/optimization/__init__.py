"""
Módulo de Otimização e Análise - CSPBench

Este módulo fornece um sistema completo de otimização de hiperparâmetros e análise
de sensibilidade para algoritmos CSP, integrando ferramentas avançadas como Optuna
e SALib para maximizar a performance dos algoritmos.

Arquitetura:
    O módulo implementa uma arquitetura modular com:
    - Sistema de otimização baseado em Optuna
    - Análise de sensibilidade via SALib
    - Visualizações interativas avançadas
    - Configuração flexível de espaços de parâmetros
    - Execução paralela de otimizações

Funcionalidades Principais:
    - Otimização automática de hiperparâmetros
    - Análise de sensibilidade de parâmetros
    - Estudos comparativos entre algoritmos
    - Visualizações de espaços de parâmetros
    - Relatórios detalhados de resultados
    - Paralelização de experimentos

Módulos Incluídos:
    - optuna_optimizer: Otimização baseada em Optuna
    - sensitivity_analyzer: Análise de sensibilidade via SALib
    - batch_optimizer: Otimização em lote para múltiplos datasets
    - batch_sensitivity: Análise de sensibilidade em lote
    - visualization: Visualizações avançadas de resultados

Métodos de Otimização:
    - **TPE (Tree-structured Parzen Estimator)**: Padrão do Optuna
    - **Random Search**: Busca aleatória
    - **Grid Search**: Busca em grade
    - **CMA-ES**: Evolution Strategy
    - **Bayesian Optimization**: Otimização Bayesiana

Métodos de Análise de Sensibilidade:
    - **Morris**: Análise de screening inicial
    - **Sobol**: Índices de sensibilidade global
    - **FAST**: Fourier Amplitude Sensitivity Test
    - **Delta**: Análise baseada em momentos

Exemplo de Uso:
    ```python
    from src.optimization import (
        optimize_algorithm,
        analyze_algorithm_sensitivity,
        OptimizationVisualizer
    )

    # Otimização de hiperparâmetros
    result = optimize_algorithm(
        algorithm_name="BLF-GA",
        sequences=sequences,
        alphabet=alphabet,
        n_trials=100,
        timeout=300
    )

    # Análise de sensibilidade
    sensitivity = analyze_algorithm_sensitivity(
        algorithm_name="BLF-GA",
        sequences=sequences,
        alphabet=alphabet,
        n_samples=1000,
        method='sobol'
    )

    # Visualizar resultados
    viz = OptimizationVisualizer(result)
    viz.plot_optimization_history()
    viz.plot_parameter_importance()
    ```

Configuração de Espaços de Parâmetros:
    ```python
    param_space = {
        'population_size': ['int', 50, 300],
        'max_generations': ['int', 100, 500],
        'mutation_rate': ['uniform', 0.01, 0.3],
        'crossover_rate': ['uniform', 0.6, 0.9]
    }
    ```

Integração com Batch System:
    - Suporte a configurações YAML
    - Execução de múltiplos estudos
    - Comparação automática de resultados
    - Relatórios consolidados

Performance:
    - Execução paralela de trials
    - Cache de resultados intermediários
    - Pruning automático de trials ruins
    - Otimização de tempo de execução

Autor: CSPBench Development Team
Data: 2024
"""

from .optuna_optimizer import (
    OptimizationConfig,
    OptimizationResult,
    OptunaOptimizer,
    create_optimization_study,
    optimize_algorithm,
)
from .sensitivity_analyzer import (
    SensitivityAnalyzer,
    SensitivityConfig,
    SensitivityResult,
    analyze_algorithm_sensitivity,
    create_parameter_space,
)
from .visualization import OptimizationVisualizer, SensitivityVisualizer

__all__ = [
    # Otimização
    "optimize_algorithm",
    "OptunaOptimizer",
    "OptimizationConfig",
    "OptimizationResult",
    "create_optimization_study",
    # Análise de sensibilidade
    "analyze_algorithm_sensitivity",
    "SensitivityAnalyzer",
    "SensitivityConfig",
    "SensitivityResult",
    "create_parameter_space",
    # Visualizações
    "OptimizationVisualizer",
    "SensitivityVisualizer",
]
