"""
Módulo de otimização e análise de sensibilidade para CSP-BLFGA.

Este módulo integra o Optuna para otimização de hiperparâmetros e o SALib
para análise de sensibilidade, fornecendo ferramentas para melhorar o
desempenho dos algoritmos CSP.

Módulos:
    optuna_optimizer: Otimização de hiperparâmetros usando Optuna
    sensitivity_analyzer: Análise de sensibilidade usando SALib
    visualization: Visualizações para resultados de otimização e sensibilidade

Exemplos de uso:

    # Otimização de hiperparâmetros
    from src.optimization import optimize_algorithm

    result = optimize_algorithm(
        algorithm_name="BLF-GA",
        sequences=sequences,
        alphabet=alphabet,
        n_trials=100
    )

    # Análise de sensibilidade
    from src.optimization import analyze_algorithm_sensitivity

    result = analyze_algorithm_sensitivity(
        algorithm_name="BLF-GA",
        sequences=sequences,
        alphabet=alphabet,
        n_samples=1000
    )
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
