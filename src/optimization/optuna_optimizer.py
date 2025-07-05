"""
Módulo de otimização de hiperparâmetros usando Optuna.

Classes:
    OptunaOptimizer: Otimizador principal usando Optuna
    OptimizationConfig: Configuração da otimização
    OptimizationResult: Resultado da otimização

Funções:
    optimize_algorithm: Otimiza hiperparâmetros de um algoritmo
    create_optimization_study: Cria estudo de otimização
"""

import logging
import time
from dataclasses import dataclass
from typing import Any

import optuna
from optuna.pruners import MedianPruner
from optuna.samplers import TPESampler

from algorithms.base import global_registry
from src.core.exec.algorithm_executor import AlgorithmExecutor
from src.ui.cli.console_manager import console

logger = logging.getLogger(__name__)


@dataclass
class OptimizationConfig:
    """Configuração para otimização de hiperparâmetros."""

    algorithm_name: str
    dataset_sequences: list[str]
    alphabet: str
    n_trials: int = 100
    timeout_per_trial: int = 300
    n_jobs: int = 1
    study_name: str | None = None
    direction: str = "minimize"  # minimize ou maximize
    sampler: str = "TPE"  # TPE, Random, CmaEs
    pruner: str = "median"  # median, hyperband, none

    def __post_init__(self):
        """Validação da configuração."""
        if self.algorithm_name not in global_registry:
            raise ValueError(f"Algoritmo não encontrado: {self.algorithm_name}")

        if self.direction not in ["minimize", "maximize"]:
            raise ValueError("direction deve ser 'minimize' ou 'maximize'")


@dataclass
class OptimizationResult:
    """Resultado da otimização."""

    best_params: dict[str, Any]
    best_value: float
    n_trials: int
    study_name: str
    optimization_time: float
    all_trials: list[dict]

    def summary(self) -> str:
        """Retorna resumo da otimização."""
        return (
            f"Otimização Concluída:\n"
            f"  Melhor valor: {self.best_value:.6f}\n"
            f"  Melhores parâmetros: {self.best_params}\n"
            f"  Trials executados: {self.n_trials}\n"
            f"  Tempo total: {self.optimization_time:.1f}s"
        )


class OptunaOptimizer:
    """Otimizador de hiperparâmetros usando Optuna."""

    def __init__(self, config: OptimizationConfig):
        self.config = config
        self.study: optuna.Study | None = None
        self.algorithm_class = global_registry[config.algorithm_name]

        # Configurar sampler
        if config.sampler == "TPE":
            sampler = TPESampler(seed=42)
        elif config.sampler == "Random":
            sampler = optuna.samplers.RandomSampler(seed=42)
        elif config.sampler == "CmaEs":
            sampler = optuna.samplers.CmaEsSampler(seed=42)
        else:
            sampler = TPESampler(seed=42)

        # Configurar pruner
        if config.pruner == "median":
            pruner = MedianPruner(n_startup_trials=5, n_warmup_steps=10)
        elif config.pruner == "hyperband":
            pruner = optuna.pruners.HyperbandPruner()
        else:
            pruner = optuna.pruners.NopPruner()

        # Criar estudo
        study_name = config.study_name or f"{config.algorithm_name}_optimization"
        self.study = optuna.create_study(
            direction=config.direction, sampler=sampler, pruner=pruner, study_name=study_name
        )

        logger.info(f"Criado estudo de otimização: {study_name}")

    def _suggest_parameters(self, trial: optuna.Trial) -> dict[str, Any]:
        """Define espaço de busca de hiperparâmetros baseado no algoritmo."""

        params = {}

        if self.config.algorithm_name == "BLF-GA":
            # Parâmetros específicos do BLF-GA
            params.update(
                {
                    "population_size": trial.suggest_int("population_size", 20, 200, step=10),
                    "max_generations": trial.suggest_int("max_generations", 100, 1000, step=50),
                    "crossover_rate": trial.suggest_float("crossover_rate", 0.6, 0.95, step=0.05),
                    "mutation_rate": trial.suggest_float("mutation_rate", 0.01, 0.3, step=0.01),
                    "elite_size": trial.suggest_int("elite_size", 1, 10),
                    "tournament_size": trial.suggest_int("tournament_size", 2, 8),
                }
            )

        elif self.config.algorithm_name == "CSC":
            # Parâmetros específicos do CSC
            params.update(
                {
                    "max_iterations": trial.suggest_int("max_iterations", 100, 2000, step=100),
                    "improvement_threshold": trial.suggest_float("improvement_threshold", 1e-6, 1e-3, log=True),
                    "restart_threshold": trial.suggest_int("restart_threshold", 10, 100, step=10),
                }
            )

        elif self.config.algorithm_name == "DP-CSP":
            # Parâmetros específicos do DP-CSP
            params.update(
                {
                    "beam_width": trial.suggest_int("beam_width", 10, 100, step=10),
                    "pruning_factor": trial.suggest_float("pruning_factor", 0.1, 0.9, step=0.1),
                }
            )

        elif self.config.algorithm_name == "H3-CSP":
            # Parâmetros específicos do H3-CSP
            params.update(
                {
                    "max_depth": trial.suggest_int("max_depth", 5, 20),
                    "branching_factor": trial.suggest_int("branching_factor", 2, 10),
                    "heuristic_weight": trial.suggest_float("heuristic_weight", 0.1, 2.0, step=0.1),
                }
            )

        # Parâmetros comuns se o algoritmo suportar
        common_params = {}
        if hasattr(self.algorithm_class, "supports_seed") and self.algorithm_class.supports_seed:
            common_params["seed"] = trial.suggest_int("seed", 1, 1000000)

        params.update(common_params)
        return params

    def _objective(self, trial: optuna.Trial) -> float:
        """Função objetivo para otimização."""

        try:
            # Sugerir parâmetros
            params = self._suggest_parameters(trial)

            # Criar instância do algoritmo com parâmetros
            algorithm_instance = self.algorithm_class(self.config.dataset_sequences, self.config.alphabet, **params)

            # Configurar timeout se suportado
            if hasattr(algorithm_instance, "set_timeout"):
                algorithm_instance.set_timeout(self.config.timeout_per_trial)

            # Executar algoritmo
            executor = AlgorithmExecutor(self.config.timeout_per_trial)
            result, _, _ = executor.execute_with_timeout(algorithm_instance)

            if result is None or "distancia" not in result:
                return float("inf")

            distance = result["distancia"]

            # Reportar valor intermediário para pruning
            trial.report(distance, step=0)

            # Verificar se deve ser podado
            if trial.should_prune():
                raise optuna.TrialPruned()

            return distance

        except optuna.TrialPruned:
            raise
        except Exception as e:
            logger.warning(f"Erro no trial {trial.number}: {e}")
            return float("inf")

    def optimize(self, show_progress: bool = True) -> OptimizationResult:
        """Executa otimização de hiperparâmetros."""

        if not self.study:
            raise ValueError("Study não foi inicializado corretamente")

        if show_progress:
            console.print(f"\n🔍 Iniciando otimização de {self.config.algorithm_name}")
            console.print(
                f"📊 Dataset: n={len(self.config.dataset_sequences)}, L={len(self.config.dataset_sequences[0])}"
            )
            console.print(f"🎯 Trials: {self.config.n_trials}")
            console.print(f"⏱️  Timeout por trial: {self.config.timeout_per_trial}s")

        start_time = time.time()

        # Callback para mostrar progresso
        def progress_callback(study: optuna.Study, trial):
            if show_progress and trial.number % 10 == 0:
                best_value = study.best_value if study.best_trial else "N/A"
                console.print(f"  Trial {trial.number}/{self.config.n_trials}: melhor={best_value}")

        # Executar otimização
        try:
            self.study.optimize(
                self._objective,
                n_trials=self.config.n_trials,
                n_jobs=self.config.n_jobs,
                callbacks=[progress_callback] if show_progress else None,
                show_progress_bar=False,  # Usamos nosso próprio progresso
            )
        except KeyboardInterrupt:
            console.print("\n⚠️  Otimização interrompida pelo usuário")

        optimization_time = time.time() - start_time

        # Coletar resultados
        all_trials = []
        for trial in self.study.trials:
            trial_data = {
                "number": trial.number,
                "value": trial.value,
                "params": trial.params,
                "state": trial.state.name,
                "datetime_start": trial.datetime_start,
                "datetime_complete": trial.datetime_complete,
                "duration": trial.duration.total_seconds() if trial.duration else None,
            }
            all_trials.append(trial_data)

        result = OptimizationResult(
            best_params=self.study.best_params if self.study.best_trial else {},
            best_value=self.study.best_value if self.study.best_trial else float("inf"),
            n_trials=len(self.study.trials),
            study_name=self.study.study_name,
            optimization_time=optimization_time,
            all_trials=all_trials,
        )

        if show_progress:
            console.print(f"\n✅ {result.summary()}")

        return result

    def get_optimization_history(self) -> list[dict]:
        """Retorna histórico de otimização."""
        if not self.study:
            return []

        history = []
        for trial in self.study.trials:
            if trial.value is not None:
                history.append({"trial": trial.number, "value": trial.value, "params": trial.params.copy()})

        return history


def optimize_algorithm(
    algorithm_name: str,
    sequences: list[str],
    alphabet: str,
    n_trials: int = 100,
    timeout_per_trial: int = 300,
    show_progress: bool = True,
) -> OptimizationResult:
    """
    Função conveniente para otimizar um algoritmo.

    Args:
        algorithm_name: Nome do algoritmo a otimizar
        sequences: Sequências do dataset
        alphabet: Alfabeto usado
        n_trials: Número de trials
        timeout_per_trial: Timeout por trial em segundos
        show_progress: Mostrar progresso

    Returns:
        Resultado da otimização
    """

    config = OptimizationConfig(
        algorithm_name=algorithm_name,
        dataset_sequences=sequences,
        alphabet=alphabet,
        n_trials=n_trials,
        timeout_per_trial=timeout_per_trial,
    )

    optimizer = OptunaOptimizer(config)
    return optimizer.optimize(show_progress=show_progress)


def create_optimization_study(
    algorithm_name: str, study_name: str | None = None, direction: str = "minimize"
) -> optuna.Study:
    """
    Cria um estudo de otimização Optuna.

    Args:
        algorithm_name: Nome do algoritmo
        study_name: Nome do estudo (opcional)
        direction: Direção da otimização

    Returns:
        Estudo Optuna criado
    """

    name = study_name or f"{algorithm_name}_study"

    return optuna.create_study(direction=direction, sampler=TPESampler(seed=42), pruner=MedianPruner(), study_name=name)
