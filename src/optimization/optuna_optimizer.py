"""
Módulo de otimização de hiperparâmetros usando Optuna.

Este módulo fornece funcionalidades para otimizar hiperparâmetros de algoritmos
CSP usando o framework Optuna, incluindo diferentes estratégias de amostragem,
poda de trials e salvamento de resultados.

Classes:
    OptimizationConfig: Configuração para otimização
    OptimizationResult: Resultado da otimização
    OptunaOptimizer: Classe principal para otimização

Funções:
    optimize_algorithm: Função principal para otimizar algoritmos
    create_optimization_study: Cria um estudo Optuna
"""

import json
import logging
import os
import time
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict, List, Optional, Union

import optuna
from optuna.pruners import MedianPruner, SuccessiveHalvingPruner
from optuna.samplers import CmaEsSampler, RandomSampler, TPESampler
from tqdm import tqdm

from algorithms.base import global_registry

logger = logging.getLogger(__name__)


@dataclass
class OptimizationConfig:
    """Configuração para otimização de hiperparâmetros."""

    n_trials: int = 100
    timeout: Optional[float] = None
    timeout_per_trial: Optional[float] = 60.0
    direction: str = "minimize"  # "minimize" ou "maximize"
    sampler: str = "TPE"  # "TPE", "CmaEs", "Random"
    pruner: str = "Median"  # "Median", "SuccessiveHalving", None
    study_name: Optional[str] = None
    storage: Optional[str] = None
    load_if_exists: bool = True
    show_progress: bool = True
    n_startup_trials: int = 10
    n_ei_candidates: int = 24
    seed: Optional[int] = None


@dataclass
class OptimizationResult:
    """Resultado da otimização de hiperparâmetros."""

    best_params: Dict[str, Any]
    best_value: float
    n_trials: int
    study_name: str
    optimization_time: float
    all_trials: List[Dict[str, Any]]
    study: optuna.Study

    def __post_init__(self):
        """Processa trials após inicialização."""
        if not self.all_trials:
            self.all_trials = []
            for trial in self.study.trials:
                trial_data = {
                    "number": trial.number,
                    "value": trial.value,
                    "params": trial.params,
                    "state": trial.state.name,
                    "datetime_start": trial.datetime_start,
                    "datetime_complete": trial.datetime_complete,
                }
                if trial.datetime_start and trial.datetime_complete:
                    trial_data["duration"] = (
                        trial.datetime_complete - trial.datetime_start
                    ).total_seconds()
                self.all_trials.append(trial_data)


class OptunaOptimizer:
    """Otimizador de hiperparâmetros usando Optuna."""

    def __init__(self, config: OptimizationConfig):
        self.config = config
        self.study: Optional[optuna.Study] = None
        self.algorithm_class: Optional[type] = None
        self.sequences: List[str] = []
        self.alphabet: str = ""

    def _create_sampler(self) -> optuna.samplers.BaseSampler:
        """Cria sampler baseado na configuração."""
        if self.config.sampler == "TPE":
            return TPESampler(
                n_startup_trials=self.config.n_startup_trials,
                n_ei_candidates=self.config.n_ei_candidates,
                seed=self.config.seed,
            )
        elif self.config.sampler == "CmaEs":
            return CmaEsSampler(seed=self.config.seed)
        elif self.config.sampler == "Random":
            return RandomSampler(seed=self.config.seed)
        else:
            logger.warning(f"Sampler desconhecido: {self.config.sampler}. Usando TPE.")
            return TPESampler(seed=self.config.seed)

    def _create_pruner(self) -> Optional[optuna.pruners.BasePruner]:
        """Cria pruner baseado na configuração."""
        if self.config.pruner == "Median":
            return MedianPruner(n_startup_trials=self.config.n_startup_trials)
        elif self.config.pruner == "SuccessiveHalving":
            return SuccessiveHalvingPruner()
        else:
            return None

    def _get_parameter_space(self, algorithm_name: str) -> Dict[str, Any]:
        """Obtém espaço de parâmetros para o algoritmo."""
        # Mapear algoritmos para seus espaços de parâmetros
        parameter_spaces = {
            "BLF-GA": {
                "pop_size": ("int", 30, 200),
                "max_gens": ("int", 50, 500),
                "cross_prob": ("float", 0.6, 0.95),
                "mut_prob": ("float", 0.01, 0.3),
                "elite_rate": ("float", 0.01, 0.15),
                "tournament_k": ("int", 2, 8),
                "immigrant_freq": ("int", 5, 20),
                "immigrant_ratio": ("float", 0.1, 0.4),
                "diversity_threshold": ("float", 0.2, 0.8),
                "no_improve_patience": ("float", 0.1, 0.5),
                "restart_patience": ("int", 10, 50),
                "restart_ratio": ("float", 0.2, 0.6),
                "crossover_type": (
                    "categorical",
                    ["one_point", "uniform", "blend_blocks"],
                ),
                "mutation_type": (
                    "categorical",
                    ["multi", "inversion", "transposition"],
                ),
                "mutation_multi_n": ("int", 1, 5),
                "refinement_type": (
                    "categorical",
                    ["greedy", "swap", "insertion", "2opt"],
                ),
                "refine_elites": ("categorical", ["best", "all"]),
                "refine_iter_limit": ("int", 50, 200),
                "niching": ("categorical", [True, False]),
                "niching_radius": ("int", 2, 8),
                "mutation_adapt_factor": ("float", 1.5, 3.0),
                "mutation_adapt_duration": ("int", 3, 10),
                "disable_elitism_gens": ("int", 3, 10),
            },
            "CSC": {
                "max_iter": ("int", 100, 1000),
                "patience": ("int", 10, 100),
                "min_improvement": ("float", 1e-6, 1e-3),
                "random_restart": ("categorical", [True, False]),
                "restart_patience": ("int", 20, 100),
                "max_restarts": ("int", 1, 10),
            },
            "H3-CSP": {
                "beam_width": ("int", 5, 50),
                "max_iterations": ("int", 50, 500),
                "diversity_factor": ("float", 0.1, 0.9),
                "local_search": ("categorical", [True, False]),
                "local_search_iters": ("int", 10, 100),
                "restart_threshold": ("int", 10, 100),
                "max_restarts": ("int", 1, 10),
            },
            "DP-CSP": {
                "max_depth": ("int", 5, 20),
                "pruning_threshold": ("float", 0.1, 0.9),
                "use_heuristic": ("categorical", [True, False]),
                "memory_limit": ("int", 100, 1000),
            },
        }

        return parameter_spaces.get(algorithm_name, {})

    def _suggest_parameters(
        self, trial: optuna.Trial, algorithm_name: str
    ) -> Dict[str, Any]:
        """Sugere parâmetros para o trial."""
        param_space = self._get_parameter_space(algorithm_name)
        params = {}

        for param_name, param_config in param_space.items():
            param_type = param_config[0]

            if param_type == "int":
                low, high = param_config[1], param_config[2]
                params[param_name] = trial.suggest_int(param_name, low, high)
            elif param_type == "float":
                low, high = param_config[1], param_config[2]
                params[param_name] = trial.suggest_float(param_name, low, high)
            elif param_type == "categorical":
                choices = param_config[1]
                params[param_name] = trial.suggest_categorical(param_name, choices)
            elif param_type == "loguniform":
                low, high = param_config[1], param_config[2]
                params[param_name] = trial.suggest_float(
                    param_name, low, high, log=True
                )

        return params

    def _objective(self, trial: optuna.Trial, algorithm_name: str) -> float:
        """Função objetivo para otimização."""
        try:
            # Sugerir parâmetros
            params = self._suggest_parameters(trial, algorithm_name)

            # Criar instância do algoritmo
            if self.algorithm_class is None:
                raise ValueError(
                    f"Classe do algoritmo não encontrada: {algorithm_name}"
                )

            # Criar algoritmo com parâmetros sugeridos
            algorithm = self.algorithm_class(self.sequences, self.alphabet)

            # Aplicar parâmetros sugeridos
            if hasattr(algorithm, "set_params"):
                algorithm.set_params(**params)
            elif hasattr(algorithm, "config"):
                # Atualizar configuração do algoritmo
                if hasattr(algorithm.config, "update"):
                    algorithm.config.update(params)
                else:
                    for key, value in params.items():
                        if hasattr(algorithm.config, key):
                            setattr(algorithm.config, key, value)
            else:
                # Tentar definir parâmetros diretamente
                for key, value in params.items():
                    if hasattr(algorithm, key):
                        setattr(algorithm, key, value)

            # Executar algoritmo
            start_time = time.time()
            result = algorithm.run()
            execution_time = time.time() - start_time

            # Extrair center e distance do resultado
            if isinstance(result, tuple):
                if len(result) == 2:
                    center, distance = result
                elif len(result) == 3:
                    center, distance, metadata = result
                else:
                    raise ValueError(
                        f"Resultado inesperado do algoritmo: {len(result)} valores"
                    )
            else:
                raise ValueError(f"Tipo de resultado inesperado: {type(result)}")

            # Verificar timeout por trial
            if (
                self.config.timeout_per_trial
                and execution_time > self.config.timeout_per_trial
            ):
                raise optuna.TrialPruned(
                    f"Trial excedeu timeout: {execution_time:.2f}s"
                )

            # Armazenar informações adicionais
            trial.set_user_attr("execution_time", execution_time)
            trial.set_user_attr("center", center)
            trial.set_user_attr("distance", distance)

            # Retornar valor baseado na direção
            if self.config.direction == "minimize":
                return distance
            else:
                return -distance

        except Exception as e:
            logger.error(f"Erro no trial {trial.number}: {e}")
            # Retornar valor que indica falha
            if self.config.direction == "minimize":
                return float("inf")
            else:
                return float("-inf")

    def optimize(
        self, algorithm_name: str, sequences: List[str], alphabet: str
    ) -> OptimizationResult:
        """Executa otimização de hiperparâmetros."""
        logger.info(f"Iniciando otimização para {algorithm_name}")

        # Verificar se algoritmo existe
        if algorithm_name not in global_registry:
            raise ValueError(f"Algoritmo não encontrado: {algorithm_name}")

        self.algorithm_class = global_registry[algorithm_name]
        self.sequences = sequences
        self.alphabet = alphabet

        # Criar sampler e pruner
        sampler = self._create_sampler()
        pruner = self._create_pruner()

        # Criar estudo
        study_name = (
            self.config.study_name or f"optimize_{algorithm_name}_{int(time.time())}"
        )

        self.study = optuna.create_study(
            direction=self.config.direction,
            sampler=sampler,
            pruner=pruner,
            study_name=study_name,
            storage=self.config.storage,
            load_if_exists=self.config.load_if_exists,
        )

        # Executar otimização
        start_time = time.time()

        if self.config.show_progress:
            # Usar tqdm para mostrar progresso
            with tqdm(total=self.config.n_trials, desc="Otimizando") as pbar:

                def callback(study, trial):
                    pbar.set_description(f"Trial {trial.number}: {trial.value:.6f}")
                    pbar.update(1)

                self.study.optimize(
                    lambda trial: self._objective(trial, algorithm_name),
                    n_trials=self.config.n_trials,
                    timeout=self.config.timeout,
                    callbacks=[callback],
                )
        else:
            self.study.optimize(
                lambda trial: self._objective(trial, algorithm_name),
                n_trials=self.config.n_trials,
                timeout=self.config.timeout,
            )

        optimization_time = time.time() - start_time

        logger.info(f"Otimização concluída em {optimization_time:.2f}s")
        logger.info(f"Melhor valor: {self.study.best_value}")
        logger.info(f"Melhores parâmetros: {self.study.best_params}")

        # Criar resultado
        result = OptimizationResult(
            best_params=self.study.best_params,
            best_value=self.study.best_value,
            n_trials=len(self.study.trials),
            study_name=study_name,
            optimization_time=optimization_time,
            all_trials=[],
            study=self.study,
        )

        return result


def create_optimization_study(
    algorithm_name: str,
    direction: str = "minimize",
    sampler: str = "TPE",
    pruner: str = "Median",
    study_name: Optional[str] = None,
    storage: Optional[str] = None,
    load_if_exists: bool = True,
    seed: Optional[int] = None,
) -> optuna.Study:
    """Cria um estudo Optuna para otimização."""

    # Criar sampler
    if sampler == "TPE":
        sampler_obj = TPESampler(seed=seed)
    elif sampler == "CmaEs":
        sampler_obj = CmaEsSampler(seed=seed)
    elif sampler == "Random":
        sampler_obj = RandomSampler(seed=seed)
    else:
        logger.warning(f"Sampler desconhecido: {sampler}. Usando TPE.")
        sampler_obj = TPESampler(seed=seed)

    # Criar pruner
    if pruner == "Median":
        pruner_obj = MedianPruner()
    elif pruner == "SuccessiveHalving":
        pruner_obj = SuccessiveHalvingPruner()
    else:
        pruner_obj = None

    # Nome do estudo
    if not study_name:
        study_name = f"optimize_{algorithm_name}_{int(time.time())}"

    # Criar estudo
    study = optuna.create_study(
        direction=direction,
        sampler=sampler_obj,
        pruner=pruner_obj,
        study_name=study_name,
        storage=storage,
        load_if_exists=load_if_exists,
    )

    return study


def optimize_algorithm(
    algorithm_name: str,
    sequences: List[str],
    alphabet: str,
    n_trials: int = 100,
    timeout: Optional[float] = None,
    timeout_per_trial: Optional[float] = 60.0,
    direction: str = "minimize",
    sampler: str = "TPE",
    pruner: str = "Median",
    study_name: Optional[str] = None,
    storage: Optional[str] = None,
    load_if_exists: bool = True,
    show_progress: bool = True,
    seed: Optional[int] = None,
) -> OptimizationResult:
    """
    Otimiza hiperparâmetros de um algoritmo CSP.

    Args:
        algorithm_name: Nome do algoritmo a ser otimizado
        sequences: Lista de sequências do dataset
        alphabet: Alfabeto usado nas sequências
        n_trials: Número de trials para otimização
        timeout: Timeout total para otimização (segundos)
        timeout_per_trial: Timeout por trial (segundos)
        direction: Direção da otimização ("minimize" ou "maximize")
        sampler: Tipo de sampler ("TPE", "CmaEs", "Random")
        pruner: Tipo de pruner ("Median", "SuccessiveHalving", None)
        study_name: Nome do estudo
        storage: URL de armazenamento (opcional)
        load_if_exists: Carregar estudo existente se existir
        show_progress: Mostrar barra de progresso
        seed: Seed para reprodutibilidade

    Returns:
        OptimizationResult: Resultado da otimização
    """

    # Criar configuração
    config = OptimizationConfig(
        n_trials=n_trials,
        timeout=timeout,
        timeout_per_trial=timeout_per_trial,
        direction=direction,
        sampler=sampler,
        pruner=pruner,
        study_name=study_name,
        storage=storage,
        load_if_exists=load_if_exists,
        show_progress=show_progress,
        seed=seed,
    )

    # Criar otimizador
    optimizer = OptunaOptimizer(config)

    # Executar otimização
    result = optimizer.optimize(algorithm_name, sequences, alphabet)

    return result


def run_optimization_with_dataset_selection():
    """
    Executa otimização com seleção interativa de dataset.
    """
    try:
        from src.ui.cli.menu import (
            configure_optimization_params,
            select_dataset_for_optimization,
            select_optimization_algorithm,
        )

        # Selecionar dataset
        sequences, alphabet, dataset_info = select_dataset_for_optimization()

        # Selecionar algoritmo
        algorithm_name = select_optimization_algorithm()

        # Configurar parâmetros
        config_dict = configure_optimization_params()

        # Executar otimização
        print(f"\n🚀 Iniciando otimização do {algorithm_name}...")
        print(
            f"📊 Dataset: {dataset_info.get('type', 'N/A')} - {len(sequences)} sequências"
        )
        print(f"🔬 Trials: {config_dict.get('n_trials', 100)}")

        result = optimize_algorithm(
            algorithm_name=algorithm_name,
            sequences=sequences,
            alphabet=alphabet,
            n_trials=config_dict.get("n_trials", 100),
            timeout_per_trial=config_dict.get("timeout_per_trial", 60),
            direction=config_dict.get("direction", "minimize"),
            sampler=config_dict.get("sampler", "TPE"),
            pruner=config_dict.get("pruner", "Median"),
            show_progress=True,
        )

        # Exibir resultados
        print(f"\n✅ Otimização concluída!")
        print(f"🎯 Melhor valor: {result.best_value:.6f}")
        print(f"⏱️ Tempo total: {result.optimization_time:.2f} segundos")
        print(f"📈 Trials realizados: {result.n_trials}")

        # Salvar resultados
        os.makedirs("outputs/reports", exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"outputs/reports/optimization_{algorithm_name}_{timestamp}.json"

        report = {
            "algorithm": algorithm_name,
            "dataset_info": dataset_info,
            "best_params": result.best_params,
            "best_value": result.best_value,
            "n_trials": result.n_trials,
            "optimization_time": result.optimization_time,
            "study_name": result.study_name,
            "timestamp": timestamp,
        }

        with open(filename, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, default=str, ensure_ascii=False)

        print(f"💾 Relatório salvo em: {filename}")

        # Opção de gerar gráficos
        plots_input = input("Gerar gráficos de visualização? (s/N): ").strip()
        if plots_input.lower() in ["s", "sim", "y", "yes"]:
            try:
                from src.optimization.visualization import OptimizationVisualizer

                visualizer = OptimizationVisualizer(result)

                # Gráfico de histórico
                history_path = f"outputs/reports/optimization_{algorithm_name}_{timestamp}_history.png"
                visualizer.plot_optimization_history(
                    save_path=history_path, interactive=False
                )

                # Gráfico de importância
                importance_path = f"outputs/reports/optimization_{algorithm_name}_{timestamp}_importance.png"
                visualizer.plot_parameter_importance(
                    save_path=importance_path, interactive=False
                )

                print(f"📊 Gráficos salvos em outputs/reports/")

            except Exception as e:
                print(f"⚠️ Erro ao gerar gráficos: {e}")

    except Exception as e:
        print(f"❌ Erro na otimização: {e}")
        import traceback

        traceback.print_exc()
