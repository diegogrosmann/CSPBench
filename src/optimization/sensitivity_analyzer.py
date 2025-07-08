"""
Módulo de análise de sensibilidade para algoritmos CSP usando SALib.

Este módulo fornece funcionalidades para analisar a sensibilidade dos algoritmos
CSP aos seus hiperparâmetros, utilizando diferentes métodos de análise como
Sobol, Morris e FAST.

Classes:
    SensitivityConfig: Configuração para análise de sensibilidade
    SensitivityResult: Resultado da análise de sensibilidade
    SensitivityAnalyzer: Classe principal para análise

Funções:
    analyze_algorithm_sensitivity: Função principal para análise
    create_parameter_space: Cria espaço de parâmetros para análise
"""

import logging
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from SALib.analyze import fast, morris, sobol
from SALib.sample import fast_sampler
from SALib.sample import morris as morris_sample
from SALib.sample import saltelli
from tqdm import tqdm

from algorithms.base import global_registry

logger = logging.getLogger(__name__)


@dataclass
class SensitivityConfig:
    """Configuração para análise de sensibilidade."""

    n_samples: int = 1000
    method: str = "sobol"  # "sobol", "morris", "fast"
    timeout_per_sample: Optional[float] = 60.0
    show_progress: bool = True
    seed: Optional[int] = None

    # Parâmetros específicos para cada método
    sobol_calc_second_order: bool = True
    morris_num_levels: int = 4
    morris_grid_jump: int = 2
    fast_m: int = 4


@dataclass
class SensitivityResult:
    """Resultado da análise de sensibilidade."""

    parameter_names: List[str]
    method: str
    n_samples: int
    analysis_time: float

    # Índices de sensibilidade (específicos para cada método)
    first_order: Optional[Dict[str, float]] = None
    second_order: Optional[Dict[str, float]] = None
    total_order: Optional[Dict[str, float]] = None

    # Para Morris
    mu: Optional[Dict[str, float]] = None
    mu_star: Optional[Dict[str, float]] = None
    sigma: Optional[Dict[str, float]] = None

    # Para FAST
    main_effect: Optional[Dict[str, float]] = None

    # Dados brutos
    problem: Optional[Dict[str, Any]] = None
    samples: Optional[np.ndarray] = None
    outputs: Optional[np.ndarray] = None


class SensitivityAnalyzer:
    """Analisador de sensibilidade para algoritmos CSP."""

    def __init__(self, config: SensitivityConfig):
        self.config = config
        self.algorithm_class: Optional[type] = None
        self.sequences: List[str] = []
        self.alphabet: str = ""

    def _get_parameter_space(self, algorithm_name: str) -> Dict[str, Any]:
        """Obtém espaço de parâmetros para o algoritmo."""
        # Mapear algoritmos para seus espaços de parâmetros
        parameter_spaces = {
            "BLF-GA": {
                "pop_size": [30, 200],
                "max_gens": [50, 500],
                "cross_prob": [0.6, 0.95],
                "mut_prob": [0.01, 0.3],
                "elite_rate": [0.01, 0.15],
                "tournament_k": [2, 8],
                "immigrant_freq": [5, 20],
                "immigrant_ratio": [0.1, 0.4],
                "diversity_threshold": [0.2, 0.8],
                "no_improve_patience": [0.1, 0.5],
                "restart_patience": [10, 50],
                "restart_ratio": [0.2, 0.6],
                "mutation_multi_n": [1, 5],
                "refine_iter_limit": [50, 200],
                "niching_radius": [2, 8],
                "mutation_adapt_factor": [1.5, 3.0],
                "mutation_adapt_duration": [3, 10],
                "disable_elitism_gens": [3, 10],
            },
            "CSC": {
                "max_iter": [100, 1000],
                "patience": [10, 100],
                "min_improvement": [1e-6, 1e-3],
                "restart_patience": [20, 100],
                "max_restarts": [1, 10],
            },
            "H3-CSP": {
                "beam_width": [5, 50],
                "max_iterations": [50, 500],
                "diversity_factor": [0.1, 0.9],
                "local_search_iters": [10, 100],
                "restart_threshold": [10, 100],
                "max_restarts": [1, 10],
            },
            "DP-CSP": {
                "max_depth": [5, 20],
                "pruning_threshold": [0.1, 0.9],
                "memory_limit": [100, 1000],
            },
        }

        return parameter_spaces.get(algorithm_name, {})

    def _create_problem(self, algorithm_name: str) -> Dict[str, Any]:
        """Cria problema SALib para análise."""
        param_space = self._get_parameter_space(algorithm_name)

        if not param_space:
            raise ValueError(f"Espaço de parâmetros não definido para {algorithm_name}")

        problem = {
            "num_vars": len(param_space),
            "names": list(param_space.keys()),
            "bounds": list(param_space.values()),
        }

        return problem

    def _generate_samples(self, problem: Dict[str, Any]) -> np.ndarray:
        """Gera amostras baseado no método escolhido."""
        if self.config.method == "sobol":
            return saltelli.sample(
                problem,
                self.config.n_samples,
                calc_second_order=self.config.sobol_calc_second_order,
                seed=self.config.seed,
            )
        elif self.config.method == "morris":
            return morris_sample.sample(
                problem,
                self.config.n_samples,
                num_levels=self.config.morris_num_levels,
                grid_jump=self.config.morris_grid_jump,
                seed=self.config.seed,
            )
        elif self.config.method == "fast":
            return fast_sampler.sample(
                problem,
                self.config.n_samples,
                M=self.config.fast_m,
                seed=self.config.seed,
            )
        else:
            raise ValueError(f"Método desconhecido: {self.config.method}")

    def _evaluate_samples(
        self, samples: np.ndarray, problem: Dict[str, Any]
    ) -> np.ndarray:
        """Avalia amostras executando o algoritmo."""
        outputs = []
        param_names = problem["names"]

        progress_bar = None
        if self.config.show_progress:
            progress_bar = tqdm(total=len(samples), desc="Avaliando amostras")

        for i, sample in enumerate(samples):
            try:
                # Converter amostra para dicionário de parâmetros
                params = {name: value for name, value in zip(param_names, sample)}

                # Criar instância do algoritmo
                algorithm = self.algorithm_class(self.sequences, self.alphabet)

                # Aplicar parâmetros
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
                center, distance = algorithm.run()
                execution_time = time.time() - start_time

                # Verificar timeout
                if (
                    self.config.timeout_per_sample
                    and execution_time > self.config.timeout_per_sample
                ):
                    logger.warning(
                        f"Amostra {i} excedeu timeout: {execution_time:.2f}s"
                    )
                    distance = float("inf")

                outputs.append(distance)

            except Exception as e:
                logger.error(f"Erro na amostra {i}: {e}")
                outputs.append(float("inf"))

            if progress_bar:
                progress_bar.update(1)

        if progress_bar:
            progress_bar.close()

        return np.array(outputs)

    def _analyze_results(
        self, problem: Dict[str, Any], outputs: np.ndarray
    ) -> Dict[str, Any]:
        """Analisa resultados usando método escolhido."""
        if self.config.method == "sobol":
            return sobol.analyze(
                problem,
                outputs,
                calc_second_order=self.config.sobol_calc_second_order,
                seed=self.config.seed,
            )
        elif self.config.method == "morris":
            return morris.analyze(problem, self.samples, outputs, seed=self.config.seed)
        elif self.config.method == "fast":
            return fast.analyze(
                problem, outputs, M=self.config.fast_m, seed=self.config.seed
            )
        else:
            raise ValueError(f"Método desconhecido: {self.config.method}")

    def analyze(
        self, algorithm_name: str, sequences: List[str], alphabet: str
    ) -> SensitivityResult:
        """Executa análise de sensibilidade."""
        logger.info(f"Iniciando análise de sensibilidade para {algorithm_name}")

        # Verificar se algoritmo existe
        if algorithm_name not in global_registry:
            raise ValueError(f"Algoritmo não encontrado: {algorithm_name}")

        self.algorithm_class = global_registry[algorithm_name]
        self.sequences = sequences
        self.alphabet = alphabet

        # Criar problema
        problem = self._create_problem(algorithm_name)
        logger.info(f"Analisando {problem['num_vars']} parâmetros: {problem['names']}")

        # Gerar amostras
        logger.info(f"Gerando amostras usando método {self.config.method}")
        samples = self._generate_samples(problem)
        self.samples = samples
        logger.info(f"Geradas {len(samples)} amostras")

        # Avaliar amostras
        start_time = time.time()
        outputs = self._evaluate_samples(samples, problem)

        # Filtrar valores inválidos
        valid_indices = np.isfinite(outputs)
        if not np.all(valid_indices):
            logger.warning(f"Removendo {np.sum(~valid_indices)} amostras inválidas")
            samples = samples[valid_indices]
            outputs = outputs[valid_indices]

        # Executar análise
        logger.info("Executando análise de sensibilidade")
        analysis_results = self._analyze_results(problem, outputs)

        analysis_time = time.time() - start_time

        # Criar resultado
        result = SensitivityResult(
            parameter_names=problem["names"],
            method=self.config.method,
            n_samples=len(outputs),
            analysis_time=analysis_time,
            problem=problem,
            samples=samples,
            outputs=outputs,
        )

        # Extrair índices específicos para cada método
        if self.config.method == "sobol":
            result.first_order = {
                name: value
                for name, value in zip(problem["names"], analysis_results["S1"])
            }
            result.total_order = {
                name: value
                for name, value in zip(problem["names"], analysis_results["ST"])
            }
            if self.config.sobol_calc_second_order:
                # S2 é uma matriz, precisamos processar os índices de segunda ordem
                s2_matrix = analysis_results["S2"]
                result.second_order = {}
                for i in range(len(problem["names"])):
                    for j in range(i + 1, len(problem["names"])):
                        param_pair = f"{problem['names'][i]}-{problem['names'][j]}"
                        result.second_order[param_pair] = s2_matrix[i, j]

        elif self.config.method == "morris":
            result.mu = {
                name: value
                for name, value in zip(problem["names"], analysis_results["mu"])
            }
            result.mu_star = {
                name: value
                for name, value in zip(problem["names"], analysis_results["mu_star"])
            }
            result.sigma = {
                name: value
                for name, value in zip(problem["names"], analysis_results["sigma"])
            }

        elif self.config.method == "fast":
            result.first_order = {
                name: value
                for name, value in zip(problem["names"], analysis_results["S1"])
            }
            result.total_order = {
                name: value
                for name, value in zip(problem["names"], analysis_results["ST"])
            }

        logger.info(f"Análise concluída em {analysis_time:.2f}s")

        return result


def create_parameter_space(algorithm_name: str) -> Dict[str, List[float]]:
    """
    Cria espaço de parâmetros para análise de sensibilidade.

    Args:
        algorithm_name: Nome do algoritmo

    Returns:
        Dict com nomes dos parâmetros e seus intervalos
    """
    analyzer = SensitivityAnalyzer(SensitivityConfig())
    return analyzer._get_parameter_space(algorithm_name)


def analyze_algorithm_sensitivity(
    algorithm_name: str,
    sequences: List[str],
    alphabet: str,
    n_samples: int = 1000,
    method: str = "sobol",
    timeout_per_sample: Optional[float] = 60.0,
    show_progress: bool = True,
    seed: Optional[int] = None,
    **kwargs,
) -> SensitivityResult:
    """
    Analisa sensibilidade de um algoritmo CSP aos seus hiperparâmetros.

    Args:
        algorithm_name: Nome do algoritmo a ser analisado
        sequences: Lista de sequências do dataset
        alphabet: Alfabeto usado nas sequências
        n_samples: Número de amostras para análise
        method: Método de análise ("sobol", "morris", "fast")
        timeout_per_sample: Timeout por amostra (segundos)
        show_progress: Mostrar barra de progresso
        seed: Seed para reprodutibilidade
        **kwargs: Parâmetros adicionais específicos do método

    Returns:
        SensitivityResult: Resultado da análise
    """

    # Criar configuração
    config = SensitivityConfig(
        n_samples=n_samples,
        method=method,
        timeout_per_sample=timeout_per_sample,
        show_progress=show_progress,
        seed=seed,
    )

    # Aplicar parâmetros específicos do método
    if method == "sobol":
        config.sobol_calc_second_order = kwargs.get("calc_second_order", True)
    elif method == "morris":
        config.morris_num_levels = kwargs.get("num_levels", 4)
        config.morris_grid_jump = kwargs.get("grid_jump", 2)
    elif method == "fast":
        config.fast_m = kwargs.get("M", 4)

    # Criar analisador
    analyzer = SensitivityAnalyzer(config)

    # Executar análise
    result = analyzer.analyze(algorithm_name, sequences, alphabet)

    return result
