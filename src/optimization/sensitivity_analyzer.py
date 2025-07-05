"""
Módulo de análise de sensibilidade usando SALib.

Classes:
    SensitivityAnalyzer: Analisador principal usando SALib
    SensitivityConfig: Configuração da análise
    SensitivityResult: Resultado da análise

Funções:
    analyze_algorithm_sensitivity: Analisa sensibilidade de um algoritmo
    create_parameter_space: Cria espaço de parâmetros para análise
"""

import logging
import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from SALib.analyze import sobol
from SALib.sample import sobol as sobol_sample

from algorithms.base import global_registry
from src.core.exec.algorithm_executor import AlgorithmExecutor
from src.ui.cli.console_manager import console

logger = logging.getLogger(__name__)


@dataclass
class SensitivityConfig:
    """Configuração para análise de sensibilidade."""

    algorithm_name: str
    dataset_sequences: list[str]
    alphabet: str
    n_samples: int = 1000
    timeout_per_sample: int = 60
    method: str = "sobol"  # sobol, morris, fast
    calc_second_order: bool = True

    def __post_init__(self):
        """Validação da configuração."""
        if self.algorithm_name not in global_registry:
            raise ValueError(f"Algoritmo não encontrado: {self.algorithm_name}")

        if self.method not in ["sobol", "morris", "fast"]:
            raise ValueError("method deve ser 'sobol', 'morris' ou 'fast'")


@dataclass
class SensitivityResult:
    """Resultado da análise de sensibilidade."""

    method: str
    parameter_names: list[str]
    first_order: dict[str, float]
    total_order: dict[str, float]
    second_order: dict[str, dict[str, float]] | None
    confidence_intervals: dict[str, dict[str, float]]
    n_samples: int
    analysis_time: float

    def summary(self) -> str:
        """Retorna resumo da análise."""
        lines = [
            f"Análise de Sensibilidade ({self.method.upper()}):",
            f"  Samples: {self.n_samples}",
            f"  Tempo: {self.analysis_time:.1f}s",
            "",
            "Índices de Primeira Ordem (S1):",
        ]

        # Ordenar por importância
        sorted_first = sorted(self.first_order.items(), key=lambda x: x[1], reverse=True)
        for param, value in sorted_first:
            lines.append(f"  {param:<20}: {value:.4f}")

        lines.append("\nÍndices Totais (ST):")
        sorted_total = sorted(self.total_order.items(), key=lambda x: x[1], reverse=True)
        for param, value in sorted_total:
            lines.append(f"  {param:<20}: {value:.4f}")

        return "\n".join(lines)

    def get_most_important_parameters(self, n: int = 5) -> list[tuple[str, float]]:
        """Retorna os n parâmetros mais importantes baseado nos índices totais."""
        sorted_params = sorted(self.total_order.items(), key=lambda x: x[1], reverse=True)
        return sorted_params[:n]


class SensitivityAnalyzer:
    """Analisador de sensibilidade usando SALib."""

    def __init__(self, config: SensitivityConfig):
        self.config = config
        self.algorithm_class = global_registry[config.algorithm_name]
        self.parameter_space = self._create_parameter_space()

        logger.info(f"Criado analisador de sensibilidade para {config.algorithm_name}")

    def _create_parameter_space(self) -> dict[str, Any]:
        """Cria espaço de parâmetros baseado no algoritmo."""

        if self.config.algorithm_name == "BLF-GA":
            return {
                "num_vars": 6,
                "names": [
                    "population_size",
                    "max_generations",
                    "crossover_rate",
                    "mutation_rate",
                    "elite_size",
                    "tournament_size",
                ],
                "bounds": [
                    [20, 200],  # population_size
                    [100, 1000],  # max_generations
                    [0.6, 0.95],  # crossover_rate
                    [0.01, 0.3],  # mutation_rate
                    [1, 10],  # elite_size
                    [2, 8],  # tournament_size
                ],
            }

        elif self.config.algorithm_name == "CSC":
            return {
                "num_vars": 3,
                "names": ["max_iterations", "improvement_threshold", "restart_threshold"],
                "bounds": [
                    [100, 2000],  # max_iterations
                    [1e-6, 1e-3],  # improvement_threshold (log scale)
                    [10, 100],  # restart_threshold
                ],
            }

        elif self.config.algorithm_name == "DP-CSP":
            return {
                "num_vars": 2,
                "names": ["beam_width", "pruning_factor"],
                "bounds": [[10, 100], [0.1, 0.9]],  # beam_width  # pruning_factor
            }

        elif self.config.algorithm_name == "H3-CSP":
            return {
                "num_vars": 3,
                "names": ["max_depth", "branching_factor", "heuristic_weight"],
                "bounds": [[5, 20], [2, 10], [0.1, 2.0]],  # max_depth  # branching_factor  # heuristic_weight
            }

        # Algoritmo genérico
        return {"num_vars": 1, "names": ["seed"], "bounds": [[1, 1000000]]}

    def _evaluate_model(self, parameter_sets: np.ndarray) -> np.ndarray:
        """Avalia o modelo para cada conjunto de parâmetros."""

        results = []
        total_sets = len(parameter_sets)

        console.print(f"🔬 Avaliando {total_sets} conjuntos de parâmetros...")

        for i, params in enumerate(parameter_sets):
            if i % 100 == 0:
                console.print(f"  Progresso: {i}/{total_sets} ({i/total_sets*100:.1f}%)")

            try:
                # Converter parâmetros para dicionário
                param_dict = {}
                for j, name in enumerate(self.parameter_space["names"]):
                    value = params[j]

                    # Converter para tipos apropriados
                    if name in [
                        "population_size",
                        "max_generations",
                        "elite_size",
                        "tournament_size",
                        "max_iterations",
                        "restart_threshold",
                        "beam_width",
                        "max_depth",
                        "branching_factor",
                        "seed",
                    ]:
                        value = int(value)
                    elif name == "improvement_threshold":
                        value = float(10**value)  # Log scale
                    else:
                        value = float(value)

                    param_dict[name] = value

                # Criar instância do algoritmo
                algorithm_instance = self.algorithm_class(
                    self.config.dataset_sequences, self.config.alphabet, **param_dict
                )

                # Executar algoritmo
                executor = AlgorithmExecutor(self.config.timeout_per_sample)
                result, _, _ = executor.execute_with_timeout(algorithm_instance)

                if result is None or "distancia" not in result:
                    distance = float("inf")
                else:
                    distance = result["distancia"]

                results.append(distance)

            except Exception as e:
                logger.warning(f"Erro na avaliação {i}: {e}")
                results.append(float("inf"))

        return np.array(results)

    def analyze(self, show_progress: bool = True) -> SensitivityResult:
        """Executa análise de sensibilidade."""

        if show_progress:
            console.print(f"\n🔬 Iniciando análise de sensibilidade de {self.config.algorithm_name}")
            console.print(
                f"📊 Dataset: n={len(self.config.dataset_sequences)}, L={len(self.config.dataset_sequences[0])}"
            )
            console.print(f"🎯 Samples: {self.config.n_samples}")
            console.print(f"📋 Parâmetros: {', '.join(self.parameter_space['names'])}")

        start_time = time.time()

        try:
            # Gerar samples usando Sobol
            if self.config.method == "sobol":
                param_values = sobol_sample.sample(
                    self.parameter_space, self.config.n_samples, calc_second_order=self.config.calc_second_order
                )
            else:
                # Para outros métodos, usar sobol como fallback
                param_values = sobol_sample.sample(self.parameter_space, self.config.n_samples)

            # Avaliar modelo
            model_outputs = self._evaluate_model(param_values)

            # Remover valores infinitos
            valid_mask = np.isfinite(model_outputs)
            if not np.any(valid_mask):
                raise ValueError("Todos os resultados são inválidos")

            valid_outputs = model_outputs[valid_mask]
            valid_params = param_values[valid_mask]

            if show_progress:
                console.print(f"📈 {len(valid_outputs)}/{len(model_outputs)} amostras válidas")

            # Análise de sensibilidade
            if self.config.method == "sobol":
                si = sobol.analyze(self.parameter_space, valid_outputs, calc_second_order=self.config.calc_second_order)
            else:
                # Fallback para Sobol
                si = sobol.analyze(self.parameter_space, valid_outputs)

            # Construir resultado
            first_order = {name: float(si["S1"][i]) for i, name in enumerate(self.parameter_space["names"])}

            total_order = {name: float(si["ST"][i]) for i, name in enumerate(self.parameter_space["names"])}

            # Intervalos de confiança
            confidence_intervals = {}
            if "S1_conf" in si:
                confidence_intervals["S1"] = {
                    name: float(si["S1_conf"][i]) for i, name in enumerate(self.parameter_space["names"])
                }
            if "ST_conf" in si:
                confidence_intervals["ST"] = {
                    name: float(si["ST_conf"][i]) for i, name in enumerate(self.parameter_space["names"])
                }

            # Segunda ordem se disponível
            second_order = None
            if self.config.calc_second_order and "S2" in si:
                second_order = {}
                names = self.parameter_space["names"]
                for i, name1 in enumerate(names):
                    second_order[name1] = {}
                    for j, name2 in enumerate(names):
                        if i != j:
                            second_order[name1][name2] = float(si["S2"][i, j])

        except KeyboardInterrupt:
            console.print("\n⚠️  Análise interrompida pelo usuário")
            raise
        except Exception as e:
            console.print(f"\n❌ Erro na análise: {e}")
            raise

        analysis_time = time.time() - start_time

        result = SensitivityResult(
            method=self.config.method,
            parameter_names=self.parameter_space["names"],
            first_order=first_order,
            total_order=total_order,
            second_order=second_order,
            confidence_intervals=confidence_intervals,
            n_samples=len(valid_outputs),
            analysis_time=analysis_time,
        )

        if show_progress:
            console.print(f"\n✅ {result.summary()}")

        return result


def analyze_algorithm_sensitivity(
    algorithm_name: str,
    sequences: list[str],
    alphabet: str,
    n_samples: int = 1000,
    timeout_per_sample: int = 60,
    show_progress: bool = True,
) -> SensitivityResult:
    """
    Função conveniente para análise de sensibilidade.

    Args:
        algorithm_name: Nome do algoritmo
        sequences: Sequências do dataset
        alphabet: Alfabeto usado
        n_samples: Número de amostras
        timeout_per_sample: Timeout por amostra
        show_progress: Mostrar progresso

    Returns:
        Resultado da análise
    """

    config = SensitivityConfig(
        algorithm_name=algorithm_name,
        dataset_sequences=sequences,
        alphabet=alphabet,
        n_samples=n_samples,
        timeout_per_sample=timeout_per_sample,
    )

    analyzer = SensitivityAnalyzer(config)
    return analyzer.analyze(show_progress=show_progress)


def create_parameter_space(algorithm_name: str) -> dict[str, Any]:
    """
    Cria espaço de parâmetros para um algoritmo.

    Args:
        algorithm_name: Nome do algoritmo

    Returns:
        Dicionário com espaço de parâmetros
    """

    config = SensitivityConfig(algorithm_name=algorithm_name, dataset_sequences=["ATCG"], alphabet="ATCG")  # Dummy

    analyzer = SensitivityAnalyzer(config)
    return analyzer.parameter_space
