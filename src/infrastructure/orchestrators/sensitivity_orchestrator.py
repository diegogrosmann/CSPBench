"""
Orchestrator para Análise de Sensibilidade

Implementa orquestração de análises de sensibilidade usando SALib.
"""

import logging
import time
from typing import Any, Dict, List, Optional

import numpy as np
from SALib.analyze import delta, fast, morris, sobol
from SALib.sample import saltelli
from SALib.sample.fast_sampler import sample as fast_sample
from SALib.sample.morris import sample as morris_sample

from src.domain import Dataset
from src.domain.errors import SensitivityExecutionError


class SensitivityOrchestrator:
    """Orchestrator para análises de sensibilidade com SALib."""

    def __init__(self, algorithm_registry, executor):
        """
        Inicializa orchestrator de sensibilidade.

        Args:
            algorithm_registry: Registry de algoritmos
            executor: Executor para execução de algoritmos
        """
        self._algorithm_registry = algorithm_registry
        self._executor = executor
        self._logger = logging.getLogger(__name__)

    def execute_sensitivity_analysis(
        self,
        algorithm_name: str,
        dataset: Dataset,
        sensitivity_config: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Executa análise de sensibilidade completa.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset para análise
            sensitivity_config: Configuração da análise

        Returns:
            Dict[str, Any]: Resultados da análise de sensibilidade

        Raises:
            SensitivityExecutionError: Se erro na execução
        """
        try:
            start_time = time.time()
            self._logger.info(
                "Iniciando análise de sensibilidade para %s", algorithm_name
            )

            # Extrair configurações
            analysis_method = sensitivity_config.get("analysis_method", "morris")
            parameters = sensitivity_config["parameters"]
            n_samples = sensitivity_config.get("n_samples", 1000)
            repetitions = sensitivity_config.get("repetitions_per_sample", 3)
            output_metrics = sensitivity_config.get(
                "output_metrics", ["distance", "execution_time"]
            )

            # Criar problema SALib
            problem_definition = self._create_salib_problem(parameters)

            # Gerar amostras
            samples = self._generate_samples(
                problem_definition, analysis_method, n_samples, sensitivity_config
            )

            self._logger.info("Geradas %d amostras para análise", len(samples))

            # Executar algoritmo para cada amostra
            results = self._execute_samples(
                algorithm_name, dataset, samples, problem_definition, repetitions
            )

            # Análise de sensibilidade
            sensitivity_results = self._analyze_sensitivity(
                problem_definition,
                analysis_method,
                results,
                output_metrics,
                sensitivity_config,
                samples,  # Passar amostras para análise delta
            )

            # Compilar resultados finais
            execution_time = time.time() - start_time
            final_results = {
                "algorithm": algorithm_name,
                "dataset_info": {
                    "size": dataset.size,
                    "length": len(dataset.sequences[0]) if dataset.sequences else 0,
                    "alphabet": dataset.alphabet,
                },
                "analysis_method": analysis_method,
                "n_samples": len(samples),
                "repetitions_per_sample": repetitions,
                "parameters_analyzed": list(parameters.keys()),
                "sensitivity_indices": sensitivity_results,
                "execution_time": execution_time,
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            }

            self._logger.info(
                "Análise de sensibilidade concluída em %.2fs", execution_time
            )
            return final_results

        except Exception as e:
            self._logger.error("Erro na análise de sensibilidade: %s", e)
            raise SensitivityExecutionError(
                f"Falha na análise de sensibilidade: {e}"
            ) from e

    def _create_salib_problem(self, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """
        Cria definição de problema para SALib.

        Args:
            parameters: Configuração dos parâmetros

        Returns:
            Dict[str, Any]: Definição do problema SALib
        """
        problem_def = {"num_vars": len(parameters), "names": [], "bounds": []}

        for param_name, param_config in parameters.items():
            problem_def["names"].append(param_name)

            param_type = param_config.get("type", "float")

            if param_type in ["integer", "float"]:
                bounds = param_config["bounds"]
                problem_def["bounds"].append(bounds)
            elif param_type == "categorical":
                # Para categóricos, mapear para índices
                values = param_config["values"]
                problem_def["bounds"].append([0, len(values) - 1])
            else:
                raise ValueError(f"Tipo de parâmetro não suportado: {param_type}")

        return problem_def

    def _generate_samples(
        self,
        problem_def: Dict[str, Any],
        analysis_method: str,
        n_samples: int,
        sensitivity_config: Dict[str, Any],
    ) -> np.ndarray:
        """
        Gera amostras para análise de sensibilidade.

        Args:
            problem_def: Definição do problema
            analysis_method: Método de análise
            n_samples: Número de amostras
            sensitivity_config: Configuração completa

        Returns:
            np.ndarray: Amostras geradas
        """
        if analysis_method == "morris":
            morris_config = sensitivity_config.get("method_config", {})
            optimal_trajectories = morris_config.get("optimal_trajectories")

            # Se optimal_trajectories for None ou inválido, não passar o parâmetro
            kwargs = {
                "N": n_samples,
                "num_levels": morris_config.get("num_levels", 4),
            }

            if optimal_trajectories is not None:
                # Validar range permitido pelo SALib
                if (
                    isinstance(optimal_trajectories, int)
                    and 2 <= optimal_trajectories <= 8
                ):
                    kwargs["optimal_trajectories"] = optimal_trajectories
                else:
                    self._logger.warning(
                        "optimal_trajectories deve ser um inteiro entre 2 e 8. Usando amostragem padrão."
                    )

            return morris_sample(problem_def, **kwargs)
        elif analysis_method == "sobol":
            # Sobol precisa de amostras Saltelli
            return saltelli.sample(problem_def, N=n_samples)
        elif analysis_method == "fast":
            return fast_sample(problem_def, N=n_samples)
        elif analysis_method == "delta":
            # Delta usa amostragem aleatória simples
            samples = np.random.random((n_samples, problem_def["num_vars"]))
            # Escalar para bounds
            for i, bounds in enumerate(problem_def["bounds"]):
                samples[:, i] = samples[:, i] * (bounds[1] - bounds[0]) + bounds[0]
            return samples
        else:
            raise ValueError(f"Método de análise não suportado: {analysis_method}")

    def _execute_samples(
        self,
        algorithm_name: str,
        dataset: Dataset,
        samples: np.ndarray,
        problem_def: Dict[str, Any],
        repetitions: int,
    ) -> Dict[str, List[float]]:
        """
        Executa algoritmo para todas as amostras.

        Args:
            algorithm_name: Nome do algoritmo
            dataset: Dataset
            samples: Amostras de parâmetros
            problem_def: Definição do problema
            repetitions: Repetições por amostra

        Returns:
            Dict[str, List[float]]: Resultados por métrica
        """
        results = {"distance": [], "execution_time": [], "fitness_calls": []}

        for i, sample in enumerate(samples):
            if i % 100 == 0:
                self._logger.info("Executando amostra %d/%d", i + 1, len(samples))

            # Converter amostra para parâmetros
            params = self._sample_to_params(sample, problem_def)

            # Executar com repetições
            sample_results = {"distance": [], "execution_time": [], "fitness_calls": []}

            for rep in range(repetitions):
                try:
                    result = self._executor.execute_single(
                        algorithm_name, dataset, params, timeout=600
                    )

                    sample_results["distance"].append(result["max_distance"])
                    sample_results["execution_time"].append(result["execution_time"])
                    sample_results["fitness_calls"].append(
                        result.get("metadata", {}).get("fitness_calls", 0)
                    )

                except (RuntimeError, ValueError, TimeoutError) as e:
                    self._logger.warning("Erro na execução %d: %s", rep, e)
                    # Usar valores de penalidade para falhas
                    sample_results["distance"].append(float("inf"))
                    sample_results["execution_time"].append(600.0)  # timeout
                    sample_results["fitness_calls"].append(0)

            # Agregar repetições (média)
            for metric in results.keys():
                valid_values = [v for v in sample_results[metric] if not np.isinf(v)]
                if valid_values:
                    results[metric].append(np.mean(valid_values))
                else:
                    results[metric].append(float("inf"))

        return results

    def _sample_to_params(
        self, sample: np.ndarray, problem_def: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Converte amostra numérica para parâmetros do algoritmo.

        Args:
            sample: Amostra numérica
            problem_def: Definição do problema

        Returns:
            Dict[str, Any]: Parâmetros do algoritmo
        """
        params = {}
        for name, value in zip(problem_def["names"], sample):
            # Converter para tipo apropriado baseado no nome do parâmetro
            if name in [
                "pop_size",
                "max_gens",
                "min_blocks",
                "max_blocks",
                "l_div",
                "min_block_size",
                "block_size",
                "k_candidates",
                "local_iters",
                "max_d",
                "warn_threshold",
                "rediv_freq",
                "immigrant_freq",
            ]:
                # Parâmetros inteiros
                params[name] = int(round(value))
            elif name in [
                "cross_prob",
                "mut_prob",
                "elite_rate",
                "initial_blocks",
                "immigrant_ratio",
                "no_improve_patience",
                "d_factor",
                "noise",
            ]:
                # Parâmetros float
                params[name] = float(value)
            elif name in [
                "crossover_type",
                "mutation_type",
                "selection_type",
                "refinement_type",
                "tie_break",
            ]:
                # Parâmetros categóricos - mapear índice para valor real
                # Definir mapeamentos conhecidos
                categorical_maps = {
                    "crossover_type": ["one_point", "uniform", "blend_blocks"],
                    "mutation_type": ["multi", "inversion", "transposition"],
                    "selection_type": ["tournament", "roulette", "ranking"],
                    "refinement_type": ["greedy", "swap", "insertion", "2opt"],
                    "tie_break": ["lex", "random", "first"],
                }
                if name in categorical_maps:
                    idx = int(round(value)) % len(categorical_maps[name])
                    params[name] = categorical_maps[name][idx]
                else:
                    params[name] = int(round(value))
            else:
                # Default para float
                params[name] = float(value)

        return params

    def _analyze_sensitivity(
        self,
        problem_def: Dict[str, Any],
        analysis_method: str,
        results: Dict[str, List[float]],
        output_metrics: List[str],
        sensitivity_config: Dict[str, Any],
        samples: Optional[np.ndarray] = None,
    ) -> Dict[str, Dict[str, Any]]:
        """
        Executa análise de sensibilidade nos resultados.

        Args:
            problem_def: Definição do problema
            analysis_method: Método de análise
            results: Resultados das execuções
            output_metrics: Métricas a analisar
            sensitivity_config: Configuração completa

        Returns:
            Dict[str, Dict[str, Any]]: Índices de sensibilidade
        """
        sensitivity_indices = {}

        for metric in output_metrics:
            if metric not in results:
                continue

            Y = np.array(results[metric])
            # Remover valores infinitos
            valid_mask = np.isfinite(Y)
            if not np.any(valid_mask):
                continue

            Y_clean = Y[valid_mask]

            try:
                if analysis_method == "morris":
                    # Morris precisa das amostras X e outputs Y
                    if samples is not None and len(samples) == len(Y):
                        X_clean = samples[valid_mask]
                        Si = morris.analyze(
                            problem_def,
                            X_clean,
                            Y_clean,
                            conf_level=0.95,
                            print_to_console=False,
                            num_levels=4,
                        )
                        sensitivity_indices[metric] = {
                            "mu": Si["mu"].tolist(),
                            "mu_star": Si["mu_star"].tolist(),
                            "sigma": Si["sigma"].tolist(),
                            "mu_star_conf": Si["mu_star_conf"].tolist(),
                            "parameter_names": problem_def["names"],
                        }
                    else:
                        self._logger.warning(
                            "Morris analysis skipped for %s: samples not available",
                            metric,
                        )
                        sensitivity_indices[metric] = {
                            "error": "Samples not available for Morris analysis"
                        }

                elif analysis_method == "sobol":
                    sobol_config = sensitivity_config.get("method_config", {})
                    Si = sobol.analyze(
                        problem_def,
                        Y_clean,
                        calc_second_order=sobol_config.get("calc_second_order", True),
                        conf_level=sobol_config.get("conf_level", 0.95),
                        print_to_console=False,
                    )
                    sensitivity_indices[metric] = {
                        "S1": Si["S1"].tolist(),
                        "S1_conf": Si["S1_conf"].tolist(),
                        "ST": Si["ST"].tolist(),
                        "ST_conf": Si["ST_conf"].tolist(),
                        "parameter_names": problem_def["names"],
                    }
                    if sobol_config.get("calc_second_order", True) and "S2" in Si:
                        sensitivity_indices[metric]["S2"] = Si["S2"].tolist()
                        sensitivity_indices[metric]["S2_conf"] = Si["S2_conf"].tolist()

                elif analysis_method == "fast":
                    Si = fast.analyze(problem_def, Y_clean, print_to_console=False)
                    sensitivity_indices[metric] = {
                        "S1": Si["S1"].tolist(),
                        "ST": Si["ST"].tolist(),
                        "parameter_names": problem_def["names"],
                    }

                elif analysis_method == "delta":
                    # Delta method precisa de samples originais
                    if samples is not None and len(samples) == len(Y):
                        X_clean = samples[valid_mask]
                        Si = delta.analyze(
                            problem_def, X_clean, Y_clean, print_to_console=False
                        )
                        sensitivity_indices[metric] = {
                            "delta": Si["delta"].tolist(),
                            "delta_conf": Si["delta_conf"].tolist(),
                            "S1": Si["S1"].tolist(),
                            "S1_conf": Si["S1_conf"].tolist(),
                            "parameter_names": problem_def["names"],
                        }
                    else:
                        self._logger.warning(
                            "Delta analysis skipped for %s: samples not available",
                            metric,
                        )
                        sensitivity_indices[metric] = {
                            "error": "Samples not available for delta analysis"
                        }

            except (ValueError, RuntimeError) as e:
                self._logger.warning(
                    "Erro na análise %s para %s: %s", analysis_method, metric, e
                )
                sensitivity_indices[metric] = {"error": str(e)}

        return sensitivity_indices
