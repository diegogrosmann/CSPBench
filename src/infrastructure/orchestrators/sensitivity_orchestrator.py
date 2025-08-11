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

    def __init__(self, algorithm_registry, executor, monitoring_service=None):
        """
        Inicializa orchestrator de sensibilidade.

        Args:
            algorithm_registry: Registry de algoritmos
            executor: Executor para execução de algoritmos
            monitoring_service: Serviço de monitoramento (opcional)
        """
        self._algorithm_registry = algorithm_registry
        self._executor = executor
        self._monitoring_service = monitoring_service
        self._logger = logging.getLogger(__name__)

    def execute_sensitivity_analysis(
        self,
        algorithm_name: str,
        dataset: Dataset,
        sensitivity_config: Dict[str, Any],
        *,
        task_index: int = 1,
        total_tasks: int = 1,
        dataset_index: int = 1,
        total_datasets: int = 1,
        config_index: int = 1,
        total_configs: int = 1,
        algorithm_index: int = 1,
        total_algorithms: int = 1,
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

            # Saída inicial padronizada (similar ao estilo de execução) e eventos de monitoramento
            dataset_name = getattr(dataset, "name", getattr(dataset, "metadata", {}).get("name", "Dataset"))
            config_name = sensitivity_config.get("config_name", "Default Configuration")
            if self._monitoring_service:
                from src.application.monitoring.progress_events import TaskType, UnifiedPhase, DisplayEvent
                # Iniciar task de análise se ainda não ativa
                try:
                    self._monitoring_service.start_task(
                        TaskType.SENSITIVITY,
                        f"Sensitivity Analysis - {algorithm_name}",
                        {"dataset": dataset_name, "algorithm": algorithm_name},
                    )
                except Exception:
                    pass
                # Evento inicial unificado, com todos os índices hierárquicos
                try:
                    self._monitoring_service.emit_event(
                        DisplayEvent(
                            phase=UnifiedPhase.ANALYSIS,
                            message=f"Starting sensitivity: {algorithm_name} on {dataset_name}",
                            dataset_id=dataset_name,
                            algorithm_name=algorithm_name,
                            progress=0.0,
                            payload={
                                "method": sensitivity_config.get("analysis_method", "morris"),
                                "task_index": task_index,
                                "total_tasks": total_tasks,
                                "task_name": "Sensitivity Analysis",
                                "dataset_index": dataset_index,
                                "total_datasets": total_datasets,
                                "config_index": config_index,
                                "total_configs": total_configs,
                                "algorithm_index": algorithm_index,
                                "total_algorithms": total_algorithms,
                            },
                        )
                    )
                except Exception:
                    pass
            else:
                print("\n🔄 Sensitivity Analysis")
                print(f"  📊 Dataset: {dataset_name}")
                print(f"  ⚙️  Algorithm: {algorithm_name}")
                print(f"  🧪 Method: {sensitivity_config.get('analysis_method', 'morris')}")
                print("  ------------------------------------------------------------")

            self._logger.debug("Starting sensitivity analysis: %s | dataset=%s", algorithm_name, dataset_name)
            self._logger.debug("Full sensitivity_config=%s", sensitivity_config)

            # Extrair configurações
            analysis_method = sensitivity_config.get("analysis_method", "morris")
            parameters = sensitivity_config["parameters"]
            self._logger.debug("parameters=%s (type=%s)", parameters, type(parameters))

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
                algorithm_name,
                dataset,
                samples,
                problem_definition,
                repetitions,
                dataset_name=dataset_name,
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
                "total_time": execution_time,  # alias usado em summaries
                "avg_time_per_sample": execution_time / max(1, len(samples)),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            }

            if self._monitoring_service:
                try:
                    from src.application.monitoring.progress_events import UnifiedPhase, DisplayEvent
                    self._monitoring_service.emit_event(
                        DisplayEvent(
                            phase=UnifiedPhase.ANALYSIS,
                            message=f"Completed sensitivity: {algorithm_name}",
                            dataset_id=dataset_name,
                            algorithm_name=algorithm_name,
                            progress=1.0,
                            payload={
                                "n_samples": len(samples),
                                "parameters": list(parameters.keys()),
                                "method": analysis_method,
                                "execution_time": execution_time,
                                "task_index": task_index,
                                "total_tasks": total_tasks,
                                "task_name": "Sensitivity Analysis",
                                "dataset_index": dataset_index,
                                "total_datasets": total_datasets,
                                "config_index": config_index,
                                "total_configs": total_configs,
                                "algorithm_index": algorithm_index,
                                "total_algorithms": total_algorithms,
                            },
                        )
                    )
                except Exception:
                    pass

            self._logger.info(
                "Análise de sensibilidade concluída em %.2fs", execution_time
            )
            return final_results

        except Exception as e:
            import traceback
            self._logger.error(
                "Erro na análise de sensibilidade: %s\n%s", e, traceback.format_exc()
            )
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
        self._logger.debug(
            "Creating SALib problem from parameters (%d)", len(parameters)
        )

        problem_def = {"num_vars": len(parameters), "names": [], "bounds": []}

        for param_name, param_config in parameters.items():
            self._logger.debug(
                "Processing parameter %s -> %s", param_name, param_config
            )

            problem_def["names"].append(param_name)

            param_type = param_config.get("type", "float")
            self._logger.debug("param_type=%s", param_type)

            if param_type in ["integer", "float"]:
                bounds = param_config["bounds"]
                self._logger.debug("bounds=%s", bounds)
                problem_def["bounds"].append(bounds)
            elif param_type == "categorical":
                # Para categóricos, mapear para índices
                values = param_config["values"]
                self._logger.debug("categorical values=%s", values)
                problem_def["bounds"].append([0, len(values) - 1])
            else:
                self._logger.debug(
                    "Unknown param_type %s, defaulting to float bounds", param_type
                )
                bounds = param_config.get("bounds", [0.0, 1.0])
                problem_def["bounds"].append(bounds)

        self._logger.debug("Final problem_def=%s", problem_def)
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
        dataset_name: str = "",
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
            # Barra de progresso simples (não usa overwrite agressivo para evitar conflito com display principal)
            progress = (i + 1) / len(samples) * 100
            if (i % max(1, len(samples)//50) == 0) or (i + 1 == len(samples)):
                bar_width = 30
                filled = int(progress / 100 * bar_width)
                bar = "█" * filled + "░" * (bar_width - filled)
                if not self._monitoring_service:
                    print(
                        f"    🔬 Samples: [{bar}] {progress:5.1f}% ({i+1}/{len(samples)})",
                        end="\r" if i + 1 < len(samples) else "\n",
                        flush=True,
                    )
                else:
                    try:
                        from src.application.monitoring.progress_events import UnifiedPhase, DisplayEvent
                        self._monitoring_service.emit_event(
                            DisplayEvent(
                                phase=UnifiedPhase.ANALYSIS,
                                message="sampling",
                                dataset_id=dataset_name,
                                algorithm_name=algorithm_name,
                                progress=(i + 1) / len(samples),
                                payload={
                                    "current_sample": i + 1,
                                    "total_samples": len(samples),
                                    "repetitions": repetitions,
                                },
                            )
                        )
                    except Exception:
                        pass
            if i % 100 == 0:
                self._logger.info("Executando amostra %d/%d", i + 1, len(samples))

            # Converter amostra para parâmetros
            params = self._sample_to_params(sample, problem_def)

            # Executar com repetições (potencialmente paralelas)
            try:
                # Usar o método público de repetições paralelas do executor
                sample_results_list = []

                for rep in range(repetitions):
                    try:
                        result = self._executor.execute_single(
                            algorithm_name, dataset, params, timeout=600
                        )
                        sample_results_list.append(result)
                    except (RuntimeError, ValueError, TimeoutError) as e:
                        self._logger.warning("Erro na execução %d: %s", rep, e)
                        sample_results_list.append(None)

                # Processar resultados das repetições
                sample_results = {
                    "distance": [],
                    "execution_time": [],
                    "fitness_calls": [],
                }

                for result in sample_results_list:
                    if result is not None:
                        sample_results["distance"].append(result["max_distance"])
                        sample_results["execution_time"].append(
                            result["execution_time"]
                        )
                        sample_results["fitness_calls"].append(
                            result.get("metadata", {}).get("fitness_calls", 0)
                        )
                    else:
                        # Resultado com falha
                        sample_results["distance"].append(float("inf"))
                        sample_results["execution_time"].append(600.0)  # timeout
                        sample_results["fitness_calls"].append(0)

                # Agregar repetições (média)
                for metric in results.keys():
                    valid_values = [
                        v for v in sample_results[metric] if not np.isinf(v)
                    ]
                    if valid_values:
                        results[metric].append(np.mean(valid_values))
                    else:
                        results[metric].append(float("inf"))

            except Exception as e:
                self._logger.warning("Erro na execução da amostra %d: %s", i, e)
                # Usar valores de penalidade para falhas completas da amostra
                for metric in results.keys():
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
        self._logger.debug(
            "Converting sample to params names=%s values=%s", problem_def["names"], sample
        )

        for name, value in zip(problem_def["names"], sample):
            self._logger.debug("Processing parameter %s=%s", name, value)
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
                "crossover_method",  # Add crossover_method here
                "crossover_type",
                "mutation_type",
                "selection_type",
                "refinement_type",
                "tie_break",
            ]:
                # Parâmetros categóricos - mapear índice para valor real
                # Definir mapeamentos conhecidos
                categorical_maps = {
                    "crossover_method": [
                        "one_point",
                        "uniform",
                        "blend_blocks",
                    ],  # Add crossover_method mapping
                    "crossover_type": ["one_point", "uniform", "blend_blocks"],
                    "mutation_type": ["multi", "inversion", "transposition"],
                    "selection_type": ["tournament", "roulette", "ranking"],
                    "refinement_type": ["greedy", "swap", "insertion", "2opt"],
                    "tie_break": ["lex", "random", "first"],
                }
                if name in categorical_maps:
                    idx = int(round(value)) % len(categorical_maps[name])
                    params[name] = categorical_maps[name][idx]
                    self._logger.debug(
                        "Categorical %s idx=%d -> %s", name, idx, params[name]
                    )
                else:
                    params[name] = int(round(value))
                    self._logger.debug(
                        "Unknown categorical %s treated as int -> %s", name, params[name]
                    )
            else:
                # Default para float
                params[name] = float(value)
                self._logger.debug("Default float %s -> %s", name, params[name])

        self._logger.debug("Final params=%s", params)
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
