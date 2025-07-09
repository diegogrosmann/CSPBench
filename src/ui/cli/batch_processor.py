"""
Processador unificado de batch para execução, otimização e análise de sensibilidade.

Este módulo implementa a lógica unificada para processar arquivos de configuração
YAML padronizados e executar diferentes tipos de tarefas.
"""

import logging
import os
import time
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from src.ui.cli.batch_config_extractor import BatchConfigExtractor
from src.ui.cli.console_manager import console

logger = logging.getLogger(__name__)


class UnifiedBatchProcessor:
    """Processador unificado para todas as tarefas de batch."""

    def __init__(self, config_path: str, silent: bool = False):
        """
        Inicializa o processador.

        Args:
            config_path: Caminho para o arquivo de configuração
            silent: Se True, executa em modo silencioso
        """
        self.config_path = config_path
        self.silent = silent
        self.extractor = BatchConfigExtractor(config_path)
        self.console_print = (
            console.print if not silent else lambda *args, **kwargs: None
        )

    def process(self) -> Dict[str, Any]:
        """
        Processa a configuração de batch baseado no tipo de tarefa.

        Returns:
            Resultados da execução
        """
        batch_info = self.extractor.get_batch_info()
        task_type = self.extractor.get_task_type()

        if not self.silent:
            self.console_print(
                f"🚀 Iniciando batch: {batch_info.get('nome', 'Sem nome')}"
            )
            self.console_print(f"📄 Descrição: {batch_info.get('descricao', 'N/A')}")
            self.console_print(f"🔧 Tipo de tarefa: {task_type}")
            self.console_print(
                f"⏰ Timeout global: {self.extractor.get_global_timeout()}s"
            )

        # Processar baseado no tipo de tarefa
        if task_type == "execution":
            return self._process_execution()
        elif task_type == "sensitivity":
            return self._process_sensitivity()
        elif task_type == "optimization":
            return self._process_optimization()
        else:
            raise ValueError(f"Tipo de tarefa não suportado: {task_type}")

    def _process_execution(self) -> Dict[str, Any]:
        """Processa execuções de algoritmos."""
        execution_configs = self.extractor.get_execution_configs()
        algorithms = self.extractor.get_algorithms()

        self.console_print(f"🧮 Algoritmos: {algorithms}")
        self.console_print(f"📊 Execuções: {len(execution_configs)}")

        all_results = {}

        for exec_idx, exec_config in enumerate(execution_configs, 1):
            exec_name = exec_config.get("nome", f"Execução {exec_idx}")

            if not self.silent:
                self.console_print(f"\n{'='*60}")
                self.console_print(
                    f"🚀 Execução {exec_idx}/{len(execution_configs)}: {exec_name}"
                )
                self.console_print(f"{'='*60}")

            try:
                # Extrair configurações da execução
                dataset_id = exec_config["dataset"]
                num_execs = exec_config.get("runs_per_algorithm_per_base", 1)
                num_bases = exec_config.get("num_bases", 1)
                timeout = exec_config.get("timeout", 300)

                # Gerar dataset
                seqs, dataset_params = self.extractor.generate_dataset(
                    dataset_id, silent=self.silent
                )
                alphabet = "".join(sorted(set("".join(seqs))))

                if not self.silent:
                    self.console_print(
                        f"📊 Dataset: n={len(seqs)}, L={len(seqs[0])}, |Σ|={len(alphabet)}"
                    )
                    self.console_print(f"⏰ Timeout: {timeout}s")

                # Executar algoritmos
                results = self._execute_algorithms(
                    algorithms=algorithms,
                    seqs=seqs,
                    alphabet=alphabet,
                    num_execs=num_execs,
                    timeout=timeout,
                    dataset_params=dataset_params,
                )

                all_results[exec_name] = {
                    "config": exec_config,
                    "dataset_params": dataset_params,
                    "results": results,
                }

                if not self.silent:
                    self.console_print(f"✅ Execução {exec_idx} concluída com sucesso!")

            except Exception as e:
                logger.exception(f"Erro na execução {exec_idx}: {e}")
                if not self.silent:
                    self.console_print(f"❌ Erro na execução {exec_idx}: {e}")

                all_results[exec_name] = {"config": exec_config, "error": str(e)}

        return all_results

    def _process_sensitivity(self) -> Dict[str, Any]:
        """Processa análises de sensibilidade."""
        sensitivity_configs = self.extractor.get_sensitivity_configs()

        self.console_print(f"🔬 Análises de sensibilidade: {len(sensitivity_configs)}")

        all_results = {}

        for analysis_idx, analysis_config in enumerate(sensitivity_configs, 1):
            analysis_name = analysis_config.get("nome", f"Análise {analysis_idx}")

            if not self.silent:
                self.console_print(f"\n{'='*60}")
                self.console_print(
                    f"🔬 Análise {analysis_idx}/{len(sensitivity_configs)}: {analysis_name}"
                )
                self.console_print(f"{'='*60}")

            try:
                # Extrair configurações da análise
                dataset_ids = analysis_config["datasets"]
                n_samples = analysis_config.get("n_samples", 100)
                timeout_per_sample = analysis_config.get("timeout_per_sample", 60)
                method = analysis_config.get("method", "morris")
                param_space = analysis_config.get("param_space", {})

                # Processar cada dataset
                dataset_results = {}
                for dataset_id in dataset_ids:
                    seqs, dataset_params = self.extractor.generate_dataset(
                        dataset_id, silent=self.silent
                    )
                    alphabet = "".join(sorted(set("".join(seqs))))

                    if not self.silent:
                        self.console_print(
                            f"📊 Dataset {dataset_id}: n={len(seqs)}, L={len(seqs[0])}"
                        )

                    # Executar análise para cada algoritmo no espaço de parâmetros
                    for algorithm_name, param_names in param_space.items():
                        if not self.silent:
                            self.console_print(
                                f"🔬 Analisando {algorithm_name} ({method})..."
                            )

                        # Executar análise de sensibilidade
                        result = self._run_sensitivity_analysis(
                            algorithm_name=algorithm_name,
                            seqs=seqs,
                            alphabet=alphabet,
                            param_names=param_names,
                            n_samples=n_samples,
                            timeout_per_sample=timeout_per_sample,
                            method=method,
                        )

                        dataset_results[f"{dataset_id}_{algorithm_name}"] = result

                all_results[analysis_name] = {
                    "config": analysis_config,
                    "results": dataset_results,
                }

                if not self.silent:
                    self.console_print(
                        f"✅ Análise {analysis_idx} concluída com sucesso!"
                    )

            except Exception as e:
                logger.exception(f"Erro na análise {analysis_idx}: {e}")
                if not self.silent:
                    self.console_print(f"❌ Erro na análise {analysis_idx}: {e}")

                all_results[analysis_name] = {
                    "config": analysis_config,
                    "error": str(e),
                }

        return all_results

    def _process_optimization(self) -> Dict[str, Any]:
        """Processa otimizações."""
        optimization_configs = self.extractor.get_optimization_configs()

        self.console_print(f"🚀 Otimizações: {len(optimization_configs)}")

        all_results = {}

        for opt_idx, opt_config in enumerate(optimization_configs, 1):
            opt_name = opt_config.get("nome", f"Otimização {opt_idx}")

            if not self.silent:
                self.console_print(f"\n{'='*60}")
                self.console_print(
                    f"🚀 Otimização {opt_idx}/{len(optimization_configs)}: {opt_name}"
                )
                self.console_print(f"{'='*60}")

            try:
                # Extrair configurações da otimização
                dataset_ids = opt_config["datasets"]
                n_trials = opt_config.get("n_trials", 100)
                timeout_per_trial = opt_config.get("timeout_per_trial", 60)
                param_space = opt_config.get("param_space", {})

                # Processar cada dataset
                dataset_results = {}
                for dataset_id in dataset_ids:
                    seqs, dataset_params = self.extractor.generate_dataset(
                        dataset_id, silent=self.silent
                    )
                    alphabet = "".join(sorted(set("".join(seqs))))

                    if not self.silent:
                        self.console_print(
                            f"📊 Dataset {dataset_id}: n={len(seqs)}, L={len(seqs[0])}"
                        )

                    # Executar otimização para cada algoritmo no espaço de parâmetros
                    for algorithm_name, param_config in param_space.items():
                        if not self.silent:
                            self.console_print(f"🚀 Otimizando {algorithm_name}...")

                        # Executar otimização
                        result = self._run_optimization(
                            algorithm_name=algorithm_name,
                            seqs=seqs,
                            alphabet=alphabet,
                            param_config=param_config,
                            n_trials=n_trials,
                            timeout_per_trial=timeout_per_trial,
                            opt_config=opt_config,
                        )

                        dataset_results[f"{dataset_id}_{algorithm_name}"] = result

                all_results[opt_name] = {
                    "config": opt_config,
                    "results": dataset_results,
                }

                if not self.silent:
                    self.console_print(
                        f"✅ Otimização {opt_idx} concluída com sucesso!"
                    )

            except Exception as e:
                logger.exception(f"Erro na otimização {opt_idx}: {e}")
                if not self.silent:
                    self.console_print(f"❌ Erro na otimização {opt_idx}: {e}")

                all_results[opt_name] = {"config": opt_config, "error": str(e)}

        return all_results

    def _execute_algorithms(
        self,
        algorithms: List[str],
        seqs: List[str],
        alphabet: str,
        num_execs: int,
        timeout: int,
        dataset_params: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Executa algoritmos usando interface curses ou modo tradicional."""
        use_curses = self.extractor.should_use_curses() and not self.silent

        if use_curses:
            try:
                from src.ui.curses_integration import CursesExecutionMonitor

                monitor = CursesExecutionMonitor(max_workers=4, timeout=timeout)
                results = monitor.execute_algorithms(
                    algorithm_names=algorithms,
                    seqs=seqs,
                    alphabet=alphabet,
                    num_execs=num_execs,
                    dataset_params=dataset_params,
                )
                return results

            except Exception as e:
                logger.warning(
                    f"Erro na interface curses: {e}. Usando modo tradicional."
                )
                use_curses = False

        if not use_curses:
            return self._execute_algorithms_traditional(
                algorithms, seqs, alphabet, num_execs, timeout, dataset_params
            )

        # Fallback caso nenhum caminho seja executado
        return {}

    def _execute_algorithms_traditional(
        self,
        algorithms: List[str],
        seqs: List[str],
        alphabet: str,
        num_execs: int,
        timeout: int,
        dataset_params: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Executa algoritmos em modo tradicional."""
        from algorithms.base import global_registry
        from src.core.interfaces import TaskStatus, create_executor

        results = {}
        executor = create_executor(timeout_seconds=timeout, max_workers=4)

        try:
            for alg_name in algorithms:
                if alg_name not in global_registry:
                    if not self.silent:
                        self.console_print(f"❌ Algoritmo {alg_name} não encontrado!")
                    continue

                AlgClass = global_registry[alg_name]
                is_deterministic = getattr(AlgClass, "is_deterministic", False)
                actual_num_execs = 1 if is_deterministic else num_execs

                if not self.silent:
                    if is_deterministic:
                        self.console_print(
                            f"  🔒 {alg_name} é determinístico - executando apenas 1 vez"
                        )
                    else:
                        self.console_print(
                            f"  🎲 {alg_name} é não-determinístico - executando {actual_num_execs} vezes"
                        )

                alg_results = []

                for i in range(actual_num_execs):
                    if not self.silent:
                        if actual_num_execs == 1:
                            self.console_print(f"  Executando {alg_name}")
                        else:
                            self.console_print(
                                f"  Executando {alg_name} - Run {i+1}/{actual_num_execs}"
                            )

                    instance = AlgClass(seqs, alphabet)
                    handle = executor.submit(instance)

                    # Aguardar conclusão
                    while executor.poll(handle) == TaskStatus.RUNNING:
                        time.sleep(0.1)

                    result = executor.result(handle)
                    alg_results.append(result)

                results[alg_name] = alg_results

        finally:
            if hasattr(executor, "shutdown"):
                executor.shutdown(wait=True)

        return results

    def _run_sensitivity_analysis(
        self,
        algorithm_name: str,
        seqs: List[str],
        alphabet: str,
        param_names: List[str],
        n_samples: int,
        timeout_per_sample: int,
        method: str,
    ) -> Dict[str, Any]:
        """Executa análise de sensibilidade."""
        try:
            from src.optimization.sensitivity_analyzer import (
                analyze_algorithm_sensitivity,
            )

            result = analyze_algorithm_sensitivity(
                algorithm_name=algorithm_name,
                sequences=seqs,
                alphabet=alphabet,
                n_samples=n_samples,
                timeout_per_sample=timeout_per_sample,
                show_progress=not self.silent,
                method=method,
            )

            return {
                "success": True,
                "method": result.method,
                "parameter_names": result.parameter_names,
                "first_order": result.first_order,
                "total_order": result.total_order,
                "mu_star": result.mu_star,
                "n_samples": result.n_samples,
                "analysis_time": result.analysis_time,
            }

        except Exception as e:
            logger.exception(f"Erro na análise de sensibilidade: {e}")
            return {"success": False, "error": str(e)}

    def _run_optimization(
        self,
        algorithm_name: str,
        seqs: List[str],
        alphabet: str,
        param_config: Dict[str, Any],
        n_trials: int,
        timeout_per_trial: int,
        opt_config: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Executa otimização."""
        try:
            from src.optimization.optuna_optimizer import optimize_algorithm

            result = optimize_algorithm(
                algorithm_name=algorithm_name,
                sequences=seqs,
                alphabet=alphabet,
                n_trials=n_trials,
                timeout_per_trial=timeout_per_trial,
                show_progress=not self.silent,
            )

            return {
                "success": True,
                "best_value": result.best_value,
                "best_params": result.best_params,
                "n_trials": result.n_trials,
                "study_name": result.study_name,
                "optimization_time": result.optimization_time,
            }

        except Exception as e:
            logger.exception(f"Erro na otimização: {e}")
            return {"success": False, "error": str(e)}
