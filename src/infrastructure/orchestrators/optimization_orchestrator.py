"""
Optimization Orchestrator - Orquestração de Otimização com Optuna

Coordena a execução de otimização de hiperparâmetros usando Optuna,
incluindo salvamento incremental, relatórios avançados e sistema de recovery.
"""

import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import optuna
from optuna.pruners import MedianPruner, SuccessiveHalvingPruner
from optuna.samplers import CmaEsSampler, RandomSampler, TPESampler

from src.application.ports import AlgorithmRegistry, DatasetRepository
from src.domain import CSPAlgorithm, Dataset
from src.infrastructure.logging_config import get_logger
from src.infrastructure.orchestrators.optimization_report_generator import (
    OptimizationReportGenerator,
)


class OptimizationOrchestrator:
    """Orquestrador para otimização de hiperparâmetros com Optuna."""

    def __init__(
        self,
        algorithm_registry: AlgorithmRegistry,
        dataset_repository: DatasetRepository,
        config: Dict[str, Any],
    ):
        """
        Inicializa o orquestrador.

        Args:
            algorithm_registry: Registry de algoritmos
            dataset_repository: Repositório de datasets
            config: Configuração completa do batch
        """
        self.algorithm_registry = algorithm_registry
        self.dataset_repository = dataset_repository
        self.config = config
        self.logger = get_logger(__name__)

        # Configuração de otimização
        self.optimization_config = config.get("optimization", {})
        self.export_config = config.get("export", {})
        self.monitoring_config = config.get("monitoring", {})
        self.resources_config = config.get("resources", {})

        # Configuração do estudo
        self.study_name = self.optimization_config.get(
            "study_name", "optimization_study"
        )
        self.direction = self.optimization_config.get("direction", "minimize")
        self.n_trials = int(
            os.getenv("N_TRIALS", self.optimization_config.get("n_trials", 100))
        )
        self.timeout_per_trial = self.optimization_config.get("timeout_per_trial", 300)

        # Configuração de salvamento usando SessionManager
        from pathlib import Path

        import yaml

        from src.infrastructure.session_manager import SessionManager

        # Carregar configuração principal
        config_path = Path("config/settings.yaml")
        if config_path.exists():
            with open(config_path, "r", encoding="utf-8") as f:
                settings = yaml.safe_load(f)
        else:
            settings = {}

        session_manager = SessionManager(settings)

        # Se não existe sessão, criar uma nova
        if not session_manager.get_session_folder():
            session_manager.create_session()

        # Usar base_result_dir + session_folder para otimização
        base_result_dir = (
            settings.get("infrastructure", {})
            .get("result", {})
            .get("base_result_dir", "./outputs/results")
        )
        session_folder = session_manager.get_session_folder()

        if session_folder:
            self.destination = os.path.join(base_result_dir, session_folder)
        else:
            self.destination = os.path.join(base_result_dir, "default_session")

        # Override se especificado na configuração
        if "destination" in self.export_config:
            self.destination = self.export_config["destination"]

        self.partial_results_file = Path(self.destination) / "partial_results.json"
        self.checkpoint_interval = self.monitoring_config.get("checkpointing", {}).get(
            "interval", 10
        )

        # Resultados e estado
        self.partial_results = []
        self.best_params = None
        self.best_value = None
        self.trial_count = 0
        self.start_time = None

        # Criar diretório de destino
        Path(self.destination).mkdir(parents=True, exist_ok=True)

        # Gerador de relatórios
        self.report_generator = OptimizationReportGenerator(
            destination=self.destination, config=config
        )

    def run_optimization(self) -> Dict[str, Any]:
        """
        Executa otimização completa.

        Returns:
            Dict[str, Any]: Resultados da otimização
        """
        try:
            self.logger.info("Iniciando otimização com Optuna")
            self.start_time = time.time()

            # Configurar Optuna
            study = self._create_study()

            # Carregar dataset
            dataset = self._load_dataset()

            # Carregar algoritmo
            algorithm_class = self._load_algorithm()

            # Definir função objetivo
            def objective(trial: optuna.trial.Trial) -> float:
                return self._objective_function(trial, algorithm_class, dataset)

            # Configurar callbacks
            callbacks = self._setup_callbacks()

            # Executar otimização
            study.optimize(
                objective,
                n_trials=self.n_trials,
                timeout=self.timeout_per_trial * self.n_trials,
                callbacks=callbacks,
                n_jobs=self.resources_config.get("parallel", {}).get("n_jobs", 1),
                show_progress_bar=self.monitoring_config.get("progress", {}).get(
                    "show_progress_bar", True
                ),
            )

            # Processar resultados finais
            results = self._process_final_results(study)

            # Gerar relatórios
            if self.export_config.get("enabled", True):
                self._generate_reports(study, results)

            self.logger.info(f"Otimização concluída: {len(study.trials)} trials")
            return results

        except Exception as e:
            self.logger.error(f"Erro durante otimização: {e}")
            raise

    def _create_study(self) -> optuna.Study:
        """Cria estudo Optuna com configurações específicas."""
        # Configurar sampler
        sampler_config = self.optimization_config.get("optuna_config", {})
        sampler_name = sampler_config.get("sampler", "TPESampler")

        if sampler_name == "TPESampler":
            sampler = TPESampler(
                n_startup_trials=sampler_config.get("n_startup_trials", 10),
                n_ei_candidates=sampler_config.get("n_ei_candidates", 24),
                multivariate=sampler_config.get("multivariate", False),
            )
        elif sampler_name == "RandomSampler":
            sampler = RandomSampler(seed=sampler_config.get("seed", None))
        elif sampler_name == "CmaEsSampler":
            sampler = CmaEsSampler(
                n_startup_trials=sampler_config.get("n_startup_trials", 1),
                restart_strategy=sampler_config.get("restart_strategy", "ipop"),
            )
        else:
            sampler = TPESampler()

        # Configurar pruner
        pruner_name = sampler_config.get("pruner", "MedianPruner")

        if pruner_name == "MedianPruner":
            pruner = MedianPruner(
                n_startup_trials=sampler_config.get("n_startup_trials", 5),
                n_warmup_steps=sampler_config.get("n_warmup_steps", 10),
                interval_steps=sampler_config.get("interval_steps", 5),
            )
        elif pruner_name == "SuccessiveHalvingPruner":
            pruner = SuccessiveHalvingPruner(
                min_resource=sampler_config.get("min_resource", 1),
                reduction_factor=sampler_config.get("reduction_factor", 4),
                min_early_stopping_rate=sampler_config.get(
                    "min_early_stopping_rate", 0
                ),
            )
        else:
            pruner = MedianPruner()

        # Configurar storage se especificado
        storage_url = os.getenv("OPTUNA_STORAGE", sampler_config.get("storage", None))

        # Criar estudo
        study = optuna.create_study(
            study_name=os.getenv("OPTUNA_STUDY_NAME", self.study_name),
            direction=self.direction,
            sampler=sampler,
            pruner=pruner,
            storage=storage_url,
            load_if_exists=self.monitoring_config.get("checkpointing", {}).get(
                "recovery", True
            ),
        )

        self.logger.info(f"Estudo criado: {study.study_name}")
        return study

    def _load_dataset(self) -> Dataset:
        """Carrega dataset para otimização."""
        dataset_id = self.config.get("dataset")
        if not dataset_id:
            raise ValueError("Dataset não especificado")

        if not self.dataset_repository.exists(dataset_id):
            raise ValueError(f"Dataset não encontrado: {dataset_id}")

        return self.dataset_repository.load(dataset_id)

    def _load_algorithm(self) -> type[CSPAlgorithm]:
        """Carrega classe do algoritmo."""
        algorithm_name = self.config.get("algorithm")
        if not algorithm_name:
            raise ValueError("Algoritmo não especificado")

        if not self.algorithm_registry.algorithm_exists(algorithm_name):
            raise ValueError(f"Algoritmo não encontrado: {algorithm_name}")

        return self.algorithm_registry.get_algorithm(algorithm_name)

    def _objective_function(
        self,
        trial: optuna.trial.Trial,
        algorithm_class: type[CSPAlgorithm],
        dataset: Dataset,
    ) -> float:
        """
        Função objetivo para otimização.

        Args:
            trial: Trial do Optuna
            algorithm_class: Classe do algoritmo
            dataset: Dataset para teste

        Returns:
            float: Valor da função objetivo
        """
        try:
            # Gerar parâmetros baseado na configuração
            params = self._generate_trial_params(trial)

            # Executar algoritmo usando a interface correta
            # A interface CSPAlgorithm espera: __init__(strings, alphabet, **params) e run()
            algorithm_instance = algorithm_class(
                strings=dataset.sequences, alphabet=dataset.alphabet, **params
            )

            best_string, max_distance, metadata = algorithm_instance.run()

            # Salvar resultado parcial
            trial_result = {
                "trial_number": trial.number,
                "params": params,
                "max_distance": max_distance,
                "best_string": best_string,
                "metadata": metadata,
                "timestamp": time.time(),
                "trial_id": trial.number,
                "value": max_distance,
                "state": "COMPLETE",
            }

            self.partial_results.append(trial_result)

            # Atualizar melhores resultados
            if (
                self.best_value is None
                or (self.direction == "minimize" and max_distance < self.best_value)
                or (self.direction == "maximize" and max_distance > self.best_value)
            ):
                self.best_value = max_distance
                self.best_params = params.copy()

            # Salvar progresso incremental
            self.trial_count += 1
            if self.trial_count % self.checkpoint_interval == 0:
                self._save_partial_results()

            return max_distance

        except Exception as e:
            self.logger.error(f"Erro no trial {trial.number}: {e}")

            # Registrar erro
            trial_result = {
                "trial_number": trial.number,
                "params": self._generate_trial_params(trial),
                "error": str(e),
                "timestamp": time.time(),
                "trial_id": trial.number,
                "state": "FAIL",
            }
            self.partial_results.append(trial_result)

            # Retornar valor de penalidade
            return float("inf") if self.direction == "minimize" else float("-inf")

    def _generate_trial_params(self, trial: optuna.trial.Trial) -> Dict[str, Any]:
        """Gera parâmetros do trial baseado na configuração."""
        params = {}
        parameters_config = self.optimization_config.get("parameters", {})

        for param_name, param_config in parameters_config.items():
            param_type = param_config.get("type", "uniform")

            if param_type == "categorical":
                params[param_name] = trial.suggest_categorical(
                    param_name, param_config["choices"]
                )
            elif param_type == "int":
                params[param_name] = trial.suggest_int(
                    param_name,
                    param_config["low"],
                    param_config["high"],
                    step=param_config.get("step", 1),
                )
            elif param_type == "float":
                params[param_name] = trial.suggest_float(
                    param_name,
                    param_config["low"],
                    param_config["high"],
                    step=param_config.get("step", None),
                )
            elif param_type == "uniform":
                params[param_name] = trial.suggest_uniform(
                    param_name, param_config["low"], param_config["high"]
                )
            elif param_type == "loguniform":
                params[param_name] = trial.suggest_loguniform(
                    param_name, param_config["low"], param_config["high"]
                )
            else:
                self.logger.warning(f"Tipo de parâmetro desconhecido: {param_type}")

        return params

    def _setup_callbacks(self) -> List:
        """Configura callbacks do Optuna."""
        callbacks = []

        # Callback de progresso
        progress_config = self.monitoring_config.get("progress", {})
        if progress_config.get("log_interval", 0) > 0:

            def progress_callback(study, trial):
                if trial.number % progress_config["log_interval"] == 0:
                    self.logger.info(f"Trial {trial.number}: {trial.value}")

            callbacks.append(progress_callback)

        # Callback de early stopping
        early_stopping_config = progress_config.get("early_stopping", {})
        if early_stopping_config.get("enabled", False):

            def early_stopping_callback(study, trial):
                if len(study.trials) >= early_stopping_config["patience"]:
                    recent_trials = study.trials[-early_stopping_config["patience"] :]
                    if self.direction == "minimize":
                        recent_best = min(
                            t.value for t in recent_trials if t.value is not None
                        )
                        if (
                            study.best_value is not None
                            and recent_best - study.best_value
                            > early_stopping_config.get("min_improvement", 0.001)
                        ):
                            study.stop()
                    else:
                        recent_best = max(
                            t.value for t in recent_trials if t.value is not None
                        )
                        if (
                            study.best_value is not None
                            and study.best_value - recent_best
                            > early_stopping_config.get("min_improvement", 0.001)
                        ):
                            study.stop()

            callbacks.append(early_stopping_callback)

        return callbacks

    def _save_partial_results(self):
        """Salva resultados parciais em arquivo JSON."""
        try:
            partial_data = {
                "study_name": self.study_name,
                "direction": self.direction,
                "n_trials_completed": self.trial_count,
                "n_trials_total": self.n_trials,
                "best_value": self.best_value,
                "best_params": self.best_params,
                "start_time": self.start_time,
                "last_update": time.time(),
                "trials": self.partial_results,
            }

            with open(self.partial_results_file, "w", encoding="utf-8") as f:
                json.dump(partial_data, f, indent=2, ensure_ascii=False)

            self.logger.debug(
                f"Resultados parciais salvos: {self.partial_results_file}"
            )

        except Exception as e:
            self.logger.error(f"Erro ao salvar resultados parciais: {e}")

    def _process_final_results(self, study: optuna.Study) -> Dict[str, Any]:
        """Processa resultados finais da otimização."""
        end_time = time.time()
        total_time = end_time - (self.start_time or end_time)

        results = {
            "study_name": study.study_name,
            "direction": study.direction.name,
            "n_trials": len(study.trials),
            "best_value": study.best_value,
            "best_params": study.best_params,
            "best_trial": study.best_trial.number if study.best_trial else None,
            "total_time": total_time,
            "start_time": self.start_time,
            "end_time": end_time,
            "algorithm": self.config.get("algorithm"),
            "dataset": self.config.get("dataset"),
            "trials": [
                {
                    "number": trial.number,
                    "value": trial.value,
                    "params": trial.params,
                    "state": trial.state.name,
                    "datetime_start": (
                        trial.datetime_start.isoformat()
                        if trial.datetime_start
                        else None
                    ),
                    "datetime_complete": (
                        trial.datetime_complete.isoformat()
                        if trial.datetime_complete
                        else None
                    ),
                    "duration": (
                        trial.duration.total_seconds() if trial.duration else None
                    ),
                    "user_attrs": trial.user_attrs,
                    "system_attrs": trial.system_attrs,
                }
                for trial in study.trials
            ],
        }

        # Salvar resultados finais
        final_results_file = Path(self.destination) / "final_results.json"
        with open(final_results_file, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

        self.logger.info(f"Resultados finais salvos: {final_results_file}")
        return results

    def _generate_reports(self, study: optuna.Study, results: Dict[str, Any]):
        """Gera relatórios avançados da otimização."""
        try:
            self.logger.info("Gerando relatórios de otimização")

            # Gerar todos os tipos de relatórios
            self.report_generator.generate_all_reports(study, results)

            self.logger.info("Relatórios gerados com sucesso")

        except Exception as e:
            self.logger.error(f"Erro ao gerar relatórios: {e}")
