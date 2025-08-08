"""
Fakes (Test Doubles) para testes dos serviços de aplicação.

Implementações simples das portas para uso em testes unitários.
"""

from typing import Any, Dict, List, Optional
from unittest.mock import Mock

from src.domain import CSPAlgorithm, Dataset, SyntheticDatasetGenerator


class MockAlgorithmClass(CSPAlgorithm):
    """Mock algorithm class for testing."""

    name = "MockAlgorithm"
    default_params = {}

    def run(self, dataset: Dataset, **kwargs) -> tuple[str, int, dict]:
        """Mock run method."""
        # Return a simple mock result
        return "ACGT" * (len(dataset.sequences[0]) // 4), 2, {"iterations": 10}

    def _report_progress(self, message: str) -> None:
        """Mock progress reporting."""
        pass


class FakeDatasetRepository:
    """Implementação fake do DatasetRepository para testes."""

    def __init__(self):
        self._datasets: Dict[str, Dataset] = {}
        # Adicionar alguns datasets de exemplo
        self._datasets["test_dataset"] = SyntheticDatasetGenerator.generate_random(
            n=10, length=20, alphabet="ACGT", seed=42
        )
        self._datasets["small_dataset"] = SyntheticDatasetGenerator.generate_random(
            n=5, length=10, alphabet="ACGT", seed=123
        )

    def save(self, dataset: Dataset, name: str) -> str:
        """Salva dataset com nome específico."""
        self._datasets[name] = dataset
        return name

    def load(self, identifier: str) -> Dataset:
        """Carrega dataset por identificador."""
        if identifier not in self._datasets:
            from src.domain import DatasetNotFoundError

            raise DatasetNotFoundError(f"Dataset não encontrado: {identifier}")
        return self._datasets[identifier]

    def list_available(self) -> List[str]:
        """Lista datasets disponíveis."""
        return list(self._datasets.keys())

    def exists(self, identifier: str) -> bool:
        """Verifica se dataset existe."""
        return identifier in self._datasets

    def delete(self, identifier: str) -> bool:
        """Remove dataset."""
        if identifier in self._datasets:
            del self._datasets[identifier]
            return True
        return False


class FakeAlgorithmRegistry:
    """Implementação fake do AlgorithmRegistry para testes."""

    def __init__(self):
        self._algorithms = {
            "baseline": "BaselineAlgorithm",
            "test_algo": "TestAlgorithm",
            "BLF-GA": "BLFGAAlgorithm",
            "CSC": "CSCAlgorithm",
        }

    def get_algorithm(self, name: str) -> type[CSPAlgorithm]:
        """Obtém classe de algoritmo por nome."""
        if name not in self._algorithms:
            from src.domain import AlgorithmNotFoundError

            raise AlgorithmNotFoundError(f"Algorithm not found: {name}")
        # Return a mock class that inherits from CSPAlgorithm
        return MockAlgorithmClass

    def list_algorithms(self) -> List[str]:
        """List available algorithms."""
        return list(self._algorithms.keys())

    def register_algorithm(self, algorithm_class) -> None:
        """Registra novo algoritmo."""
        name = getattr(algorithm_class, "name", algorithm_class.__name__)
        self._algorithms[name] = algorithm_class.__name__

    def algorithm_exists(self, name: str) -> bool:
        """Verifica se algoritmo existe."""
        return name in self._algorithms

    def get_algorithm_metadata(self, name: str) -> Dict[str, Any]:
        """Get algorithm metadata."""
        if name not in self._algorithms:
            from src.domain import AlgorithmNotFoundError

            raise AlgorithmNotFoundError(f"Algorithm not found: {name}")

        return {
            "name": name,
            "class": self._algorithms[name],
            "is_deterministic": True,
            "supports_parallel": False,
        }


class FakeExportPort:
    """Implementação fake do ExportPort para testes."""

    def __init__(self):
        self._exported_files = []

    def export(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """Export data in specified format."""
        file_path = f"{destination}.{format_type}"
        self._exported_files.append(
            {"path": file_path, "format": format_type, "data": results}
        )
        return file_path

    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """Exporta resultados em formato específico."""
        file_path = f"{destination}.{format_type}"
        self._exported_files.append(
            {"path": file_path, "format": format_type, "data": results}
        )
        return file_path

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """Exporta resultados de batch."""
        file_path = f"{destination}_batch.{format_type}"
        self._exported_files.append(
            {"path": file_path, "format": format_type, "data": batch_results}
        )
        return file_path

    def get_supported_formats(self) -> List[str]:
        """Lista formatos suportados."""
        return ["json", "csv", "xlsx", "parquet", "pickle"]

    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """Exporta resultados de otimização."""
        file_path = f"{destination}_optimization.json"
        self._exported_files.append(
            {"path": file_path, "format": "json", "data": optimization_data}
        )
        return file_path

    def get_exported_files(self) -> List[Dict[str, Any]]:
        """Retorna lista de arquivos exportados (para testes)."""
        return self._exported_files.copy()


class FakeExecutorPort:
    """Implementação fake do ExecutorPort para testes."""

    def __init__(self):
        self._execution_count = 0
        self._fail_on_algorithm = None  # Para simular falhas

    def set_fail_on_algorithm(self, algorithm_name: str):
        """Configura para falhar em algoritmo específico (para testes)."""
        self._fail_on_algorithm = algorithm_name

    def execute_batch(
        self,
        batch_config: Dict[str, Any],
        monitoring_service=None,
        resources_config=None,
    ) -> List[Dict[str, Any]]:
        """Executa batch de experimentos."""
        results = []

        # Support new batch structure
        if "execution" in batch_config:
            # New structure
            executions = batch_config["execution"].get("executions", [])
            for execution in executions:
                for dataset_id in execution.get("datasets", []):
                    for algo_config_id in execution.get("algorithms", []):
                        for repetition in range(execution.get("repetitions", 1)):
                            try:
                                # Simulate dataset loading
                                dataset = SyntheticDatasetGenerator.generate_random(
                                    n=10, length=20, alphabet="ACGT", seed=42
                                )

                                self._execution_count += 1

                                if algo_config_id == self._fail_on_algorithm:
                                    from src.domain import AlgorithmExecutionError

                                    raise AlgorithmExecutionError(
                                        f"Falha simulada para {algo_config_id}"
                                    )

                                result = {
                                    "algorithm_name": algo_config_id,
                                    "dataset_id": dataset_id,
                                    "dataset_size": dataset.size,
                                    "result": "ACGT" * (dataset.length // 4),
                                    "distance": 2,
                                    "execution_time": 1.5,
                                    "status": "success",
                                    "metadata": {"repetition": repetition + 1},
                                    "execution_name": execution["nome"],
                                }
                                results.append(result)

                            except Exception as e:
                                results.append(
                                    {
                                        "algorithm_name": algo_config_id,
                                        "dataset_id": dataset_id,
                                        "status": "error",
                                        "error": str(e),
                                        "execution_name": execution["nome"],
                                    }
                                )
        else:
            # Legacy structure
            for exp in batch_config.get("experiments", []):
                try:
                    # Simular carregamento de dataset
                    dataset = SyntheticDatasetGenerator.generate_random(
                        n=10, length=20, alphabet="ACGT", seed=42
                    )

                    # Simular execução direta em vez de chamar execute_single
                    self._execution_count += 1

                    if exp["algorithm"] == self._fail_on_algorithm:
                        from src.domain import AlgorithmExecutionError

                        raise AlgorithmExecutionError(
                            f"Falha simulada para {exp['algorithm']}"
                        )

                    result = {
                        "algorithm": exp["algorithm"],
                        "dataset_size": dataset.size,
                        "result": "ACGT" * (dataset.length // 4),
                        "distance": 2,
                        "execution_time": 1.5,
                        "status": "success",
                        "metadata": {"params": exp.get("params", {})},
                        "dataset": exp["dataset"],
                    }
                    results.append(result)

                except Exception as e:
                    results.append(
                        {
                            "algorithm": exp["algorithm"],
                            "dataset": exp["dataset"],
                            "status": "error",
                            "error": str(e),
                        }
                    )

        return results

    def execute_optimization(
        self,
        algorithm_name: str,
        dataset: Dataset,
        optimization_config: Dict[str, Any],
        monitoring_service=None,
        config_index: int = 1,
        total_configs: int = 1,
        dataset_index: int = 1,
        total_datasets: int = 1,
    ) -> Dict[str, Any]:
        """Executa otimização de hiperparâmetros."""
        return {
            "algorithm": algorithm_name,
            "dataset_size": dataset.size,
            "best_params": {"param1": 10, "param2": 0.5},
            "best_score": 3.2,
            "optimization_history": [
                {"params": {"param1": 5, "param2": 0.3}, "score": 4.1},
                {"params": {"param1": 8, "param2": 0.4}, "score": 3.8},
                {"params": {"param1": 10, "param2": 0.5}, "score": 3.2},
            ],
            "trials": 50,
            "status": "completed",
        }

    def execute_sensitivity_analysis(
        self, algorithm_name: str, dataset: Dataset, sensitivity_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Executa análise de sensibilidade."""
        return {
            "algorithm": algorithm_name,
            "dataset_size": dataset.size,
            "parameters_analyzed": list(
                sensitivity_config.get("parameters", {}).keys()
            ),
            "sensitivity_results": {
                "param1": {"sensitivity_index": 0.7, "impact": "high"},
                "param2": {"sensitivity_index": 0.3, "impact": "low"},
            },
            "status": "completed",
        }

    def get_execution_status(self, execution_id: str) -> str:
        """Obtém status de execução."""
        return "completed"

    def cancel_execution(self, execution_id: str) -> bool:
        """Cancela execução em andamento."""
        return True

    def get_execution_count(self) -> int:
        """Retorna número de execuções (para testes)."""
        return self._execution_count


class FakeMonitoringService:
    """Fake monitoring service for testing."""

    def __init__(self):
        self._is_monitoring = False
        self._task_type = None
        self._task_name = None
        self._interface = "simple"
        self._update_interval = 3
        self._errors = []
        self._progress_messages = []

    def start_monitoring(self, task_type, task_name: str):
        """Start monitoring a task."""
        self._is_monitoring = True
        self._task_type = task_type
        self._task_name = task_name

    def finish_monitoring(self, results: dict):
        """Finish monitoring."""
        self._is_monitoring = False

    def show_error(self, error_message: str):
        """Show an error message."""
        self._errors.append(error_message)

    def show_progress(self, message: str):
        """Show progress message."""
        self._progress_messages.append(message)

    def set_interface(self, interface: str):
        """Set monitoring interface."""
        self._interface = interface

    def set_update_interval(self, interval: int):
        """Set update interval."""
        self._update_interval = interval

    # Test helper methods
    def is_monitoring(self) -> bool:
        """Check if currently monitoring."""
        return self._is_monitoring

    def get_task_type(self):
        """Get current task type."""
        return self._task_type

    def get_task_name(self) -> str:
        """Get current task name."""
        return self._task_name

    def get_errors(self) -> List[str]:
        """Get list of errors."""
        return self._errors.copy()

    def get_progress_messages(self) -> List[str]:
        """Get list of progress messages."""
        return self._progress_messages.copy()

    def get_interface(self) -> str:
        """Get current interface."""
        return self._interface

    def get_update_interval(self) -> int:
        """Get current update interval."""
        return self._update_interval


class MockAlgorithmClass(CSPAlgorithm):
    """Classe mock para simular algoritmos nos testes."""

    name = "MockAlgorithm"
    default_params = {"param1": 10}
    is_deterministic = True
    supports_internal_parallel = False

    def __init__(self, strings, alphabet, **params):
        super().__init__(strings, alphabet, **params)

    def run(self):
        """Simula execução do algoritmo."""
        return "ACGT", 2, {"iterations": 100}
