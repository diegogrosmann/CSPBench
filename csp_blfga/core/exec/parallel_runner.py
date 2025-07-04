"""
Executor paralelo para múltiplos algoritmos CSP.

Classes:
    ParallelRunner: Gerencia execução paralela de algoritmos com tqdm.

Funções:
    execute_algorithms_parallel(...): Executa múltiplos algoritmos em paralelo.
"""

import logging
import time
from typing import Any

from algorithms.base import global_registry
from csp_blfga.core.exec.algorithm_executor import ModernParallelExecutor
from csp_blfga.core.exec.runner import ProgressTracker
from csp_blfga.utils.resource_monitor import force_garbage_collection

logger = logging.getLogger(__name__)


class ParallelRunner:
    """
    Gerenciador de execução paralela de algoritmos CSP.

    Permite executar múltiplos algoritmos simultaneamente com
    monitoramento de progresso e controle de recursos.
    """

    def __init__(self, max_workers: int | None = None, timeout: int = 300):
        """
        Inicializa o runner paralelo.

        Args:
            max_workers: Número máximo de workers paralelos
            timeout: Timeout por algoritmo em segundos
        """
        self.executor = ModernParallelExecutor(max_workers, timeout)
        self.timeout = timeout

    def execute_algorithms_parallel(
        self,
        algorithm_names: list[str],
        seqs: list[str],
        alphabet: str,
        console=None,
        baseline_val: float | None = None,
    ) -> dict[str, dict[str, Any]]:
        """
        Executa múltiplos algoritmos em paralelo.

        Args:
            algorithm_names: Lista de nomes dos algoritmos
            seqs: Sequências de entrada
            alphabet: Alfabeto utilizado
            console: Console manager para output
            baseline_val: Valor baseline para comparação

        Returns:
            Dicionário com resultados de cada algoritmo
        """
        logger.info(f"Iniciando execução paralela de {len(algorithm_names)} algoritmos")

        # Força garbage collection antes de começar
        force_garbage_collection()

        # Preparar tarefas
        tasks = []
        valid_algorithms = []

        for alg_name in algorithm_names:
            if alg_name not in global_registry:
                if console:
                    console.print(
                        f"⚠️ Algoritmo '{alg_name}' não encontrado no registry"
                    )
                continue

            alg_class = global_registry[alg_name]
            tasks.append(
                {
                    "alg_class": alg_class,
                    "strings": seqs,
                    "alphabet": alphabet,
                    "params": {},
                    "name": alg_name,
                }
            )
            valid_algorithms.append(alg_name)

        if not tasks:
            if console:
                console.print("❌ Nenhum algoritmo válido encontrado")
            return {}

        # Criar tracker de progresso
        progress_tracker = ProgressTracker(
            f"Executando {len(tasks)} algoritmos", total=len(tasks), console=console
        )
        progress_tracker.start()

        results = {}

        # Executar em paralelo
        start_time = time.time()

        try:
            parallel_results = self.executor.execute_algorithm_batch(tasks)
            execution_time = time.time() - start_time

            # Processar resultados
            for i, result in enumerate(parallel_results):
                alg_name = valid_algorithms[i]

                if result.get("success", False):
                    # Sucesso
                    dist_status = ""
                    if (
                        baseline_val is not None
                        and result.get("distance", float("inf")) <= baseline_val
                    ):
                        dist_status = " ✓"

                    progress_tracker.update(
                        1,
                        f"{alg_name}: dist={result.get('distance', 'inf')}{dist_status}",
                    )

                    results[alg_name] = {
                        "success": True,
                        "center": result.get("center", ""),
                        "distance": result.get("distance", float("inf")),
                        "metadata": result.get("metadata", {}),
                        "tempo": result.get("execution_time", 0.0),
                        "iteracoes": result.get("metadata", {}).get("iteracoes", 0),
                        "erro": None,
                    }
                else:
                    # Falha
                    error_msg = result.get("error", "Erro desconhecido")
                    progress_tracker.update(1, f"{alg_name}: {error_msg}")

                    results[alg_name] = {
                        "success": False,
                        "center": None,
                        "distance": float("inf"),
                        "metadata": {},
                        "tempo": result.get("execution_time", 0.0),
                        "iteracoes": 0,
                        "erro": error_msg,
                    }

        except Exception as e:
            logger.error(f"Erro na execução paralela: {e}")
            progress_tracker.finish(f"Erro: {e}")
            return {}

        # Finalizar progresso
        successful_algs = len([r for r in results.values() if r["success"]])
        total_algs = len(results)

        final_msg = f"{successful_algs}/{total_algs} algoritmos executados com sucesso"
        progress_tracker.finish(final_msg)

        return results


def execute_algorithms_parallel(
    algorithm_names: list[str],
    seqs: list[str],
    alphabet: str,
    console=None,
    baseline_val: float | None = None,
    max_workers: int | None = None,
    timeout: int = 300,
) -> dict[str, dict[str, Any]]:
    """
    Função de conveniência para execução paralela de algoritmos.

    Args:
        algorithm_names: Lista de nomes dos algoritmos
        seqs: Sequências de entrada
        alphabet: Alfabeto utilizado
        console: Console manager para output
        baseline_val: Valor baseline para comparação
        max_workers: Número máximo de workers paralelos
        timeout: Timeout por algoritmo em segundos

    Returns:
        Dicionário com resultados de cada algoritmo
    """
    runner = ParallelRunner(max_workers=max_workers, timeout=timeout)
    return runner.execute_algorithms_parallel(
        algorithm_names=algorithm_names,
        seqs=seqs,
        alphabet=alphabet,
        console=console,
        baseline_val=baseline_val,
    )
