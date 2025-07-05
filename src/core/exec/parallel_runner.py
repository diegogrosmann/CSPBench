"""
Executor paralelo para múltiplos algoritmos CSP.

Este módulo implementa execução paralela de algoritmos CSP com controle
de recursos, monitoramento de progresso e encerramento gracioso.

Classes:
    ParallelRunner: Gerencia execução paralela de algoritmos com tqdm.

Funções:
    execute_algorithms_parallel(...): Executa múltiplos algoritmos em paralelo.

Funcionalidades:
    - Execução paralela com ProcessPoolExecutor
    - Monitoramento de progresso com tqdm
    - Controle de recursos e timeouts
    - Encerramento gracioso via sinais
    - Análise de resultados em tempo real
"""

import logging
import time
from typing import Any

from algorithms.base import global_registry
from src.core.exec.algorithm_executor import ParallelAlgorithmExecutor
from src.core.exec.runner import ProgressTracker
from src.utils.resource_monitor import force_garbage_collection
from src.utils.signal_manager import is_interrupted, register_shutdown_callback

logger = logging.getLogger(__name__)


class ParallelRunner:
    """
    Gerenciador de execução paralela de algoritmos CSP.

    Esta classe implementa execução paralela robusta com:
    - Controle de recursos e monitoramento
    - Integração com sistema de sinais para encerramento gracioso
    - Progresso em tempo real com tqdm
    - Análise de resultados automática

    Attributes:
        max_workers: Número máximo de workers paralelos
        timeout: Timeout por algoritmo em segundos
        _shutdown_requested: Flag de shutdown solicitado
        _current_executor: Referência ao executor atual

    Examples:
        >>> runner = ParallelRunner(max_workers=4, timeout=300)
        >>> results = runner.execute_algorithms_parallel(
        ...     algorithm_names=['BLF-GA', 'Baseline'],
        ...     seqs=sequences,
        ...     alphabet='ACGT'
        ... )
    """

    def __init__(self, max_workers: int | None = None, timeout: int = 300):
        """
        Inicializa o runner paralelo.

        Args:
            max_workers: Número máximo de workers paralelos
            timeout: Timeout por algoritmo em segundos
        """
        self.max_workers = max_workers
        self.timeout = timeout
        self._shutdown_requested = False
        self._current_executor = None

        # Registrar callback de shutdown para encerrar pool graciosamente
        register_shutdown_callback(self._shutdown_callback)

    def _shutdown_callback(self):
        """Callback para encerramento gracioso do pool de processos."""
        logger.info("ParallelRunner: Iniciando encerramento gracioso")
        self._shutdown_requested = True

        if self._current_executor:
            logger.info("ParallelRunner: Encerrando pool de processos...")
            try:
                if hasattr(self._current_executor, "executor") and self._current_executor.executor:
                    self._current_executor.executor.shutdown(wait=False)
                logger.info("ParallelRunner: Pool de processos encerrado")
            except Exception as e:
                logger.error(f"Erro ao encerrar pool: {e}")

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
                    console.print(f"⚠️ Algoritmo '{alg_name}' não encontrado no registry")
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
        progress_tracker = ProgressTracker(f"Executando {len(tasks)} algoritmos", total=len(tasks), console=console)
        progress_tracker.start()

        results = {}

        # Executar em paralelo
        start_time = time.time()

        try:
            # Verificar se shutdown foi solicitado
            if self._shutdown_requested or is_interrupted():
                logger.warning("Shutdown solicitado - cancelando execução")
                if console:
                    console.print("⚠️ Execução cancelada pelo usuário")
                return {}

            # Usar context manager para garantir limpeza
            with ParallelAlgorithmExecutor(max_workers=self.max_workers, timeout_seconds=self.timeout) as executor:
                # Guardar referência para callback de shutdown
                self._current_executor = executor

                parallel_results = executor.execute_parallel(tasks)
                execution_time = time.time() - start_time

            # Processar resultados
            for i, result in enumerate(parallel_results):
                # Verificar se houve interrupção durante o processamento
                if self._shutdown_requested or is_interrupted():
                    logger.warning("Interrupção detectada durante processamento dos resultados")
                    break

                alg_name = valid_algorithms[i]

                if result.get("success", False):
                    # Sucesso
                    dist_status = ""
                    if baseline_val is not None and result.get("distance", float("inf")) <= baseline_val:
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
        finally:
            # Limpar referência do executor
            self._current_executor = None

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
