"""
Utilitários para execução de algoritmos e exibição de progresso.

Classes:
    ProgressTracker: Gerencia progresso com tqdm e mensagens através do ConsoleManager.

Funções:
    execute_algorithm_runs(...): Executa múltiplas execuções de um algoritmo, coletando resultados.
"""

import logging
import time

from csp_blfga.core.exec.algorithm_executor import AlgorithmExecutor
from csp_blfga.utils.config import ALGORITHM_TIMEOUT
from csp_blfga.utils.resource_monitor import (
    check_algorithm_feasibility,
    force_garbage_collection,
)

try:
    from tqdm import tqdm

    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    print("⚠️ tqdm não disponível. Instale com: pip install tqdm")


class ProgressTracker:
    """
    Gerenciador de progresso moderno usando tqdm para barras de progresso
    e ConsoleManager para mensagens fora da barra.
    """

    def __init__(self, description: str, total: int | None = None, console=None):
        self.description = description
        self.total = total
        self.console = console
        self.pbar = None
        self.current = 0
        self.started = False

    def start(self):
        """Inicia a barra de progresso."""
        if self.started:
            return

        if TQDM_AVAILABLE and self.total is not None:
            self.pbar = tqdm(
                total=self.total,
                desc=self.description,
                unit="exec",
                leave=False,
                ncols=80,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
            )
        elif self.console:
            self.console.print(f"🔄 {self.description}...")
        else:
            print(f"🔄 {self.description}...")

        self.started = True

    def update(self, n: int = 1, message: str | None = None):
        """
        Atualiza o progresso.

        Args:
            n: Incremento no progresso
            message: Mensagem adicional para exibir
        """
        self.current += n

        if self.pbar:
            self.pbar.update(n)
            if message:
                self.pbar.set_postfix_str(message)
        elif message:
            # Usar console para mensagens quando tqdm não está disponível
            if self.console:
                self.console.print(f"  {message}")
            else:
                print(f"  {message}")

    def set_description(self, description: str):
        """Atualiza a descrição da barra de progresso."""
        self.description = description
        if self.pbar:
            self.pbar.set_description(description)

    def finish(self, final_message: str | None = None):
        """Finaliza a barra de progresso."""
        if self.pbar:
            if final_message:
                self.pbar.set_postfix_str(final_message)
            self.pbar.close()
            self.pbar = None
        elif final_message:
            if self.console:
                self.console.print(f"✅ {final_message}")
            else:
                print(f"✅ {final_message}")

        self.started = False


def execute_algorithm_runs(
    alg_name,
    AlgClass,
    seqs,
    alphabet,
    num_execs,
    baseline_val,
    console=None,
    timeout=None,
):
    """
    Executa múltiplas execuções de um algoritmo com monitoramento de progresso moderno.

    Args:
        alg_name: Nome do algoritmo
        AlgClass: Classe do algoritmo
        seqs: Sequências de entrada
        alphabet: Alfabeto utilizado
        num_execs: Número de execuções
        baseline_val: Valor baseline para comparação
        console: Instância do ConsoleManager para mensagens
        timeout: Timeout por execução

    Returns:
        Lista de resultados das execuções
    """
    logger = logging.getLogger(__name__)
    logger.debug(f"[Runner] Iniciando {alg_name} (baseline_val={baseline_val})")

    if timeout is None:
        timeout = ALGORITHM_TIMEOUT

    # Força garbage collection antes de começar
    force_garbage_collection()

    # Verificar viabilidade do algoritmo antes da execução
    n, L = len(seqs), len(seqs[0])
    is_feasible, feasibility_msg = check_algorithm_feasibility(n, L, alg_name)

    if not is_feasible:
        if console:
            console.print(f"❌ {alg_name}: INVIÁVEL - {feasibility_msg}")
        else:
            print(f"❌ {alg_name}: INVIÁVEL - {feasibility_msg}")

        return [
            {
                "tempo": 0.0,
                "iteracoes": 0,
                "distancia": float("inf"),
                "melhor_string": "",
                "erro": f"Inviável: {feasibility_msg}",
            }
        ]

    is_deterministic = getattr(AlgClass, "is_deterministic", False)
    actual_execs = 1 if is_deterministic else num_execs
    executions = []

    # Criar tracker de progresso para todas as execuções
    progress_tracker = ProgressTracker(
        f"{alg_name}", total=actual_execs, console=console
    )
    progress_tracker.start()

    for i in range(actual_execs):
        exec_description = f"{alg_name}"
        if actual_execs > 1:
            exec_description += f" ({i+1}/{actual_execs})"

        executor = AlgorithmExecutor(timeout)
        warning_holder = []

        def warning_callback(msg):
            warning_holder.append(msg)
            # T8-2: Usar console para warnings fora da barra de progresso
            if console:
                console.print(f"⚠️ {msg}")

        def progress_callback(msg: str):
            """Callback de progresso do algoritmo."""
            # T8-2: Usar console para mensagens de progresso detalhadas
            if console:
                console.print(f"  📊 {msg}")

        t0 = time.time()

        try:
            # Criar instância do algoritmo
            instance = AlgClass(seqs, alphabet)

            # Executar com timeout
            center, val, info = executor.execute_with_timeout(
                instance,
                progress_callback=progress_callback,
                warning_callback=warning_callback,
            )

            tempo_execucao = time.time() - t0

            if "erro" in info:
                if info["erro"] == "timeout" or "timeout" in info["erro"].lower():
                    # Atualizar progresso com timeout
                    progress_tracker.update(1, f"TIMEOUT ({timeout}s)")
                    executions.append(
                        {
                            "tempo": tempo_execucao,
                            "iteracoes": 0,
                            "distancia": float("inf"),
                            "melhor_string": "",
                            "erro": f"Timeout ({timeout}s)",
                        }
                    )
                    continue
                else:
                    # Erro genérico
                    progress_tracker.update(1, f"ERRO: {info['erro']}")
                    executions.append(
                        {
                            "tempo": tempo_execucao,
                            "iteracoes": info.get("iteracoes", 0),
                            "distancia": float("inf"),
                            "melhor_string": "",
                            "erro": info["erro"],
                        }
                    )
                    continue

            # Execução bem-sucedida
            dist_status = ""
            if baseline_val is not None and val <= baseline_val:
                dist_status = " ✓"

            progress_tracker.update(1, f"dist={val}{dist_status}")

            executions.append(
                {
                    "tempo": tempo_execucao,
                    "iteracoes": info.get("iteracoes", 0),
                    "distancia": val,
                    "melhor_string": center,
                    "erro": "",
                }
            )

            # Limpeza de memória a cada execução
            force_garbage_collection()

        except Exception as e:
            tempo_execucao = time.time() - t0
            progress_tracker.update(1, f"ERRO: {str(e)}")

            executions.append(
                {
                    "tempo": tempo_execucao,
                    "iteracoes": 0,
                    "distancia": float("inf"),
                    "melhor_string": "",
                    "erro": str(e),
                }
            )

    # Finalizar progresso
    successful_execs = len([e for e in executions if not e["erro"]])
    avg_dist = sum(
        e["distancia"] for e in executions if e["distancia"] != float("inf")
    ) / max(successful_execs, 1)

    final_msg = f"{successful_execs}/{actual_execs} sucessos"
    if successful_execs > 0:
        final_msg += f", avg_dist={avg_dist:.1f}"

    progress_tracker.finish(final_msg)

    return executions
