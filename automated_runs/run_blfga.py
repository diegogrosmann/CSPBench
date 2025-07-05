#!/usr/bin/env python3
"""
Ponto de entrada para execução BLFGA automatizada.
"""

import os
import sys

# ensure project root is in path for module imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import argparse
import resource
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from functools import reduce
from itertools import product
from multiprocessing import cpu_count
from operator import mul
from typing import Any

import pandas as pd
import psutil
import yaml
from datasets.dataset_synthetic import generate_dataset_with_params

from algorithms.blf_ga.algorithm import BLFGAAlgorithm
from algorithms.blf_ga.config import BLF_GA_DEFAULTS

# Configurações
CONFIG_DIR = os.path.join(os.path.dirname(__file__), "configs")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


@dataclass
class ExperimentTask:
    """Representa uma tarefa de experimento para execução paralela."""

    dataset_idx: int
    exp_idx: int
    strings: list[str]
    alphabet: str
    params: dict[str, Any]
    dataset_info: dict[str, Any]


@dataclass
class ExperimentResult:
    """Resultado de um experimento."""

    dataset_idx: int
    exp_idx: int
    dataset_n: int
    dataset_L: int
    params: dict[str, Any]
    dist: float
    tempo: float
    memoria_usada: float
    success: bool
    error_msg: str = ""


def load_yaml(path: str) -> dict:
    """Carrega arquivo YAML."""
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def limit_worker_memory(max_mem_gb: float = 1.0):
    """Limita memória de um worker específico."""
    try:
        _soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        max_bytes = int(max_mem_gb * 1024**3)
        resource.setrlimit(resource.RLIMIT_AS, (max_bytes, hard))
    except Exception as e:
        print(f"[Warning] Não foi possível limitar memória do worker: {e}")


def execute_single_experiment(task: ExperimentTask) -> ExperimentResult:
    """
    Executa um único experimento BLF-GA.

    Args:
        task: Tarefa contendo todos os dados necessários

    Returns:
        ExperimentResult: Resultado do experimento
    """
    # Obter limite de memória dos argumentos ou usar padrão
    memory_limit = float(os.environ.get("BLFGA_MEMORY_LIMIT", "1.0"))

    # Limitar memória do worker
    limit_worker_memory(max_mem_gb=memory_limit)

    try:
        # Monitorar memória inicial
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / 1024**2

        # Executar algoritmo
        alg = BLFGAAlgorithm(task.strings, task.alphabet, **task.params)
        t0 = time.time()
        _, dist, _ = alg.run_with_history()
        t1 = time.time()

        # Monitorar memória final
        mem_after = process.memory_info().rss / 1024**2
        memoria_usada = mem_after - mem_before

        return ExperimentResult(
            dataset_idx=task.dataset_idx,
            exp_idx=task.exp_idx,
            dataset_n=len(task.strings),
            dataset_L=len(task.strings[0]),
            params=task.params.copy(),
            dist=dist,
            tempo=t1 - t0,
            memoria_usada=memoria_usada,
            success=True,
        )

    except Exception as e:  # pylint: disable=broad-except
        return ExperimentResult(
            dataset_idx=task.dataset_idx,
            exp_idx=task.exp_idx,
            dataset_n=len(task.strings) if task.strings else 0,
            dataset_L=len(task.strings[0]) if task.strings else 0,
            params=task.params.copy(),
            dist=float("inf"),
            tempo=0.0,
            memoria_usada=0.0,
            success=False,
            error_msg=str(e),
        )


def calculate_optimal_workers(
    total_experiments: int,
    available_memory_gb: float,
    user_workers: int | None = None,
) -> int:
    """
    Calcula o número ótimo de workers baseado em recursos disponíveis.

    Args:
        total_experiments: Número total de experimentos
        available_memory_gb: Memória disponível em GB
        user_workers: Número de workers especificado pelo usuário (opcional)

    Returns:
        int: Número ótimo de workers
    """
    # Se o usuário especificou um número, use-o (com validação básica)
    if user_workers is not None:
        if user_workers <= 0:
            print(f"⚠️ Número de workers inválido ({user_workers}). Usando padrão.")
        elif user_workers > cpu_count() * 2:
            print(f"⚠️ Muitos workers ({user_workers}). Limitando a {cpu_count() * 2}.")
            return min(user_workers, cpu_count() * 2)
        else:
            print(f"👤 Usando {user_workers} workers (especificado pelo usuário)")
            return user_workers

    # Limites baseados em recursos (padrão: usar todos os CPUs disponíveis)
    max_workers_cpu = cpu_count()  # Usar todos os CPUs por padrão
    max_workers_memory = max(1, int(available_memory_gb // 1.5))  # 1.5GB por worker
    max_workers_experiments = min(total_experiments, 16)  # Máximo 16 workers

    optimal = min(max_workers_cpu, max_workers_memory, max_workers_experiments)
    return max(1, optimal)


def parse_arguments():
    """Parse argumentos de linha de comando."""
    parser = argparse.ArgumentParser(
        description="Execução paralela do BLF-GA para múltiplos experimentos",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  python run_blfga_parallel.py                    # Usar todos os CPUs
  python run_blfga_parallel.py --workers 4        # Usar 4 workers
  python run_blfga_parallel.py --workers auto     # Detecção automática (padrão)

Variáveis de ambiente:
  BLFGA_WORKERS=4    # Definir número de workers via env var
        """,
    )

    parser.add_argument(
        "--workers",
        "-w",
        type=str,
        default=None,
        help='Número de workers paralelos. Use "auto" para detecção automática, ou um número específico (padrão: número de CPUs)',
    )

    parser.add_argument(
        "--memory-limit",
        type=float,
        default=1.0,
        help="Limite de memória por worker em GB (padrão: 1.0)",
    )

    parser.add_argument(
        "--max-experiments",
        type=int,
        default=1000,
        help="Limite máximo de experimentos (padrão: 1000)",
    )

    return parser.parse_args()


def get_worker_count(args) -> int | None:
    """Determina o número de workers baseado em argumentos e variáveis de ambiente."""

    # 1. Verificar argumento de linha de comando
    if args.workers is not None:
        if args.workers.lower() == "auto":
            return None  # Usar detecção automática
        try:
            return int(args.workers)
        except ValueError:
            print(f"⚠️ Valor inválido para --workers: '{args.workers}'. Usando detecção automática.")
            return None

    # 2. Verificar variável de ambiente
    env_workers = os.environ.get("BLFGA_WORKERS")
    if env_workers:
        if env_workers.lower() == "auto":
            return None
        try:
            workers = int(env_workers)
            print(f"🌍 Usando {workers} workers (variável de ambiente BLFGA_WORKERS)")
            return workers
        except ValueError:
            print(f"⚠️ Valor inválido em BLFGA_WORKERS: '{env_workers}'. Usando detecção automática.")

    # 3. Padrão: usar detecção automática (todos os CPUs)
    return None


def main():
    """Função principal com execução paralela."""
    # Parse argumentos
    args = parse_arguments()
    user_workers = get_worker_count(args)

    print(f"[Sistema] PID: {os.getpid()}")
    print(f"[Sistema] CPUs disponíveis: {cpu_count()}")

    # Monitorar recursos do sistema
    memory = psutil.virtual_memory()
    available_memory_gb = memory.available / (1024**3)
    print(f"[Sistema] Memória disponível: {available_memory_gb:.1f} GB")

    # Carregar configurações
    dataset_params_config = load_yaml(os.path.join(CONFIG_DIR, "dataset.yaml"))
    blfga_param_grid = load_yaml(os.path.join(CONFIG_DIR, "blfga_grid.yaml"))
    if blfga_param_grid is None:
        blfga_param_grid = {}

    # Gerar combinações de datasets
    dataset_param_names = list(dataset_params_config.keys())
    dataset_param_values = [val if isinstance(val, list) else [val] for val in dataset_params_config.values()]
    dataset_configs = [dict(zip(dataset_param_names, values)) for values in product(*dataset_param_values)]

    print(f"[Config] Datasets: {len(dataset_configs)}")
    print(
        f"[Config] Parâmetros: n={dataset_params_config.get('n')}, "
        f"L={dataset_params_config.get('L')}, "
        f"alphabet='{dataset_params_config.get('alphabet')}'"
    )

    # Preparar tarefas de experimentos
    all_tasks = []
    total_experiments = 0

    for ds_idx, ds_params in enumerate(dataset_configs):
        # Gerar dataset
        strings, params_usados = generate_dataset_with_params(ds_params)
        print(f"[Dataset {ds_idx+1}] n={len(strings)}, L={len(strings[0])}, " f"|Σ|={len(params_usados['alphabet'])}")

        # Preparar grid BLF-GA
        param_names = list(BLF_GA_DEFAULTS.keys())
        param_values = [
            (
                blfga_param_grid.get(k, [BLF_GA_DEFAULTS[k]])
                if isinstance(blfga_param_grid.get(k, BLF_GA_DEFAULTS[k]), list)
                else [blfga_param_grid.get(k, BLF_GA_DEFAULTS[k])]
            )
            for k in param_names
        ]

        # Calcular experimentos para este dataset
        dataset_experiments = reduce(mul, [len(v) for v in param_values], 1)
        total_experiments += dataset_experiments

        print(f"[Dataset {ds_idx+1}] Experimentos: {dataset_experiments}")

        # Criar tarefas
        for exp_idx, valores in enumerate(product(*param_values)):
            params = dict(zip(param_names, valores))
            task = ExperimentTask(
                dataset_idx=ds_idx,
                exp_idx=exp_idx,
                strings=strings,
                alphabet=params_usados["alphabet"],
                params=params,
                dataset_info=params_usados,
            )
            all_tasks.append(task)

    print(f"[Total] Experimentos: {total_experiments}")

    # Verificar limites de segurança
    max_experiments_limit = args.max_experiments
    if total_experiments > max_experiments_limit:
        print(
            f"❌ Número de experimentos ({total_experiments}) excede o limite seguro ({max_experiments_limit}). "
            "Reduza o grid de parâmetros ou use --max-experiments."
        )
        return

    # Calcular workers ótimos
    optimal_workers = calculate_optimal_workers(total_experiments, available_memory_gb, user_workers)
    print(f"[Paralelização] Workers utilizados: {optimal_workers}")

    # Informações adicionais sobre a configuração
    efficiency = total_experiments / optimal_workers if optimal_workers > 0 else 0
    print(f"[Paralelização] Experimentos por worker: {efficiency:.1f}")
    print(f"[Paralelização] Memória estimada total: {optimal_workers * args.memory_limit:.1f} GB")

    # Executar experimentos em paralelo
    all_results = []
    start_time = time.time()

    # Configurar variável de ambiente para workers
    os.environ["BLFGA_MEMORY_LIMIT"] = str(args.memory_limit)

    try:
        # Importar tqdm se disponível para progress bar
        try:
            from tqdm import tqdm

            progress_bar = tqdm(total=total_experiments, desc="Experimentos")
        except ImportError:
            progress_bar = None
            print("[Info] tqdm não disponível. Instale com: pip install tqdm")

        with ProcessPoolExecutor(max_workers=optimal_workers) as executor:
            # Submeter todas as tarefas
            future_to_task = {executor.submit(execute_single_experiment, task): task for task in all_tasks}

            # Coletar resultados conforme completam
            for future in as_completed(future_to_task):
                result = future.result()
                all_results.append(result)

                if progress_bar:
                    progress_bar.update(1)
                else:
                    print(f"[Progresso] {len(all_results)}/{total_experiments} concluídos")

                if not result.success:
                    print(f"❌ Erro no experimento DS{result.dataset_idx+1}-{result.exp_idx+1}: " f"{result.error_msg}")

        if progress_bar:
            progress_bar.close()

    except KeyboardInterrupt:
        print("\n❌ Execução interrompida pelo usuário!")
        return
    except Exception as e:  # pylint: disable=broad-except
        print(f"❌ Erro na execução paralela: {e}")
        return

    total_time = time.time() - start_time
    print(f"\n✅ Execução concluída em {total_time:.1f}s")

    # Processar e salvar resultados
    if all_results:
        successful_results = [r for r in all_results if r.success]
        failed_results = [r for r in all_results if not r.success]

        print(f"[Resultados] Sucessos: {len(successful_results)}, " f"Falhas: {len(failed_results)}")

        if successful_results:
            # Converter para DataFrame
            result_data = []
            for result in successful_results:
                row = {
                    "dataset_idx": result.dataset_idx,
                    "exp_idx": result.exp_idx,
                    "dataset_n": result.dataset_n,
                    "dataset_L": result.dataset_L,
                    **result.params,
                    "dist": result.dist,
                    "tempo": result.tempo,
                    "memoria_mb": result.memoria_usada,
                }
                result_data.append(row)

            df = pd.DataFrame(result_data)

            # Salvar resultados
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            csv_path = os.path.join(RESULTS_DIR, f"blfga_parallel_{timestamp}.csv")
            df.to_csv(csv_path, index=False)

            # Estatísticas
            print("\n📊 Estatísticas:")
            print(f"   Tempo médio por experimento: {df['tempo'].mean():.2f}s")
            print(f"   Melhor distância encontrada: {df['dist'].min()}")
            print(f"   Memória média por experimento: {df['memoria_mb'].mean():.1f} MB")
            print(f"   Speedup estimado: {len(successful_results) * df['tempo'].mean() / total_time:.1f}x")

            print(f"\n💾 Resultados salvos em: {csv_path}")

        if failed_results:
            print(f"\n⚠️ {len(failed_results)} experimentos falharam")
            error_summary = {}
            for result in failed_results:
                error_type = type(Exception(result.error_msg)).__name__
                error_summary[error_type] = error_summary.get(error_type, 0) + 1
            for error_type, count in error_summary.items():
                print(f"   {error_type}: {count} ocorrências")
    else:
        print("❌ Nenhum resultado foi gerado.")


if __name__ == "__main__":
    main()
