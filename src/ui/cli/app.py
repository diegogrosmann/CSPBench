"""
Arquivo principal de execução da aplicação Closest String Problem (CSP).

Este módulo orquestra o fluxo principal da aplicação, incluindo:
    1. Seleção e leitura/geração do dataset.
    2. Seleção dos algoritmos a executar.
    3. Execução dos algoritmos selecionados.
    4. Exibição e salvamento dos resultados.

A aplicação pode ser executada em modo interativo ou automatizado (para testes).

Attributes:
    Nenhum atributo global relevante.
"""

import argparse
import logging
import os
import signal
import sys
import traceback
import uuid
from datetime import datetime

from algorithms.base import global_registry
from src.core.io.results_formatter import ResultsFormatter
from src.core.report.report_utils import print_quick_summary
from src.ui.cli.console_manager import console
from src.ui.cli.menu import (
    configure_optimization_params,
    configure_sensitivity_params,
    menu,
    select_algorithms,
    select_optimization_algorithm,
    select_sensitivity_algorithm,
)
from src.ui.cli.save_wizard import ask_save_dataset
from src.utils.config import ALGORITHM_TIMEOUT, safe_input
from src.utils.curses_console import create_console_manager
from src.utils.logging import setup_logging
from src.utils.resource_monitor import (
    check_algorithm_feasibility,
    get_safe_memory_limit,
)

RESULTS_DIR = "results"
LOGS_DIR = "logs"
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)


def signal_handler(signum, frame):
    """Handler para sinais de interrupção.

    Args:
        signum (int): Número do sinal recebido.
        frame (frame object): Frame atual de execução.
    """
    print("\n\nOperação cancelada pelo usuário. Encerrando.")
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser(description="Closest String Problem (CSP) - Execução principal")
    parser.add_argument(
        "--silent",
        action="store_true",
        help="Executa em modo silencioso (sem prints interativos)",
    )
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["synthetic", "file", "entrez", "batch"],
        help="Fonte do dataset",
    )
    parser.add_argument(
        "--algorithms",
        type=str,
        nargs="+",
        help="Algoritmos a executar (nomes separados por espaço)",
    )
    parser.add_argument("--num-execs", type=int, help="Número de execuções por algoritmo")
    parser.add_argument("--timeout", type=int, help="Timeout por execução (segundos)")
    parser.add_argument("--workers", "-w", type=int, default=4, help="Número de workers paralelos (padrão: 4)")
    parser.add_argument(
        "--console",
        type=str,
        choices=["auto", "curses", "simple"],
        default="auto",
        help="Tipo de console: auto (detecta automaticamente), curses (forçar curses), simple (console simples)",
    )
    args = parser.parse_args()

    silent = args.silent

    # Criar console manager apropriado
    if silent:
        # Em modo silencioso, sempre usar console simples
        console_manager = create_console_manager(force_simple=True)
    else:
        # Detectar ou forçar tipo de console
        if args.console == "simple":
            console_manager = create_console_manager(force_simple=True)
        elif args.console == "curses":
            console_manager = create_console_manager(force_simple=False)
        else:  # auto
            console_manager = create_console_manager(force_simple=False)

    def silent_print(*a, **k):
        pass

    p = print if not silent else silent_print
    cprint = console_manager.print if not silent else silent_print

    # Mostrar o número do processo (PID) logo no início
    p(f"[PID] Processo em execução: {os.getpid()}")

    # Configurar handlers de sinal para saída limpa
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        uid = uuid.uuid4().hex[:8]
        base_name = f"{ts}_{uid}"
        # Passa silent para setup_logging
        setup_logging(base_name, silent=silent)

        # Se modo silencioso, use defaults para qualquer parâmetro não informado
        if silent:
            if not args.dataset:
                args.dataset = "synthetic"
            if not args.algorithms:
                args.algorithms = ["Baseline"]  # Usar Baseline em vez de BLF-GA para testes
            if not args.num_execs:
                args.num_execs = 1
            if not args.timeout:
                args.timeout = ALGORITHM_TIMEOUT

        # Dataset
        params = {}
        seqs = []
        seed = None

        # Escolha do dataset
        if args.dataset:
            choice = args.dataset
            if args.dataset == "synthetic":
                from src.datasets.dataset_synthetic import generate_dataset

                seqs, p = generate_dataset(silent=silent)
                logging.debug(f"[main] Parâmetros do dataset: {p}")
                params = {"dataset_source": "1"}
                params.update(p)
            elif args.dataset == "file":
                from src.datasets.dataset_file import load_dataset

                seqs, p = load_dataset(silent=silent)
                params = {"dataset_source": "2"}
                params.update(p)
            elif args.dataset == "entrez":
                from src.datasets.dataset_entrez import fetch_dataset

                seqs, p = fetch_dataset()
                params = {"dataset_source": "3"}
                params.update(p)
            elif args.dataset == "batch":
                from src.core.exec.batch_executor import (
                    BatchExecutor,
                    select_batch_config,
                )

                # Em modo silencioso, usar o arquivo de teste
                if silent:
                    config_file = "batch_configs/teste.yaml"
                else:
                    config_file = select_batch_config()
                    if not config_file:
                        cprint("❌ Nenhum arquivo de configuração selecionado.")
                        return
                executor = BatchExecutor(config_file, workers=args.workers)
                batch_result = executor.execute_batch()
                cprint("\n✅ Execução em lote concluída!")
                cprint(f"Tempo total: {batch_result['tempo_total']:.1f}s")
                cprint(f"Taxa de sucesso: {batch_result['resumo']['taxa_sucesso']:.1f}%")
                from src.core.io.exporter import CSPExporter

                batch_dir = executor.results_dir
                import os as os_batch

                json_path = os_batch.path.join(batch_dir, "batch_results.json")
                csv_path = os_batch.path.join(batch_dir, "batch_results.csv")
                if os_batch.path.exists(json_path):
                    exporter = CSPExporter()
                    exporter.export_batch_json_to_csv(json_path, csv_path)
                    cprint(f"📄 Resultados detalhados do batch exportados para CSV: {csv_path}")
                return
            else:
                cprint("❌ Fonte de dataset inválida.")
                return
        else:
            # Fluxo interativo
            choice = menu()
            params = {"dataset_source": choice}

            if choice == "4":
                from src.core.exec.batch_executor import (
                    BatchExecutor,
                    select_batch_config,
                )

                config_file = select_batch_config()
                if not config_file:
                    cprint("❌ Nenhum arquivo de configuração selecionado.")
                    return
                try:
                    executor = BatchExecutor(config_file, workers=args.workers)
                    batch_result = executor.execute_batch()
                    cprint("\n✅ Execução em lote concluída!")
                    cprint(f"Tempo total: {batch_result['tempo_total']:.1f}s")
                    cprint(f"Taxa de sucesso: {batch_result['resumo']['taxa_sucesso']:.1f}%")
                    from src.core.io.exporter import CSPExporter

                    batch_dir = executor.results_dir
                    import os as os_batch

                    json_path = os_batch.path.join(batch_dir, "batch_results.json")
                    csv_path = os_batch.path.join(batch_dir, "batch_results.csv")
                    if os_batch.path.exists(json_path):
                        exporter = CSPExporter()
                        exporter.export_batch_json_to_csv(json_path, csv_path)
                        cprint(f"📄 Resultados detalhados do batch exportados para CSV: {csv_path}")
                except Exception as e:
                    cprint(f"❌ Erro na execução em lote: {e}")
                    return
            elif choice == "5":
                # Execução em lote com interface curses
                from src.ui.curses_integration import add_curses_batch_option_to_menu

                curses_executor = add_curses_batch_option_to_menu()
                curses_executor()
                return
            elif choice == "6":
                # Otimização de hiperparâmetros
                try:
                    run_optimization_workflow()
                except Exception as e:
                    cprint(f"❌ Erro na otimização: {e}")
                    logging.exception("Erro na otimização", exc_info=e)
                return
            elif choice == "7":
                # Análise de sensibilidade
                try:
                    run_sensitivity_workflow()
                except Exception as e:
                    cprint(f"❌ Erro na análise de sensibilidade: {e}")
                    logging.exception("Erro na análise de sensibilidade", exc_info=e)
                return
            elif choice in ["1", "2", "3"]:
                try:
                    if choice == "1":
                        from src.datasets.dataset_synthetic import generate_dataset

                        seqs, p = generate_dataset(silent=silent)
                        logging.debug(f"[main] Parâmetros do dataset: {p}")
                        params.update(p)
                        if not silent:
                            ask_save_dataset(seqs, "synthetic", p)
                    elif choice == "2":
                        from src.datasets.dataset_file import load_dataset

                        seqs, p = load_dataset(silent=silent)
                        params.update(p)
                    elif choice == "3":
                        from src.datasets.dataset_entrez import fetch_dataset

                        seqs, p = fetch_dataset()
                        params.update(p)
                        ask_save_dataset(seqs, "entrez", p)
                except Exception as exc:
                    cprint(f"Erro: {exc}")
                    logging.exception("Erro ao carregar dataset", exc_info=exc)
                    return

        # Após carregar o dataset, extrair informações extras se disponíveis
        seed = params.get("seed")
        logging.debug(f"[main] seed: {seed}")

        # Exibir distância da string base se disponível
        distancia_base = params.get("distancia_string_base")
        if distancia_base is not None:
            cprint(f"Distância da string base: {distancia_base}")

        if not seqs:
            cprint("Nenhuma sequência lida.")
            return

        alphabet = "".join(sorted(set("".join(seqs))))
        n, L = len(seqs), len(seqs[0])
        cprint(f"\nDataset: n={n}, L={L}, |Σ|={len(alphabet)}")

        # Log simplificado do dataset
        logging.debug(f"[DATASET] n={n}, L={L}, |Σ|={len(alphabet)}")
        if len(seqs) <= 5:
            logging.debug(f"[DATASET] Strings: {seqs}")
        else:
            logging.debug(f"[DATASET] {len(seqs)} strings (primeiras 2: {seqs[:2]})")

        # Verificação de recursos do sistema
        safe_memory = get_safe_memory_limit()
        cprint(f"Limite seguro de memória: {safe_memory:.1f}%")

        # Algoritmos
        if args.algorithms:
            algs = args.algorithms
        else:
            algs = select_algorithms()
        if not algs:
            cprint("Nenhum algoritmo selecionado.")
            return

        # Verificar viabilidade dos algoritmos selecionados
        viable_algs = []
        for alg_name in algs:
            is_viable, msg = check_algorithm_feasibility(n, L, alg_name)
            if is_viable:
                viable_algs.append(alg_name)
                cprint(f"✓ {alg_name}: {msg}")
            else:
                cprint(f"⚠ {alg_name}: {msg} (será pulado)")

        if not viable_algs:
            cprint("Nenhum algoritmo viável.")
            return

        # Número de execuções e timeout
        if args.num_execs:
            num_execs = args.num_execs
        else:
            runs = safe_input("\nNº execuções p/ algoritmo [3]: ")
            num_execs = int(runs) if runs.isdigit() and int(runs) > 0 else 3
        if args.timeout:
            timeout = args.timeout
        else:
            default_timeout = ALGORITHM_TIMEOUT
            timeout_input = safe_input(f"\nTimeout por execução em segundos [{default_timeout}]: ")
            timeout = int(timeout_input) if timeout_input.isdigit() and int(timeout_input) > 0 else default_timeout
        cprint(f"Timeout configurado: {timeout}s por execução")

        # Execução dos algoritmos
        cprint("\n" + "=" * 50)
        cprint("EXECUTANDO ALGORITMOS")
        cprint("=" * 50)

        # Configurar workers internos baseado no número de CPUs e workers externos
        import multiprocessing

        cpu_count = multiprocessing.cpu_count()
        external_workers = args.workers

        # Verificar se algum algoritmo suporta paralelismo interno
        has_internal_parallel = any(
            getattr(global_registry.get(alg_name), "supports_internal_parallel", False)
            for alg_name in viable_algs
            if alg_name in global_registry
        )

        if has_internal_parallel:
            # Se há paralelismo interno, use apenas 1 worker externo
            external_workers = 1
            internal_workers = max(1, cpu_count // 1)
        else:
            # Se não há paralelismo interno, use todos os workers externos
            internal_workers = max(1, cpu_count // external_workers)

        # Configurar variável de ambiente para workers internos
        os.environ["INTERNAL_WORKERS"] = str(internal_workers)

        cprint("🔧 Configuração de paralelismo:")
        cprint(f"   - CPUs disponíveis: {cpu_count}")
        cprint(f"   - Workers externos: {external_workers}")
        cprint(f"   - Workers internos: {internal_workers}")
        cprint(f"   - Paralelismo interno detectado: {has_internal_parallel}")

        # Decidir se usar execução paralela ou sequencial
        use_parallel = len(viable_algs) > 1 and num_execs == 1

        if use_parallel:
            cprint(f"🚀 Usando execução paralela para {len(viable_algs)} algoritmos")
            from src.core.exec.parallel_runner import execute_algorithms_parallel

            parallel_results = execute_algorithms_parallel(
                algorithm_names=viable_algs,
                seqs=seqs,
                alphabet=alphabet,
                console=console_manager,
                max_workers=external_workers,
                timeout=timeout,
            )

            formatter = ResultsFormatter()
            results = {}

            for alg_name, alg_data in parallel_results.items():
                # Converter resultado do formato paralelo para o formato esperado
                executions = [alg_data]  # Resultado paralelo já é uma execução
                executions[0]["seed"] = seed

                formatter.add_algorithm_results(alg_name, executions)

                if "distancia" in alg_data and alg_data["distancia"] != float("inf"):
                    dist_base = params.get("distancia_string_base", "-")
                    results[alg_name] = {
                        "dist": alg_data["distancia"],
                        "dist_base": dist_base,
                        "time": alg_data["tempo"],
                    }
                else:
                    results[alg_name] = {
                        "dist": "-",
                        "dist_base": "-",
                        "time": "-",
                    }
        else:
            cprint(f"🔄 Usando execução sequencial para {len(viable_algs)} algoritmos")
            from src.core.exec.runner import execute_algorithm_runs

            formatter = ResultsFormatter()
            results = {}

            for alg_name in viable_algs:
                if alg_name not in global_registry:
                    cprint(f"ERRO: Algoritmo '{alg_name}' não encontrado!")
                    continue

                # Log apenas início simplificado
                logging.debug(f"[ALG_EXEC] Iniciando {alg_name}")
                AlgClass = global_registry[alg_name]

                executions = execute_algorithm_runs(
                    alg_name, AlgClass, seqs, alphabet, num_execs, None, console_manager, timeout
                )

                # Log resumido das execuções
                logging.debug(f"[ALG_EXEC] {alg_name} concluído: {len(executions)} execuções")
                for _, exec_data in enumerate(executions):
                    # Não calcular mais distancia_string_base aqui, apenas usar seed
                    exec_data["seed"] = seed

                formatter.add_algorithm_results(alg_name, executions)
                valid_results = [e for e in executions if "distancia" in e and e["distancia"] != float("inf")]
                if valid_results:
                    best_exec = min(valid_results, key=lambda e: e["distancia"])
                    logging.debug(f"[ALG_EXEC] {alg_name} melhor: dist={best_exec['distancia']}")

                    # Adicionar distância da string base ao resultado
                    dist_base = params.get("distancia_string_base", "-")

                    results[alg_name] = {
                        "dist": best_exec["distancia"],
                        "dist_base": dist_base,
                        "time": best_exec["tempo"],
                    }
                else:
                    error_exec = next((e for e in executions if "erro" in e), executions[0])
                    logging.debug(f"[ALG_EXEC] {alg_name} sem resultados válidos")
                    results[alg_name] = {
                        "dist": "-",
                        "dist_base": "-",
                        "time": error_exec["tempo"],
                        "warn": error_exec.get("erro", "Erro desconhecido"),
                    }

        # Exibir resumo dos resultados
        print_quick_summary(results, console_manager)

        cprint("\n📄 Gerando relatório detalhado...")
        # Adicionar informações básicas ao formatter para o relatório
        if hasattr(formatter, "__dict__"):
            # Captura todas as strings base e suas distâncias
            base_strings_info = []
            for _, execs in formatter.results.items():
                for exec_data in execs:
                    if exec_data.get("melhor_string"):
                        base_strings_info.append(
                            {
                                "base_string": exec_data["melhor_string"],
                                "distancia_string_base": params.get("distancia_string_base", "-"),
                            }
                        )
            formatter.extra_info = {
                "seed": seed,
                "params": params,
                "dataset_strings": seqs,
                "base_strings_info": base_strings_info,
            }
            logging.debug("[main] formatter configurado")

        results_dir = "results"
        txt_path = os.path.join(results_dir, f"{base_name}.txt")
        csv_path = os.path.join(results_dir, f"{base_name}.csv")

        formatter.save_detailed_report(txt_path)

        # Salvar resultados detalhados em CSV
        formatter.export_to_csv(csv_path)
        cprint(f"📄 Resultados detalhados exportados para CSV: {csv_path}")

        # Sucesso - retornar 0 explicitamente
        if silent:
            sys.exit(0)

    except Exception as e:
        if not silent:
            console_manager.print(f"\nERRO FATAL: {e}")
            traceback.print_exc()
        else:
            logging.exception("Erro fatal durante execução", exc_info=e)
        sys.exit(1)
    finally:
        # Limpar console curses se necessário
        if hasattr(console_manager, "cleanup") and callable(getattr(console_manager, "cleanup", None)):
            console_manager.cleanup()


def generate_dataset_automated():
    """Gera dataset sintético com valores padrão para testes automatizados.

    Utiliza parâmetros fixos para garantir reprodutibilidade em testes.

    Returns:
        tuple: Lista de strings do dataset e dicionário de parâmetros utilizados.
    """
    from src.datasets.dataset_synthetic import SYNTHETIC_DEFAULTS

    n = SYNTHETIC_DEFAULTS["n"]
    L = SYNTHETIC_DEFAULTS["L"]
    alphabet = SYNTHETIC_DEFAULTS["alphabet"]
    noise = SYNTHETIC_DEFAULTS["noise"]
    fully_random = False
    seed = 42
    params = {
        "n": n,
        "L": L,
        "alphabet": alphabet,
        "noise": noise,
        "fully_random": fully_random,
        "seed": seed,
    }
    import random

    rng = random.Random(seed)
    base_string = "".join(rng.choices(alphabet, k=L))
    data = []
    for _ in range(n):
        s = list(base_string)
        num_mut = int(round(noise * L))
        mut_pos = rng.sample(range(L), num_mut) if num_mut > 0 else []
        for pos in mut_pos:
            orig = s[pos]
            alt = rng.choice([c for c in alphabet if c != orig])
            s[pos] = alt
        new_s = "".join(s)
        data.append(new_s)
    return data, params


def run_optimization_workflow():
    """Executa o workflow de otimização de hiperparâmetros."""
    from src.datasets.dataset_synthetic import generate_dataset
    from src.optimization.optuna_optimizer import optimize_algorithm

    console.print("\n=== Otimização de Hiperparâmetros ===")

    # Selecionar algoritmo
    algorithm_name = select_optimization_algorithm()
    if not algorithm_name:
        console.print("❌ Nenhum algoritmo selecionado.")
        return

    # Configurar parâmetros
    config = configure_optimization_params()

    # Gerar dataset para otimização
    console.print("\n📊 Gerando dataset para otimização...")
    seqs, _ = generate_dataset(silent=True)

    console.print(f"✅ Dataset gerado: {len(seqs)} sequências de tamanho {len(seqs[0])}")

    # Executar otimização
    alphabet = "".join(sorted(set("".join(seqs))))

    try:
        result = optimize_algorithm(
            algorithm_name=algorithm_name,
            sequences=seqs,
            alphabet=alphabet,
            n_trials=config["n_trials"],
            timeout_per_trial=config["timeout"],
            show_progress=True,
        )

        console.print("\n✅ Otimização concluída!")
        console.print(f"Melhor valor: {result.best_value}")
        console.print(f"Melhores parâmetros: {result.best_params}")

        # Salvar visualizações se solicitado
        if config["save_plots"]:
            from src.optimization.visualization import OptimizationVisualizer

            visualizer = OptimizationVisualizer(result)
            plots_dir = os.path.join(RESULTS_DIR, "optimization_plots")
            os.makedirs(plots_dir, exist_ok=True)

            # Salvar gráficos
            history_path = os.path.join(plots_dir, f"{algorithm_name}_history.png")
            importance_path = os.path.join(plots_dir, f"{algorithm_name}_importance.png")

            visualizer.plot_optimization_history(save_path=history_path)
            visualizer.plot_parameter_importance(save_path=importance_path)

            console.print(f"📊 Gráficos salvos em: {plots_dir}")

        # Salvar resultados
        import json

        results_path = os.path.join(
            RESULTS_DIR, f"optimization_{algorithm_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )
        with open(results_path, "w") as f:
            json.dump(
                {
                    "best_params": result.best_params,
                    "best_value": result.best_value,
                    "n_trials": result.n_trials,
                    "study_name": result.study_name,
                    "optimization_time": result.optimization_time,
                    "all_trials": result.all_trials,
                },
                f,
                indent=2,
                default=str,
            )

        console.print(f"💾 Resultados salvos em: {results_path}")

    except Exception as e:
        console.print(f"❌ Erro durante otimização: {e}")
        logging.exception("Erro na otimização", exc_info=e)


def run_sensitivity_workflow():
    """Executa o workflow de análise de sensibilidade."""
    from src.datasets.dataset_synthetic import generate_dataset
    from src.optimization.sensitivity_analyzer import analyze_algorithm_sensitivity

    console.print("\n=== Análise de Sensibilidade ===")

    # Selecionar algoritmo
    algorithm_name = select_sensitivity_algorithm()
    if not algorithm_name:
        console.print("❌ Nenhum algoritmo selecionado.")
        return

    # Configurar parâmetros
    config = configure_sensitivity_params()

    # Gerar dataset para análise
    console.print("\n📊 Gerando dataset para análise...")
    seqs, _ = generate_dataset(silent=True)

    console.print(f"✅ Dataset gerado: {len(seqs)} sequências de tamanho {len(seqs[0])}")

    # Executar análise de sensibilidade
    alphabet = "".join(sorted(set("".join(seqs))))

    try:
        result = analyze_algorithm_sensitivity(
            algorithm_name=algorithm_name,
            sequences=seqs,
            alphabet=alphabet,
            n_samples=config["n_samples"],
            timeout_per_sample=config.get("timeout", 60),
            show_progress=True,
        )

        console.print("\n✅ Análise de sensibilidade concluída!")
        console.print(f"Parâmetros analisados: {len(result.parameter_names)}")

        # Mostrar principais parâmetros sensíveis
        important_params = result.get_most_important_parameters(n=5)
        console.print("\n📈 Parâmetros mais sensíveis:")
        for param, value in important_params:
            console.print(f"  • {param}: {value:.4f}")

        # Salvar visualizações se solicitado
        if config["save_plots"]:
            from src.optimization.visualization import SensitivityVisualizer

            visualizer = SensitivityVisualizer(result)
            plots_dir = os.path.join(RESULTS_DIR, "sensitivity_plots")
            os.makedirs(plots_dir, exist_ok=True)

            # Salvar gráficos
            sensitivity_path = os.path.join(plots_dir, f"{algorithm_name}_sensitivity.png")
            visualizer.plot_sensitivity_indices(save_path=sensitivity_path)

            console.print(f"📊 Gráficos salvos em: {plots_dir}")

        # Salvar resultados
        import json

        results_path = os.path.join(
            RESULTS_DIR, f"sensitivity_{algorithm_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )
        with open(results_path, "w") as f:
            json.dump(
                {
                    "method": result.method,
                    "parameter_names": result.parameter_names,
                    "first_order": result.first_order,
                    "total_order": result.total_order,
                    "second_order": result.second_order,
                    "confidence_intervals": result.confidence_intervals,
                    "n_samples": result.n_samples,
                    "analysis_time": result.analysis_time,
                },
                f,
                indent=2,
                default=str,
            )

        console.print(f"💾 Resultados salvos em: {results_path}")

    except Exception as e:
        console.print(f"❌ Erro durante análise de sensibilidade: {e}")
        logging.exception("Erro na análise de sensibilidade", exc_info=e)
