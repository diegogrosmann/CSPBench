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

from datetime import datetime
import uuid
import sys
import traceback
import signal
import os
import argparse

from src.menu import menu, select_algorithms
from src.runner import execute_algorithm_runs
from src.logging_utils import setup_logging
from src.report_utils import print_quick_summary, save_detailed_report
from algorithms.base import global_registry
from datasets.dataset_utils import ask_save_dataset
from src.results_formatter import ResultsFormatter
from src.console_manager import console
from utils.resource_monitor import get_safe_memory_limit, check_algorithm_feasibility, ResourceMonitor
from utils.config import safe_input, ALGORITHM_TIMEOUT
from utils.distance import max_distance, hamming_distance
import logging

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
    parser.add_argument('--silent', action='store_true', help='Executa em modo silencioso (sem prints interativos)')
    parser.add_argument('--dataset', type=str, choices=['synthetic', 'file', 'entrez', 'batch'], help='Fonte do dataset')
    parser.add_argument('--algorithms', type=str, nargs='+', help='Algoritmos a executar (nomes separados por espaço)')
    parser.add_argument('--num-execs', type=int, help='Número de execuções por algoritmo')
    parser.add_argument('--timeout', type=int, help='Timeout por execução (segundos)')
    args = parser.parse_args()

    silent = args.silent
    def silent_print(*a, **k):
        pass
    p = print if not silent else silent_print
    cprint = console.print if not silent else silent_print

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
                args.dataset = 'synthetic'
            if not args.algorithms:
                args.algorithms = ['BLF-GA']
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
            if args.dataset == 'synthetic':
                from datasets.dataset_synthetic import generate_dataset
                seqs, p = generate_dataset(silent=silent)
                logging.debug(f"[main] Parâmetros do dataset: {p}")
                params = {'dataset_source': '1'}
                params.update(p)
            elif args.dataset == 'file':
                from datasets.dataset_file import load_dataset
                seqs, p = load_dataset(silent=silent)
                params = {'dataset_source': '2'}
                params.update(p)
            elif args.dataset == 'entrez':
                from datasets.dataset_entrez import fetch_dataset
                seqs, p = fetch_dataset()
                params = {'dataset_source': '3'}
                params.update(p)
            elif args.dataset == 'batch':
                from src.batch_executor import BatchExecutor, select_batch_config
                config_file = ''
                executor = BatchExecutor(config_file)
                batch_result = executor.execute_batch()
                cprint(f"\n✅ Execução em lote concluída!")
                cprint(f"Tempo total: {batch_result['tempo_total']:.1f}s")
                cprint(f"Taxa de sucesso: {batch_result['resumo']['taxa_sucesso']:.1f}%")
                from src.export_csv_batch import export_batch_json_to_csv
                batch_dir = executor.results_dir
                import os as os_batch
                json_path = os_batch.path.join(batch_dir, "batch_results.json")
                csv_path = os_batch.path.join(batch_dir, "batch_results.csv")
                if os_batch.path.exists(json_path):
                    export_batch_json_to_csv(json_path, csv_path)
                    cprint(f"📄 Resultados detalhados do batch exportados para CSV: {csv_path}")
                return
            else:
                cprint("❌ Fonte de dataset inválida.")
                return
        else:
            # Fluxo interativo
            choice = menu()
            params = {'dataset_source': choice}
            
            if choice == '4':
                from src.batch_executor import BatchExecutor, select_batch_config
                config_file = select_batch_config()
                if not config_file:
                    cprint("❌ Nenhum arquivo de configuração selecionado.")
                    return
                try:
                    executor = BatchExecutor(config_file)
                    batch_result = executor.execute_batch()
                    cprint(f"\n✅ Execução em lote concluída!")
                    cprint(f"Tempo total: {batch_result['tempo_total']:.1f}s")
                    cprint(f"Taxa de sucesso: {batch_result['resumo']['taxa_sucesso']:.1f}%")
                    from src.export_csv_batch import export_batch_json_to_csv
                    batch_dir = executor.results_dir
                    import os as os_batch
                    json_path = os_batch.path.join(batch_dir, "batch_results.json")
                    csv_path = os_batch.path.join(batch_dir, "batch_results.csv")
                    if os_batch.path.exists(json_path):
                        export_batch_json_to_csv(json_path, csv_path)
                        cprint(f"📄 Resultados detalhados do batch exportados para CSV: {csv_path}")
                except Exception as e:
                    cprint(f"❌ Erro na execução em lote: {e}")
                    return
            elif choice in ['1', '2', '3']:
                try:
                    if choice == '1':
                        from datasets.dataset_synthetic import generate_dataset
                        seqs, p = generate_dataset(silent=silent)
                        logging.debug(f"[main] Parâmetros do dataset: {p}")
                        params.update(p)
                        if not silent:
                            ask_save_dataset(seqs, "synthetic", p)
                    elif choice == '2':
                        from datasets.dataset_file import load_dataset
                        seqs, p = load_dataset(silent=silent)
                        params.update(p)
                    elif choice == '3':
                        from datasets.dataset_entrez import fetch_dataset
                        seqs, p = fetch_dataset()
                        params.update(p)
                        ask_save_dataset(seqs, "entrez", p)
                except Exception as exc:
                    cprint(f"Erro: {exc}")
                    logging.exception("Erro ao carregar dataset", exc_info=exc)
                    return

        # Após carregar o dataset, extrair informações extras se disponíveis
        seed = params.get('seed')
        logging.debug(f"[main] seed: {seed}")

        # Exibir distância da string base se disponível
        distancia_base = params.get('distancia_string_base')
        if distancia_base is not None:
            cprint(f"Distância da string base: {distancia_base}")

        if not seqs:
            cprint("Nenhuma sequência lida.")
            return

        alphabet = ''.join(sorted(set(''.join(seqs))))
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
        cprint("\n" + "="*50)
        cprint("EXECUTANDO ALGORITMOS")
        cprint("="*50)
        
        from src.runner import Spinner, execute_algorithm_runs

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
                alg_name, AlgClass, seqs, alphabet, num_execs, None, console, timeout
            )
            
            # Log resumido das execuções
            logging.debug(f"[ALG_EXEC] {alg_name} concluído: {len(executions)} execuções")
            for i, exec_data in enumerate(executions):
                # Não calcular mais distancia_string_base aqui, apenas usar seed
                exec_data['seed'] = seed
            
            formatter.add_algorithm_results(alg_name, executions)
            valid_results = [e for e in executions if 'distancia' in e and e['distancia'] != float('inf')]
            if valid_results:
                best_exec = min(valid_results, key=lambda e: e['distancia'])
                logging.debug(f"[ALG_EXEC] {alg_name} melhor: dist={best_exec['distancia']}")
                
                # Adicionar distância da string base ao resultado
                dist_base = params.get('distancia_string_base', '-')
                
                results[alg_name] = {
                    'dist': best_exec['distancia'],
                    'dist_base': dist_base,
                    'time': best_exec['tempo']
                }
            else:
                error_exec = next((e for e in executions if 'erro' in e), executions[0])
                logging.debug(f"[ALG_EXEC] {alg_name} sem resultados válidos")
                results[alg_name] = {
                    'dist': '-',
                    'dist_base': '-',
                    'time': error_exec['tempo'],
                    'warn': error_exec.get('erro', 'Erro desconhecido')
                }

        # Exibir resumo dos resultados
        print_quick_summary(results, console)

        cprint(f"\n📄 Gerando relatório detalhado...")
        # Adicionar informações básicas ao formatter para o relatório
        if hasattr(formatter, '__dict__'):
            # Captura todas as strings base e suas distâncias
            base_strings_info = []
            for alg_name, execs in formatter.results.items():
                for exec_data in execs:
                    if exec_data.get('melhor_string'):
                        base_strings_info.append({
                            'base_string': exec_data['melhor_string'],
                            'distancia_string_base': params.get('distancia_string_base', '-')
                        })
            formatter.extra_info = {
                'seed': seed,
                'params': params,
                'dataset_strings': seqs,
                'base_strings_info': base_strings_info
            }
            logging.debug(f"[main] formatter configurado")
        
        results_dir = "results"
        txt_path = os.path.join(results_dir, f"{base_name}.txt")
        csv_path = os.path.join(results_dir, f"{base_name}.csv")

        save_detailed_report(formatter, txt_path)

        # Salvar resultados detalhados em CSV
        from src.export_csv import export_results_to_csv
        export_results_to_csv(formatter, csv_path)
        cprint(f"📄 Resultados detalhados exportados para CSV: {csv_path}")

    except Exception as e:
        console.print(f"\nERRO FATAL: {e}")
        traceback.print_exc()
        sys.exit(1)

def generate_dataset_automated():
    """Gera dataset sintético com valores padrão para testes automatizados.

    Utiliza parâmetros fixos para garantir reprodutibilidade em testes.

    Returns:
        tuple: Lista de strings do dataset e dicionário de parâmetros utilizados.
    """
    from datasets.dataset_synthetic import SYNTHETIC_DEFAULTS, generate_dataset
    n = SYNTHETIC_DEFAULTS['n']
    L = SYNTHETIC_DEFAULTS['L']
    alphabet = SYNTHETIC_DEFAULTS['alphabet']
    noise = SYNTHETIC_DEFAULTS['noise']
    fully_random = False
    seed = 42
    params = {'n': n, 'L': L, 'alphabet': alphabet, 'noise': noise, 'fully_random': fully_random, 'seed': seed}
    import random
    rng = random.Random(seed)
    base_string = ''.join(rng.choices(alphabet, k=L))
    data = []
    for _ in range(n):
        s = list(base_string)
        num_mut = int(round(noise * L))
        mut_pos = rng.sample(range(L), num_mut) if num_mut > 0 else []
        for pos in mut_pos:
            orig = s[pos]
            alt = rng.choice([c for c in alphabet if c != orig])
            s[pos] = alt
        new_s = ''.join(s)
        data.append(new_s)
    return data, params

if __name__ == "__main__":
    main()