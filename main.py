"""
Arquivo principal de execução da aplicação Closest String Problem (CSP).

Fluxo principal:
    1. Seleção e leitura/geração do dataset.
    2. Seleção dos algoritmos a executar.
    3. Execução do Baseline e dos algoritmos selecionados.
    4. Exibição e salvamento dos resultados.

Funções:
    main(): Executa o fluxo principal da aplicação.
"""

from datetime import datetime
import uuid
import sys
import traceback
import signal

from src.menu import menu, select_algorithms
from src.runner import execute_algorithm_runs
from src.logging_utils import setup_logging
from src.report_utils import print_quick_summary, save_detailed_report
from algorithms.baseline.implementation import greedy_consensus, max_distance
from algorithms.base import global_registry
from datasets.dataset_utils import ask_save_dataset
from src.results_formatter import ResultsFormatter
from src.console_manager import console
from utils.resource_monitor import get_safe_memory_limit, check_algorithm_feasibility, ResourceMonitor
from utils.config import safe_input, ALGORITHM_TIMEOUT
import logging

def signal_handler(signum, frame):
    """Handler para sinais de interrupção."""
    print("\n\nOperação cancelada pelo usuário. Encerrando.")
    sys.exit(0)

def main():
    # Configurar handlers de sinal para saída limpa
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    try:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        uid = uuid.uuid4().hex[:8]
        base_name = f"{ts}_{uid}"
        setup_logging(base_name)

        # Dataset
        choice = menu()
        
        # Nova opção: Execução em lote
        if choice == '4':
            from src.batch_executor import BatchExecutor, select_batch_config
            
            # Menu de seleção de arquivo de configuração
            config_file = select_batch_config()
            
            if not config_file:
                console.print("❌ Nenhum arquivo de configuração selecionado.")
                return
            
            try:
                executor = BatchExecutor(config_file)
                batch_result = executor.execute_batch()
                
                console.print(f"\n✅ Execução em lote concluída!")
                console.print(f"Tempo total: {batch_result['tempo_total']:.1f}s")
                console.print(f"Taxa de sucesso: {batch_result['resumo']['taxa_sucesso']:.1f}%")
                
            except Exception as e:
                console.print(f"❌ Erro na execução em lote: {e}")
                logging.exception("Erro na execução em lote")
                return
            
            return  # Sair após execução em lote
        
        # Fluxo normal para opções 1, 2, 3
        params = {'dataset_source': choice}
        seqs = []

        try:
            if choice == '1':
                from datasets.dataset_synthetic import generate_dataset
                seqs, p = generate_dataset()
                params.update(p)
                ask_save_dataset(seqs, "synthetic", p)
            elif choice == '2':
                from datasets.dataset_file import load_dataset
                seqs, p = load_dataset()
                params.update(p)
            else:
                from datasets.dataset_entrez import fetch_dataset
                seqs, p = fetch_dataset()
                params.update(p)
                ask_save_dataset(seqs, "entrez", p)
        except Exception as exc:
            console.print(f"Erro: {exc}")
            logging.exception("Erro ao carregar dataset", exc_info=exc)
            return

        if not seqs:
            console.print("Nenhuma sequência lida.")
            return

        alphabet = ''.join(sorted(set(''.join(seqs))))
        n, L = len(seqs), len(seqs[0])
        console.print(f"\nDataset: n={n}, L={L}, |Σ|={len(alphabet)}")
        
        # Verificação de recursos do sistema
        safe_memory = get_safe_memory_limit()
        console.print(f"Limite seguro de memória: {safe_memory:.1f}%")
        
        # Algoritmos
        algs = select_algorithms()
        if not algs:
            console.print("Nenhum algoritmo selecionado além do Baseline.")

        # Verificar viabilidade dos algoritmos selecionados
        viable_algs = []
        for alg_name in algs:
            is_viable, msg = check_algorithm_feasibility(n, L, alg_name)
            if is_viable:
                viable_algs.append(alg_name)
                console.print(f"✓ {alg_name}: {msg}")
            else:
                console.print(f"⚠ {alg_name}: {msg} (será pulado)")
        
        if not viable_algs:
            console.print("Nenhum algoritmo viável. Executando apenas Baseline.")
        
        runs = safe_input("\nNº execuções p/ algoritmo [3]: ")
        num_execs = int(runs) if runs.isdigit() and int(runs) > 0 else 3

        # Configurar timeout com recomendações baseadas nos algoritmos
        default_timeout = ALGORITHM_TIMEOUT
        if 'DP-CSP' in viable_algs and n >= 8:
            default_timeout = max(ALGORITHM_TIMEOUT, 120)  # Mínimo 2 minutos para DP-CSP complexo
            console.print("⚠ DP-CSP detectado em dataset complexo - timeout mínimo aumentado")
        
        timeout_input = safe_input(f"\nTimeout por execução em segundos [{default_timeout}]: ")
        timeout = int(timeout_input) if timeout_input.isdigit() and int(timeout_input) > 0 else default_timeout
        
        console.print(f"Timeout configurado: {timeout}s por execução")
        
        # Baseline
        console.print("\n" + "="*50)
        console.print("EXECUTANDO ALGORITMOS")
        console.print("="*50)
        
        from src.runner import Spinner
        baseline_spinner = Spinner("Baseline", console)
        baseline_spinner.start()
        
        import time
        t0 = time.time()
        baseline_center = greedy_consensus(seqs, alphabet)
        baseline_val = max_distance(baseline_center, seqs)
        baseline_time = time.time() - t0
        
        baseline_spinner.stop()
        console.print(f"Baseline                 ... dist={baseline_val}, tempo={baseline_time:.3f}s")

        formatter = ResultsFormatter()
        baseline_executions = [{
            'tempo': baseline_time,
            'iteracoes': 1,
            'distancia': baseline_val,
            'melhor_string': baseline_center,
            'gap': 0.0
        }]
        formatter.add_algorithm_results("Baseline", baseline_executions)

        # Execução dos algoritmos (usando viable_algs em vez de algs)
        results = {}
        for alg_name in viable_algs:
            if alg_name not in global_registry:
                console.print(f"ERRO: Algoritmo '{alg_name}' não encontrado!")
                continue
            AlgClass = global_registry[alg_name]
            executions = execute_algorithm_runs(
                alg_name, AlgClass, seqs, alphabet, num_execs, baseline_val, console, timeout
            )
            formatter.add_algorithm_results(alg_name, executions)
            
            valid_results = [e for e in executions if 'distancia' in e and e['distancia'] != float('inf')]
            if valid_results:
                best_exec = min(valid_results, key=lambda e: e['distancia'])
                results[alg_name] = {
                    'dist': best_exec['distancia'],
                    'time': best_exec['tempo'],
                    'gap': best_exec['gap']
                }
            else:
                error_exec = next((e for e in executions if 'erro' in e), executions[0])
                results[alg_name] = {
                    'dist': '-',
                    'time': error_exec['tempo'],
                    'gap': '-',
                    'warn': error_exec.get('erro', 'Erro desconhecido')
                }

        print_quick_summary(baseline_val, baseline_time, results, console)

        console.print(f"\n📄 Gerando relatório detalhado...")
        save_detailed_report(formatter, f"{base_name}.txt")

    except Exception as e:
        console.print(f"\nERRO FATAL: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
