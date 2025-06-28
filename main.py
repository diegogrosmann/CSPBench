# main.py
"""
main.py
=======

Interface interativa para resolver o Closest String Problem (CSP)
usando algoritmos registrados dinamicamente com relatórios detalhados.

Execute:
$ python main.py
"""

from __future__ import annotations

import algorithms  # This will trigger auto-import of all algorithm packages
import logging
import sys
import time
import traceback
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
import threading
import itertools
import statistics

from algorithms.baseline.implementation import greedy_consensus, max_distance
from algorithms.base import global_registry
from utils.config import DEBUG_DEFAULT
from datasets.dataset_utils import ask_save_dataset
from src.results_formatter import ResultsFormatter


# ---------------------------------------------------------------------
# Utilitários
# ---------------------------------------------------------------------

class Spinner:
    def __init__(self, prefix: str):
        self.prefix = prefix
        self.stop_event = threading.Event()
        self.thread = None
        self.spinner = itertools.cycle(['   ', '.  ', '.. ', '...'])
    
    def start(self):
        if self.thread and self.thread.is_alive():
            return
        self.stop_event.clear()
        self.thread = threading.Thread(target=self._animate)
        self.thread.daemon = True
        self.thread.start()
    
    def _animate(self):
        while not self.stop_event.is_set():
            print(f"\r{self.prefix:<25s}...{next(self.spinner)}", end="", flush=True)
            time.sleep(0.3)
    
    def stop(self):
        self.stop_event.set()
        if self.thread and self.thread.is_alive():
            self.thread.join(timeout=1.0)
        print(f"\r{self.prefix:<25s}...   ", end="", flush=True)

def _gap(baseline: int, val: Optional[int]) -> float:
    """% de melhora – seguro p/ val=None."""
    return 100 * (baseline - val) / baseline if (baseline and val is not None) else 0.0


# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------

def setup_logging(base_name: str) -> None:
    debug_input = input("Habilitar modo debug? [s/N]: ").strip().lower()
    debug_mode = debug_input if debug_input else DEBUG_DEFAULT
    
    Path("logs").mkdir(exist_ok=True)
    log_file = Path("logs") / f"{base_name}.log"
    
    if debug_mode == 's':
        level = logging.DEBUG
        print(f"Debug ON -> {log_file}")
    else:
        level = logging.WARNING
        
    logging.basicConfig(
        level=level,
        filename=str(log_file),
        filemode='w',
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        force=True
    )
    
    for handler in logging.root.handlers[:]:
        if isinstance(handler, logging.StreamHandler) and handler.stream == sys.stdout:
            logging.root.removeHandler(handler)


# ---------------------------------------------------------------------
# Menus
# ---------------------------------------------------------------------

def menu() -> str:
    print("\n=== Closest String Problem ===")
    print("1) Gerar dataset sintético")
    print("2) Carregar dataset de arquivo")
    print("3) Baixar dataset via NCBI")
    while True:
        c = input("Escolha [1/2/3]: ").strip()
        if c in {'1', '2', '3'}:
            return c
        print("Inválido.")


def select_algorithms() -> List[str]:
    """Lista dinâmica de algoritmos registrados (excluindo Baseline)."""
    all_algs = list(global_registry.keys())
    algs = [name for name in all_algs if name != "Baseline"]
    
    print("\nAlgoritmos disponíveis:")
    print(" 0) Executar todos")
    for idx, name in enumerate(algs, 1):
        print(f" {idx}) {name}")
    
    raw = input("Escolha (ex.: 1,3 ou 0 para todos) [padrão 1]: ").strip()
    if not raw:
        return [algs[0]] if algs else []
    
    if raw == '0':
        return algs
    
    selected = []
    for part in raw.split(','):
        if part.strip().isdigit():
            i = int(part)
            if 1 <= i <= len(algs):
                selected.append(algs[i-1])
    return selected


# ---------------------------------------------------------------------
# Execução de algoritmos com resultados detalhados
# ---------------------------------------------------------------------

def execute_algorithm_runs(alg_name: str, AlgClass, seqs: List[str], alphabet: str, 
                          num_execs: int, baseline_val: int) -> List[Dict[str, Any]]:
    """Executa múltiplas rodadas de um algoritmo e coleta resultados detalhados."""
    
    is_deterministic = getattr(AlgClass, 'is_deterministic', False)
    actual_execs = 1 if is_deterministic else num_execs
    
    executions = []
    
    for i in range(actual_execs):
        exec_prefix = f"{alg_name}"
        if actual_execs > 1:
            exec_prefix += f" ({i+1}/{actual_execs})"
        
        spinner = Spinner(exec_prefix)
        result_holder = {}
        exc_holder = {}

        def run_algorithm():
            try:
                instance = AlgClass(seqs, alphabet)
                if hasattr(instance, 'set_progress_callback'):
                    def progress_callback(msg: str):
                        spinner.stop()
                        print(f"\r{exec_prefix:<25s}... {msg:<40s}", end="", flush=True)
                    instance.set_progress_callback(progress_callback)
                
                center, val = instance.run()
                result_holder['center'] = center
                result_holder['val'] = val
                
                # Interface padronizada para capturar dados adicionais
                if hasattr(instance, 'geracao'):
                    result_holder['iteracoes'] = instance.geracao
                elif hasattr(instance, 'iterations'):
                    result_holder['iteracoes'] = instance.iterations
                elif hasattr(instance, 'num_iteracoes'):
                    result_holder['iteracoes'] = instance.num_iteracoes
                else:
                    result_holder['iteracoes'] = 1
                    
            except Exception as e:
                exc_holder['exc'] = e

        t0 = time.time()
        algo_thread = threading.Thread(target=run_algorithm)
        spinner.start()
        algo_thread.start()
        
        while algo_thread.is_alive():
            time.sleep(0.1)
        
        spinner.stop()
        tempo_execucao = time.time() - t0

        print(f"\r{'':<80s}", end="")

        if 'exc' in exc_holder:
            error_msg = str(exc_holder['exc'])[:50]
            print(f"\r{exec_prefix:<25s}... ERRO: {error_msg}...")
            
            executions.append({
                'tempo': tempo_execucao,
                'iteracoes': 0,
                'distancia': float('inf'),
                'melhor_string': '',
                'erro': error_msg
            })
        else:
            center = result_holder['center']
            val = result_holder['val']
            iteracoes = result_holder.get('iteracoes', 1)
            
            gap = _gap(baseline_val, val)
            print(f"\r{exec_prefix:<25s}... dist={val}, tempo={tempo_execucao:.3f}s, gap={gap:.1f}%")
            
            executions.append({
                'tempo': tempo_execucao,
                'iteracoes': iteracoes,
                'distancia': val,
                'melhor_string': center,
                'gap': gap
            })
    
    return executions


# ---------------------------------------------------------------------
# Salvamento de resultados
# ---------------------------------------------------------------------

def save_results(params: Dict[str, Any], baseline_center: str, baseline_val: int, 
                baseline_time: float, results: Dict[str, Dict[str, Any]], 
                out_path: Path) -> None:
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)

    with out_path.open('w', encoding='utf-8') as f:
        f.write("=== Resultados Closest String Problem ===\n\n")
        f.write(f"n={len(params['seqs'])}  L={len(params['seqs'][0])}\n\n")

        f.write("## Baseline\n")
        f.write(f"dist={baseline_val}  time={baseline_time:.3f}s\n")
        f.write(f"string={baseline_center}\n\n")

        for alg, res in results.items():
            f.write(f"## {alg}\n")
            f.write(f"dist={res.get('dist','-')}  time={res.get('time',0):.3f}s  gap={res.get('gap',0):.1f}%\n")
            if res.get('warn'):
                f.write(f"(Aviso: {res['warn']})\n")

    print(f"Resultados salvos em {out_path.resolve()}")
    print(f"Logs salvos em logs/{out_path.stem}.log")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> None:
    try:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        uid = uuid.uuid4().hex[:8]
        base_name = f"{ts}_{uid}"
        setup_logging(base_name)

        # ------------ Dataset -------------------------------------------------
        choice = menu()
        params: Dict[str, Any] = {'dataset_source': choice}
        seqs: List[str] = []

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
            print(f"Erro: {exc}")
            logging.exception(exc)
            return

        if not seqs:
            print("Nenhuma sequência lida.")
            return
        params['seqs'] = seqs

        alphabet = ''.join(sorted(set(''.join(seqs))))
        print(f"\nDataset: n={len(seqs)}, L={len(seqs[0])}, |Σ|={len(alphabet)}")

        # ------------ Algoritmos ---------------------------------------------
        algs = select_algorithms()
        if not algs:
            print("Nenhum algoritmo selecionado além do Baseline.")
        
        runs = input("\nNº execuções p/ algoritmo [3]: ").strip()
        num_execs = int(runs) if runs.isdigit() and int(runs) > 0 else 3

        # ------------ Baseline (sempre executado) ----------------------------
        print("\n" + "="*50)
        print("EXECUTANDO ALGORITMOS")
        print("="*50)
        
        print(f"Baseline...", end="", flush=True)
        t0 = time.time()
        baseline_center = greedy_consensus(seqs, alphabet)
        baseline_val = max_distance(baseline_center, seqs)
        baseline_time = time.time() - t0
        print(f" dist={baseline_val}, tempo={baseline_time:.3f}s")

        # ------------ Inicializar Formatador de Resultados ------------------
        formatter = ResultsFormatter()
        
        # Adicionar baseline ao formatador
        baseline_executions = [{
            'tempo': baseline_time,
            'iteracoes': 1,
            'distancia': baseline_val,
            'melhor_string': baseline_center,
            'gap': 0.0
        }]
        formatter.add_algorithm_results("Baseline", baseline_executions)

        # ------------ Execução Dinâmica de Algoritmos -------------------------
        results = {}
        
        for alg_name in algs:
            if alg_name not in global_registry:
                print(f"ERRO: Algoritmo '{alg_name}' não encontrado!")
                continue
                
            AlgClass = global_registry[alg_name]
            
            # Executar múltiplas rodadas e coletar resultados detalhados
            executions = execute_algorithm_runs(alg_name, AlgClass, seqs, alphabet, num_execs, baseline_val)
            
            # Adicionar ao formatador
            formatter.add_algorithm_results(alg_name, executions)
            
            # Manter melhor resultado para compatibilidade
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

        # ------------ Resumo Rápido ------------------------------------------
        print("\n" + "="*50)
        print("RESUMO RÁPIDO")
        print("="*50)
        print(f"{'Algoritmo':<12} {'Melhor Dist':<12} {'Tempo Médio':<12} {'Gap%':<6}")
        print("-" * 45)
        print(f"{'Baseline':<12} {baseline_val:<12} {baseline_time:<12.3f} {'-':<6}")
        
        for alg_name, res in results.items():
            dist_str = str(res.get('dist', '-'))
            time_str = f"{res.get('time', 0):.3f}"
            gap_str = f"{res.get('gap', 0):.1f}" if isinstance(res.get('gap'), (int, float)) else "-"
            print(f"{alg_name:<12} {dist_str:<12} {time_str:<12} {gap_str:<6}")

        # ------------ Salvar Relatório Detalhado (sem exibir no console) ------
        print(f"\n📄 Gerando relatório detalhado...")
        formatter.save_detailed_report(f"{base_name}.txt")

        # ------------ Salvar Resultados Básicos ------------------------------
        results_path = Path("results") / f"{base_name}.txt"
        save_results(params, baseline_center, baseline_val, baseline_time, results, results_path)
        
    except Exception as e:
        print(f"\nERRO FATAL: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
