# main.py
"""
main.py
=======

Interface 100 % interativa (sem parâmetros CLI) para:

1) Escolher a origem do dataset:
   - Gerar sintético
   - Carregar arquivo
   - Buscar via NCBI (Bio.Entrez)

2) Executar BLF-GA e baseline e mostrar resultados.

Execute:
$ python main.py
"""

import time
from pathlib import Path
import sys
import logging
from typing import Dict, Any

from blf_ga import BLFGA, max_distance
from baseline import greedy_consensus
from config import (DEBUG_DEFAULT, SYNTHETIC_DEFAULTS, FILE_DEFAULTS, 
                   ENTREZ_DEFAULTS, BLFGA_DEFAULTS)
from dataset_utils import ask_save_dataset

# ---- Funções Auxiliares ---- #

def setup_logging():
    """Configura o logging com base na escolha do usuário."""
    debug_mode = input(f"Habilitar modo debug (salvar em debug.log)? [s/N]: ").strip().lower() or DEBUG_DEFAULT
    if debug_mode == 's':
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            filename='debug.log',
            filemode='w'
        )
        print("Modo debug habilitado. Log será salvo em debug.log")
    else:
        logging.disable(logging.CRITICAL)
        print("Modo debug desabilitado.")

def save_results(params, baseline_center, baseline_val, baseline_time, best, best_val, ga_time, gap):
    """Salva os parâmetros e resultados em um arquivo de texto."""
    out_path = Path("results.txt")
    with out_path.open('w', encoding='utf-8') as f:
        f.write("=== Resultados da Execução CSP ===\n\n")
        
        f.write("--- Parâmetros do Dataset ---\n")
        source = params.get('dataset_source')
        if source == '1':
            f.write("Fonte: Sintético\n")
            f.write(f" - Strings (n): {params.get('n')}\n")
            f.write(f" - Comprimento (L): {params.get('L')}\n")
            f.write(f" - Alfabeto: {params.get('alphabet')}\n")
            f.write(f" - Ruído: {params.get('noise')}\n")
        elif source == '2':
            f.write("Fonte: Arquivo\n")
            f.write(f" - Caminho: {params.get('filepath')}\n")
        elif source == '3':
            f.write("Fonte: NCBI\n")
            f.write(f" - Base: {params.get('db')}\n")
            f.write(f" - Termo: {params.get('term')}\n")
            f.write(f" - Registros: {params.get('n')}\n")

        f.write("\n--- Parâmetros BLF-GA ---\n")
        ga_p = params.get('blfga_params', {})
        for key, val in ga_p.items():
            f.write(f" - {key}: {val}\n")

        f.write("\n--- Resultados ---\n")
        f.write(">> Baseline (Consenso Guloso)\n")
        f.write(f" - String: {baseline_center}\n")
        f.write(f" - Distância: {baseline_val}\n")
        f.write(f" - Tempo: {baseline_time:.4f}s\n")

        f.write("\n>> BLF-GA\n")
        f.write(f" - String: {best}\n")
        f.write(f" - Distância: {best_val}\n")
        f.write(f" - Tempo: {ga_time:.4f}s\n")

        f.write("\n--- Comparativo ---\n")
        f.write(f"Melhora vs Baseline: {gap:.2f}%\n")

    print(f"\nResultados e parâmetros salvos em {out_path.resolve()}")


# ---- menu ---- #

def menu() -> str:
    print("\n=== Closest String Problem – BLF‐GA ===")
    print("1) Gerar dataset sintético")
    print("2) Carregar dataset de arquivo (.txt / .fasta)")
    print("3) Baixar dataset via NCBI (Bio.Entrez)")
    while True:
        choice = input("Escolha uma opção [1/2/3]: ").strip()
        if choice in ["1", "2", "3"]:
            logging.debug(f"Você escolheu a opção '{choice}'")
            return choice
        print("Opção inválida. Por favor, escolha 1, 2 ou 3.")

def main():
    setup_logging()
    choice = menu()
    params: Dict[str, Any] = {'dataset_source': choice}
    seqs = []

    try:
        if choice == "1":
            logging.debug("Importando generate_dataset de dataset_synthetic.py")
            from dataset_synthetic import generate_dataset
            seqs, synth_params = generate_dataset()
            params.update(synth_params)
            # Perguntar se deseja salvar o dataset sintético
            ask_save_dataset(seqs, "synthetic", synth_params)
        elif choice == "2":
            logging.debug("Importando load_dataset de dataset_file.py")
            from dataset_file import load_dataset
            seqs, file_params = load_dataset()
            params.update(file_params)
        elif choice == "3":
            logging.debug("Importando fetch_dataset de dataset_entrez.py")
            from dataset_entrez import fetch_dataset
            seqs, entrez_params = fetch_dataset()
            params.update(entrez_params)
            # Perguntar se deseja salvar o dataset do NCBI
            ask_save_dataset(seqs, "entrez", entrez_params)
    except (RuntimeError, FileNotFoundError, ValueError) as e:
        print(f"\nOcorreu um erro: {e}")
        logging.error(f"Erro durante a obtenção do dataset: {e}", exc_info=True)
        return

    if not seqs:
        print("\nNenhuma sequência foi carregada. Encerrando.")
        return

    alphabet = "".join(sorted(set("".join(seqs))))
    logging.debug(f"Sequências lidas: n={len(seqs)}, L={len(seqs[0])}, alfabeto={alphabet}")
    print(f"\nDataset carregado: {len(seqs)} sequências de tamanho {len(seqs[0])}.")


    # Solicitar número de execuções do BLF-GA
    exec_input = input("\nQuantas execuções do BLF-GA deseja realizar? [1]: ").strip()
    num_execs = int(exec_input) if exec_input.isdigit() and int(exec_input) > 0 else 1

    print(f"\nExecutando baseline (consenso guloso)…")
    t0 = time.time()
    baseline_center = greedy_consensus(seqs, alphabet)
    baseline_val = max_distance(baseline_center, seqs)
    baseline_time = time.time() - t0
    print(f"Baseline: dist={baseline_val}  time={baseline_time:.3f}s")
    logging.debug(f"Centro consenso: {baseline_center}")

    print("\nExecutando BLF-GA…")
    ga_params = BLFGA_DEFAULTS.copy()
    params['blfga_params'] = ga_params

    results = []
    for run_idx in range(num_execs):
        print(f"\n--- Execução {run_idx+1} de {num_execs} ---")
        ga = BLFGA(seqs, alphabet, **ga_params)
        t0 = time.time()
        best, best_val = ga.run()
        ga_time = time.time() - t0
        print(f"\nBLF-GA : dist={best_val}  time={ga_time:.3f}s")
        logging.debug(f"Melhor centro BLF-GA (execução {run_idx+1}): {best}")
        gap = 100 * (baseline_val - best_val) / baseline_val if baseline_val > 0 else 0.0
        results.append({
            'best': best,
            'best_val': best_val,
            'ga_time': ga_time,
            'gap': gap
        })

    # Resumo dos resultados
    print("\nResumo das execuções BLF-GA:")
    for idx, res in enumerate(results, 1):
        print(f"Execução {idx}: dist={res['best_val']}  time={res['ga_time']:.3f}s  melhora={res['gap']:.1f}%")

    # Salvar resultados da melhor execução (menor distância)
    best_result = min(results, key=lambda r: r['best_val'])
    save_results(
        params,
        baseline_center,
        baseline_val,
        baseline_time,
        best_result['best'],
        best_result['best_val'],
        best_result['ga_time'],
        best_result['gap']
    )

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nExecução interrompida pelo usuário. Saindo...")
        sys.exit(0)