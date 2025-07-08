"""
Menu interativo simples para seleção de dataset e algoritmos.

Funções:
    menu(): Exibe menu de datasets e retorna escolha do usuário.
    select_algorithms(): Exibe menu de algoritmos e retorna lista selecionada.
    select_optimization_algorithm(): Seleção de algoritmo para otimização.
    select_sensitivity_algorithm(): Seleção de algoritmo para análise de sensibilidade.
    configure_optimization_params(): Configuração de parâmetros de otimização.
    configure_sensitivity_params(): Configuração de parâmetros de análise de sensibilidade.
"""

import os

from algorithms.base import global_registry
from src.utils.config import safe_input


def menu() -> str:
    """
    Exibe o menu principal para seleção do tipo de dataset usando console simples.

    Returns:
        str: Opção escolhida pelo usuário ('1', '2', '3', '4', '5' ou '6').
    """
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return "1"  # Gerar dataset sintético

    print("\n" + "="*60)
    print("           CLOSEST STRING PROBLEM - CSP-BLFGA")
    print("="*60)
    print("")
    print("� EXECUÇÃO:")
    print("   1) Dataset sintético")
    print("   2) Dataset de arquivo")
    print("   3) Dataset do NCBI")
    print("   4) Execução em lote")
    print("")
    print("🔬 OTIMIZAÇÃO:")
    print("   5) Otimização de hiperparâmetros")
    print("   6) Análise de sensibilidade")
    print("")
    print("💡 SOBRE AS OPÇÕES:")
    print("   • Opções 1-3: Executa algoritmos em datasets individuais")
    print("   • Opção 4: Executa múltiplos algoritmos em configurações batch")
    print("   • Opção 5: Encontra os melhores parâmetros usando Optuna")
    print("   • Opção 6: Analisa sensibilidade dos parâmetros dos algoritmos")
    print("")

    while True:
        c = safe_input("Escolha uma opção [1-6]: ")
        if c in {"1", "2", "3", "4", "5", "6"}:
            return c
        print("❌ Opção inválida. Por favor, escolha uma opção entre 1 e 6.")


def select_algorithms() -> list[str]:
    """
    Exibe menu de seleção de algoritmos disponíveis usando console simples.

    Returns:
        list[str]: Lista com os nomes dos algoritmos selecionados.
    """
    all_algs = list(global_registry.keys())
    # Modo automatizado para testes
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return [all_algs[0]] if all_algs else []

    print("\nAlgoritmos disponíveis:")
    print(" 0) Executar todos")
    for idx, name in enumerate(all_algs, 1):
        print(f" {idx}) {name}")

    selected = []

    raw = safe_input("Escolha (ex.: 1,3 ou 0 para todos) [padrão 1]: ")
    if not raw:
        return [all_algs[0]] if all_algs else []
    if raw == "0":
        return all_algs
    for part in raw.split(","):
        if part.strip().isdigit():
            i = int(part)
            if 1 <= i <= len(all_algs):
                selected.append(all_algs[i - 1])
    return selected


def select_optimization_algorithm() -> str:
    """
    Exibe menu para seleção do algoritmo para otimização usando console simples.

    Returns:
        str: Nome do algoritmo selecionado.
    """
    all_algs = list(global_registry.keys())
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return all_algs[0] if all_algs else ""

    print("\nSelecione o algoritmo para otimização:")
    for idx, name in enumerate(all_algs, 1):
        print(f" {idx}) {name}")

    while True:
        choice = safe_input("Escolha o algoritmo [1]: ").strip()
        if not choice:
            return all_algs[0] if all_algs else ""

        # Tentar por índice numérico
        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(all_algs):
                return all_algs[idx - 1]

        # Tentar por nome do algoritmo
        choice_upper = choice.upper()
        for alg_name in all_algs:
            if alg_name.upper() == choice_upper:
                return alg_name

        # Casos especiais para compatibilidade
        if choice_upper in ["BLF-GA", "BLFGA", "BLF_GA"]:
            if "BLF-GA" in all_algs:
                return "BLF-GA"

        print("Opção inválida. Tente novamente.")


def select_sensitivity_algorithm() -> str:
    """
    Exibe menu para seleção do algoritmo para análise de sensibilidade usando console simples.

    Returns:
        str: Nome do algoritmo selecionado.
    """
    all_algs = list(global_registry.keys())
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return all_algs[0] if all_algs else ""

    print("\nSelecione o algoritmo para análise de sensibilidade:")
    for idx, name in enumerate(all_algs, 1):
        print(f" {idx}) {name}")

    while True:
        choice = safe_input("Escolha o algoritmo [1]: ")
        if not choice:
            return all_algs[0] if all_algs else ""

        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(all_algs):
                return all_algs[idx - 1]

        print("Opção inválida. Tente novamente.")


def configure_optimization_params() -> dict:
    """
    Configura parâmetros para otimização de hiperparâmetros usando console simples.

    Returns:
        dict: Dicionário com configurações da otimização.
    """
    print("\n=== Configuração da Otimização ===")

    # Número de trials
    n_trials_input = safe_input("Número de trials [100]: ")
    n_trials = int(n_trials_input) if n_trials_input.isdigit() else 100

    # Timeout por trial
    timeout_input = safe_input("Timeout por trial em segundos [60]: ")
    timeout = int(timeout_input) if timeout_input.isdigit() else 60

    # Direção da otimização
    print("\nDireção da otimização:")
    print("1) Minimizar")
    print("2) Maximizar")

    direction_input = safe_input("Escolha [1]: ")
    direction = "minimize" if direction_input != "2" else "maximize"

    # Salvar plots
    save_plots_input = safe_input("Salvar gráficos de visualização? (s/N): ")
    save_plots = save_plots_input.lower() in ["s", "sim", "y", "yes"]

    return {
        "n_trials": n_trials,
        "timeout": timeout,
        "direction": direction,
        "save_plots": save_plots,
    }


def configure_sensitivity_params() -> dict:
    """
    Configura parâmetros para análise de sensibilidade usando console simples.

    Returns:
        dict: Dicionário com configurações da análise.
    """
    print("\n=== Configuração da Análise de Sensibilidade ===")

    # Número de amostras
    n_samples_input = safe_input("Número de amostras [1000]: ")
    n_samples = int(n_samples_input) if n_samples_input.isdigit() else 1000

    # Método de análise
    print("\nMétodo de análise:")
    print("1) Sobol")
    print("2) Morris")
    print("3) FAST")

    method_input = safe_input("Escolha [1]: ")
    method_map = {"1": "sobol", "2": "morris", "3": "fast"}
    method = method_map.get(method_input, "sobol")

    # Timeout por amostra
    timeout_input = safe_input("Timeout por amostra em segundos [60]: ")
    timeout = int(timeout_input) if timeout_input.isdigit() else 60

    # Salvar plots
    save_plots_input = safe_input("Salvar gráficos de visualização? (s/N): ")
    save_plots = save_plots_input.lower() in ["s", "sim", "y", "yes"]

    return {
        "n_samples": n_samples,
        "method": method,
        "timeout": timeout,
        "save_plots": save_plots,
    }


def select_dataset_for_optimization() -> tuple:
    """
    Permite selecionar dataset para otimização.

    Returns:
        tuple: (sequences, alphabet, dataset_info)
    """
    print("\n=== Seleção de Dataset para Otimização ===")
    print("1) Dataset sintético (configurável)")
    print("2) Carregar de arquivo")
    print("3) Dataset via NCBI")
    print("4) Dataset sintético (padrão)")

    choice = safe_input("Escolha o tipo de dataset [4]: ")

    if choice == "1":
        # Dataset sintético configurável
        print("\n--- Configuração do Dataset Sintético ---")
        n_input = safe_input("Número de sequências [20]: ")
        n = int(n_input) if n_input.isdigit() else 20

        L_input = safe_input("Tamanho das sequências [100]: ")
        L = int(L_input) if L_input.isdigit() else 100

        print("\nAlfabetos disponíveis:")
        print("1) DNA (ACGT)")
        print("2) RNA (ACGU)")
        print("3) Proteína (20 aminoácidos)")
        print("4) Binário (AB)")
        print("5) Personalizado")

        alpha_choice = safe_input("Escolha o alfabeto [1]: ")
        alphabet_map = {
            "1": "ACGT",
            "2": "ACGU",
            "3": "ACDEFGHIKLMNPQRSTVWY",
            "4": "AB",
            "5": None,
        }
        alphabet = alphabet_map.get(alpha_choice, "ACGT")

        if alphabet is None:
            alphabet = safe_input("Digite o alfabeto personalizado: ").upper()

        noise_input = safe_input("Nível de ruído [0.1]: ")
        noise = float(noise_input) if noise_input.replace(".", "").isdigit() else 0.1

        fully_random_input = safe_input("Completamente aleatório? (s/N): ")
        fully_random = fully_random_input.lower() in ["s", "sim", "y", "yes"]

        from src.datasets.dataset_synthetic import generate_dataset_from_params

        seqs, params = generate_dataset_from_params(
            n=n, L=L, alphabet=alphabet, noise=noise, fully_random=fully_random
        )

        dataset_info = {
            "type": "synthetic_custom",
            "n": n,
            "L": L,
            "alphabet": alphabet,
            "noise": noise,
            "fully_random": fully_random,
            **params,
        }

    elif choice == "2":
        # Carregar de arquivo
        from src.datasets.dataset_file import load_dataset

        seqs, params = load_dataset(silent=False)
        alphabet = "".join(sorted(set("".join(seqs))))
        dataset_info = {"type": "file", **params}

    elif choice == "3":
        # Dataset via NCBI
        from src.datasets.dataset_entrez import fetch_dataset

        seqs, params = fetch_dataset()
        alphabet = "".join(sorted(set("".join(seqs))))
        dataset_info = {"type": "entrez", **params}

    else:
        # Dataset sintético padrão
        from src.datasets.dataset_synthetic import generate_dataset

        seqs, params = generate_dataset(silent=True)
        alphabet = "".join(sorted(set("".join(seqs))))
        dataset_info = {"type": "synthetic_default", **params}

    print(f"\n✅ Dataset carregado: {len(seqs)} sequências de tamanho {len(seqs[0])}")
    print(f"🔤 Alfabeto: {alphabet}")

    return seqs, alphabet, dataset_info


def configure_batch_optimization_params() -> dict:
    """
    Configura parâmetros para otimização em lote.

    Returns:
        dict: Configurações do batch de otimização
    """
    print("\n=== Configuração de Otimização em Lote ===")

    # Número de trials por configuração
    n_trials_input = safe_input("Número de trials por configuração [50]: ")
    n_trials = int(n_trials_input) if n_trials_input.isdigit() else 50

    # Timeout por trial
    timeout_input = safe_input("Timeout por trial em segundos [60]: ")
    timeout = int(timeout_input) if timeout_input.isdigit() else 60

    # Número de configurações de dataset diferentes
    n_configs_input = safe_input("Número de configurações de dataset [3]: ")
    n_configs = int(n_configs_input) if n_configs_input.isdigit() else 3

    # Algoritmos para testar
    print("\nAlgoritmos disponíveis para otimização:")
    from algorithms.base import global_registry

    available_algs = [
        name for name in global_registry.keys() if name not in ["Baseline"]
    ]  # Excluir baseline da otimização

    for idx, alg in enumerate(available_algs, 1):
        print(f" {idx}) {alg}")

    alg_input = safe_input("Escolha algoritmos (ex: 1,2 ou 'todos') [1]: ")

    if alg_input.lower() == "todos":
        selected_algs = available_algs
    elif "," in alg_input:
        indices = [int(x.strip()) for x in alg_input.split(",") if x.strip().isdigit()]
        selected_algs = [
            available_algs[i - 1] for i in indices if 1 <= i <= len(available_algs)
        ]
    elif alg_input.isdigit():
        idx = int(alg_input)
        selected_algs = (
            [available_algs[idx - 1]]
            if 1 <= idx <= len(available_algs)
            else [available_algs[0]]
        )
    else:
        selected_algs = [available_algs[0]]

    # Salvar resultados
    save_results_input = safe_input("Salvar resultados detalhados? (S/n): ")
    save_results = save_results_input.lower() not in ["n", "no", "nao"]

    return {
        "n_trials": n_trials,
        "timeout": timeout,
        "n_configs": n_configs,
        "algorithms": selected_algs,
        "save_results": save_results,
    }


def create_custom_optimization_config() -> dict:
    """
    Cria uma configuração de otimização personalizada interativamente.

    Returns:
        dict: Configuração completa para otimização em lote
    """
    print("\n=== Criador de Configuração de Otimização Personalizada ===")

    # Informações básicas
    print("\n--- Informações Básicas ---")
    nome = safe_input("Nome da configuração: ")
    descricao = safe_input("Descrição: ")

    # Configurações de otimização
    print("\n--- Configurações de Otimização ---")
    n_trials_input = safe_input("Número de trials por algoritmo [30]: ")
    n_trials = int(n_trials_input) if n_trials_input.isdigit() else 30

    timeout_input = safe_input("Timeout por trial em segundos [60]: ")
    timeout_per_trial = int(timeout_input) if timeout_input.isdigit() else 60

    # Seleção de algoritmos
    print("\n--- Seleção de Algoritmos ---")
    from algorithms.base import global_registry

    available_algs = [
        name for name in global_registry.keys() if name not in ["Baseline"]
    ]

    print("Algoritmos disponíveis:")
    for idx, alg in enumerate(available_algs, 1):
        print(f" {idx}) {alg}")

    alg_input = safe_input("Escolha algoritmos (ex: 1,2 ou 'todos') [1]: ")

    if alg_input.lower() == "todos":
        selected_algs = available_algs
    elif "," in alg_input:
        indices = [int(x.strip()) for x in alg_input.split(",") if x.strip().isdigit()]
        selected_algs = [
            available_algs[i - 1] for i in indices if 1 <= i <= len(available_algs)
        ]
    elif alg_input.isdigit():
        idx = int(alg_input)
        selected_algs = (
            [available_algs[idx - 1]]
            if 1 <= idx <= len(available_algs)
            else [available_algs[0]]
        )
    else:
        selected_algs = [available_algs[0]]

    print(f"Algoritmos selecionados: {', '.join(selected_algs)}")

    # Configuração de datasets
    print("\n--- Configuração de Datasets ---")
    datasets = []

    while True:
        print(f"\n--- Dataset {len(datasets) + 1} ---")
        print("Tipos disponíveis:")
        print("1) Sintético")
        print("2) Arquivo")
        print("3) NCBI")

        dataset_type_input = safe_input("Escolha o tipo de dataset [1]: ")

        if dataset_type_input == "2":
            # Dataset de arquivo
            nome_dataset = safe_input("Nome do dataset: ")
            file_path = safe_input("Caminho do arquivo: ")

            datasets.append(
                {
                    "nome": nome_dataset,
                    "tipo": "file",
                    "parametros": {"file_path": file_path},
                }
            )

        elif dataset_type_input == "3":
            # Dataset NCBI
            nome_dataset = safe_input("Nome do dataset: ")
            email = safe_input("Email para NCBI: ")
            db = safe_input("Base de dados [nucleotide]: ") or "nucleotide"
            term = safe_input("Termo de busca: ")
            n_seqs = safe_input("Número de sequências [10]: ")
            n_seqs = int(n_seqs) if n_seqs.isdigit() else 10

            datasets.append(
                {
                    "nome": nome_dataset,
                    "tipo": "entrez",
                    "parametros": {"email": email, "db": db, "term": term, "n": n_seqs},
                }
            )

        else:
            # Dataset sintético (padrão)
            nome_dataset = safe_input("Nome do dataset: ")

            n_input = safe_input("Número de sequências [20]: ")
            n = int(n_input) if n_input.isdigit() else 20

            L_input = safe_input("Tamanho das sequências [100]: ")
            L = int(L_input) if L_input.isdigit() else 100

            print("Alfabetos disponíveis:")
            print("1) DNA (ACGT)")
            print("2) RNA (ACGU)")
            print("3) Proteína (20 aminoácidos)")
            print("4) Binário (AB)")
            print("5) Personalizado")

            alpha_choice = safe_input("Escolha o alfabeto [1]: ")
            alphabet_map = {
                "1": "ACGT",
                "2": "ACGU",
                "3": "ACDEFGHIKLMNPQRSTVWY",
                "4": "AB",
                "5": None,
            }
            alphabet = alphabet_map.get(alpha_choice, "ACGT")

            if alphabet is None:
                alphabet = safe_input("Digite o alfabeto personalizado: ").upper()

            noise_input = safe_input("Nível de ruído [0.1]: ")
            noise = (
                float(noise_input) if noise_input.replace(".", "").isdigit() else 0.1
            )

            fully_random_input = safe_input("Completamente aleatório? (s/N): ")
            fully_random = fully_random_input.lower() in ["s", "sim", "y", "yes"]

            seed_input = safe_input("Semente (opcional): ")
            seed = int(seed_input) if seed_input.isdigit() else None

            datasets.append(
                {
                    "nome": nome_dataset,
                    "tipo": "synthetic",
                    "parametros": {
                        "n": n,
                        "L": L,
                        "alphabet": alphabet,
                        "noise": noise,
                        "fully_random": fully_random,
                        "seed": seed,
                    },
                }
            )

        mais_input = safe_input("Adicionar mais datasets? (s/N): ")
        if mais_input.lower() not in ["s", "sim", "y", "yes"]:
            break

    # Configurações de saída
    print("\n--- Configurações de Saída ---")
    save_plots_input = safe_input("Salvar gráficos? (S/n): ")
    save_plots = save_plots_input.lower() not in ["n", "no", "nao"]

    results_dir = safe_input(
        "Diretório de resultados [outputs/batch_optimization_custom]: "
    )
    if not results_dir:
        results_dir = "outputs/batch_optimization_custom"

    # Montar configuração final
    config = {
        "batch_info": {
            "nome": nome,
            "descricao": descricao,
            "timeout_global": timeout_per_trial
            * n_trials
            * len(selected_algs)
            * len(datasets)
            * 1.5,  # Estimativa
        },
        "optimization_config": {
            "n_trials": n_trials,
            "timeout_per_trial": timeout_per_trial,
            "direction": "minimize",
            "sampler": "TPE",
            "pruner": "Median",
            "save_plots": save_plots,
        },
        "algorithms": selected_algs,
        "datasets": datasets,
        "algorithm_params": {},
        "output": {
            "save_detailed_results": True,
            "save_plots": save_plots,
            "results_dir": results_dir,
            "plot_format": "png",
        },
    }

    return config


def save_optimization_config(config: dict, filename: str | None = None) -> str:
    """
    Salva uma configuração de otimização em arquivo YAML.

    Args:
        config: Configuração a ser salva
        filename: Nome do arquivo (opcional)

    Returns:
        str: Caminho do arquivo salvo
    """
    import os
    from datetime import datetime

    import yaml

    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"otimizacao_personalizada_{timestamp}.yaml"

    config_dir = "batch_configs"
    os.makedirs(config_dir, exist_ok=True)

    filepath = os.path.join(config_dir, filename)

    with open(filepath, "w", encoding="utf-8") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True, indent=2)

    print(f"✅ Configuração salva em: {filepath}")
    return filepath


def interactive_optimization_menu():
    """
    Menu interativo para configuração e execução de otimização.
    """
    print("\n" + "="*60)
    print("           OTIMIZAÇÃO DE HIPERPARÂMETROS")
    print("="*60)
    print("")
    print("🔍 TIPOS DE OTIMIZAÇÃO:")
    print("   1) Otimização simples (dataset customizado)")
    print("   2) Otimização em lote (arquivos YAML)")
    print("")
    print("💡 DICAS:")
    print("   • Opção 1: Escolha dataset interativamente e otimize")
    print("   • Opção 2: Use arquivos YAML para otimização em lote")
    print("   • Resultados salvos automaticamente em outputs/")
    print("")

    choice = safe_input("Escolha uma opção [1]: ")

    if choice == "1":
        # Otimização com dataset customizado
        from src.optimization.optuna_optimizer import (
            run_optimization_with_dataset_selection,
        )

        run_optimization_with_dataset_selection()

    elif choice == "2":
        # Executar otimização em lote YAML
        from src.optimization.batch_optimizer import run_yaml_optimization_batch

        run_yaml_optimization_batch()

    # Sempre retorna ao menu principal após execução
