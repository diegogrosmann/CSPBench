"""
Menu interativo para seleção de dataset e algoritmos.

Funções:
    menu(): Exibe menu de datasets e retorna escolha do usuário.
    select_algorithms(): Exibe menu de algoritmos e retorna lista selecionada.
"""

import os

from algorithms.base import global_registry
from src.ui.cli.console_manager import console
from src.utils.config import safe_input


def menu() -> str:
    """
    Exibe o menu principal para seleção do tipo de dataset.

    Returns:
        str: Opção escolhida pelo usuário ('1', '2', '3', '4' ou '5').
    """
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return "1"  # Gerar dataset sintético
    console.print("\n=== Closest String Problem ===")
    console.print("1) Gerar dataset sintético")
    console.print("2) Carregar dataset de arquivo")
    console.print("3) Baixar dataset via NCBI")
    console.print("4) Execução em lote (batch)")
    console.print("5) Execução em lote com interface curses")
    console.print("6) Otimização de hiperparâmetros")
    console.print("7) Análise de sensibilidade")
    while True:
        c = safe_input("Escolha [1/2/3/4/5/6/7]: ")
        if c in {"1", "2", "3", "4", "5", "6", "7"}:
            return c
        console.print("Inválido.")


def select_algorithms() -> list[str]:
    """
    Exibe menu de seleção de algoritmos disponíveis.

    Returns:
        list[str]: Lista com os nomes dos algoritmos selecionados.
    """
    all_algs = list(global_registry.keys())
    # Modo automatizado para testes
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return [all_algs[0]] if all_algs else []
    console.print("\nAlgoritmos disponíveis:")
    console.print(" 0) Executar todos")
    for idx, name in enumerate(all_algs, 1):
        console.print(f" {idx}) {name}")
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
    Exibe menu para seleção do algoritmo para otimização.

    Returns:
        str: Nome do algoritmo selecionado.
    """
    all_algs = list(global_registry.keys())
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return all_algs[0] if all_algs else ""

    console.print("\nSelecione o algoritmo para otimização:")
    for idx, name in enumerate(all_algs, 1):
        console.print(f" {idx}) {name}")

    while True:
        choice = safe_input("Escolha o algoritmo [1]: ")
        if not choice:
            return all_algs[0] if all_algs else ""

        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(all_algs):
                return all_algs[idx - 1]

        console.print("Opção inválida. Tente novamente.")


def select_sensitivity_algorithm() -> str:
    """
    Exibe menu para seleção do algoritmo para análise de sensibilidade.

    Returns:
        str: Nome do algoritmo selecionado.
    """
    all_algs = list(global_registry.keys())
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return all_algs[0] if all_algs else ""

    console.print("\nSelecione o algoritmo para análise de sensibilidade:")
    for idx, name in enumerate(all_algs, 1):
        console.print(f" {idx}) {name}")

    while True:
        choice = safe_input("Escolha o algoritmo [1]: ")
        if not choice:
            return all_algs[0] if all_algs else ""

        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(all_algs):
                return all_algs[idx - 1]

        console.print("Opção inválida. Tente novamente.")


def configure_optimization_params() -> dict:
    """
    Configura parâmetros para otimização.

    Returns:
        dict: Parâmetros de configuração para otimização.
    """
    console.print("\n=== Configuração da Otimização ===")

    # Número de trials
    n_trials = safe_input("Número de trials [100]: ")
    n_trials = int(n_trials) if n_trials.isdigit() else 100

    # Timeout por trial
    timeout = safe_input("Timeout por trial em segundos [60]: ")
    timeout = int(timeout) if timeout.isdigit() else 60

    # Direção da otimização
    console.print("\nDireção da otimização:")
    console.print("1) Minimizar")
    console.print("2) Maximizar")
    direction_choice = safe_input("Escolha [1]: ")
    direction = "minimize" if direction_choice != "2" else "maximize"

    # Salvar visualizações
    save_plots = safe_input("Salvar gráficos de visualização? [s/N]: ").lower()
    save_plots = save_plots in ["s", "sim", "y", "yes"]

    return {
        "n_trials": n_trials,
        "timeout": timeout,
        "direction": direction,
        "save_plots": save_plots,
    }


def configure_sensitivity_params() -> dict:
    """
    Configura parâmetros para análise de sensibilidade.

    Returns:
        dict: Parâmetros de configuração para análise de sensibilidade.
    """
    console.print("\n=== Configuração da Análise de Sensibilidade ===")

    # Número de amostras
    n_samples = safe_input("Número de amostras [1000]: ")
    n_samples = int(n_samples) if n_samples.isdigit() else 1000

    # Método de análise
    console.print("\nMétodo de análise:")
    console.print("1) Sobol")
    console.print("2) Morris")
    console.print("3) FAST")
    method_choice = safe_input("Escolha [1]: ")
    method_map = {"1": "sobol", "2": "morris", "3": "fast"}
    method = method_map.get(method_choice, "sobol")

    # Salvar visualizações
    save_plots = safe_input("Salvar gráficos de visualização? [s/N]: ").lower()
    save_plots = save_plots in ["s", "sim", "y", "yes"]

    return {
        "n_samples": n_samples,
        "method": method,
        "save_plots": save_plots,
    }
