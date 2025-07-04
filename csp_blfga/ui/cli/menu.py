"""
Menu interativo para seleção de dataset e algoritmos.

Funções:
    menu(): Exibe menu de datasets e retorna escolha do usuário.
    select_algorithms(): Exibe menu de algoritmos e retorna lista selecionada.
"""

import os

from algorithms.base import global_registry
from src.console_manager import console
from utils.config import safe_input


def menu() -> str:
    """
    Exibe o menu principal para seleção do tipo de dataset.

    Returns:
        str: Opção escolhida pelo usuário ('1', '2', '3' ou '4').
    """
    if os.environ.get("CSP_AUTOMATED_TEST") == "1":
        return "1"  # Gerar dataset sintético
    console.print("\n=== Closest String Problem ===")
    console.print("1) Gerar dataset sintético")
    console.print("2) Carregar dataset de arquivo")
    console.print("3) Baixar dataset via NCBI")
    console.print("4) Execução em lote (batch)")
    while True:
        c = safe_input("Escolha [1/2/3/4]: ")
        if c in {"1", "2", "3", "4"}:
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
