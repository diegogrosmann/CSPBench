from algorithms.base import global_registry
from src.console_manager import console
from utils.config import safe_input

"""
Menu interativo para seleção de dataset e algoritmos.

Funções:
    menu(): Exibe menu de datasets e retorna escolha.
    select_algorithms(): Exibe menu de algoritmos e retorna lista selecionada.
"""

def menu() -> str:
    console.print("\n=== Closest String Problem ===")
    console.print("1) Gerar dataset sintético")
    console.print("2) Carregar dataset de arquivo")
    console.print("3) Baixar dataset via NCBI")
    while True:
        c = safe_input("Escolha [1/2/3]: ")
        if c in {'1', '2', '3'}:
            return c
        console.print("Inválido.")

def select_algorithms() -> list[str]:
    all_algs = list(global_registry.keys())
    algs = [name for name in all_algs if name != "Baseline"]
    console.print("\nAlgoritmos disponíveis:")
    console.print(" 0) Executar todos")
    for idx, name in enumerate(algs, 1):
        console.print(f" {idx}) {name}")
    selected = []
    
    raw = safe_input("Escolha (ex.: 1,3 ou 0 para todos) [padrão 1]: ")
    if not raw:
        return [algs[0]] if algs else []
    if raw == '0':
        return algs
    for part in raw.split(','):
        if part.strip().isdigit():
            i = int(part)
            if 1 <= i <= len(algs):
                selected.append(algs[i-1])
    return selected
    if raw == '0':
        return algs
    for part in raw.split(','):
        if part.strip().isdigit():
            i = int(part)
            if 1 <= i <= len(algs):
                selected.append(algs[i-1])
    return selected
