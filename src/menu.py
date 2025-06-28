from algorithms.base import global_registry
from src.console_manager import console

def menu() -> str:
    console.print("\n=== Closest String Problem ===")
    console.print("1) Gerar dataset sintético")
    console.print("2) Carregar dataset de arquivo")
    console.print("3) Baixar dataset via NCBI")
    while True:
        try:
            c = input("Escolha [1/2/3]: ").strip()
            if c in {'1', '2', '3'}:
                return c
            console.print("Inválido.")
        except KeyboardInterrupt:
            console.print("\nOperação cancelada pelo usuário.")
            import sys
            sys.exit(0)

def select_algorithms() -> list[str]:
    all_algs = list(global_registry.keys())
    algs = [name for name in all_algs if name != "Baseline"]
    console.print("\nAlgoritmos disponíveis:")
    console.print(" 0) Executar todos")
    for idx, name in enumerate(algs, 1):
        console.print(f" {idx}) {name}")
    selected = []
    try:
        raw = input("Escolha (ex.: 1,3 ou 0 para todos) [padrão 1]: ").strip()
    except KeyboardInterrupt:
        console.print("\nOperação cancelada pelo usuário.")
        import sys
        sys.exit(0)
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
