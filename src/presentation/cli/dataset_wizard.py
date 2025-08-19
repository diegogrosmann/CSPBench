"""Dataset Wizard (CLI)

Substitui antigo módulo TUI descontinuado. Fornece uma interface
simples baseada em ``input()`` para coletar parâmetros de geração
de datasets sintéticos ou reais. A API mantém os mesmos métodos
utilizados pelo ``DatasetGenerationOrchestrator`` para minimizar
impacto em código existente:

    - show_main_menu() -> str
    - collect_synthetic_params() -> dict
    - collect_real_params() -> dict
    - generate_default_filename(kind, params) -> str
    - get_output_filename(default) -> str

Uso: instanciado internamente pelo orchestrator (datasetsave).
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Dict, Optional


def _ask(prompt: str, default: str | None = None) -> str:
    # Mantém comportamento simples: se Enter e default não é None, retorna default; caso contrário, string vazia
    suffix = f" [{default}]" if default is not None else " [None]"
    value = input(f"{prompt}{suffix}: ").strip()
    if not value and default is not None:
        return default
    return value


def _ask_int(prompt: str, default: int, min_value: int | None = 1) -> int:
    while True:
        raw = _ask(prompt, str(default))
        try:
            value = int(raw)
            if min_value is not None and value < min_value:
                print(f"⚠️  Valor deve ser >= {min_value}.")
                continue
            return value
        except ValueError:
            print("⚠️  Informe um número inteiro válido.")


def _ask_float(
    prompt: str,
    default: float,
    min_value: float | None = 0.0,
    max_value: float | None = 1.0,
) -> float:
    while True:
        raw = _ask(prompt, str(default))
        try:
            value = float(raw)
            if min_value is not None and value < min_value:
                print(f"⚠️  Valor deve ser >= {min_value}.")
                continue
            if max_value is not None and value > max_value:
                print(f"⚠️  Valor deve ser <= {max_value}.")
                continue
            return value
        except ValueError:
            print("⚠️  Informe um número numérico válido.")


def _ask_optional_int(prompt: str, default: Optional[int]) -> Optional[int]:
    while True:
        suffix = f" [{default}]" if default is not None else " [None]"
        raw = input(f"{prompt}{suffix}: ").strip()
        if raw == "":
            return default  # pode ser None
        try:
            return int(raw)
        except ValueError:
            print("⚠️  Informe um número inteiro válido ou deixe vazio para None.")


@dataclass
class DatasetWizard:
    """CLI dataset wizard substituindo implementação TUI antiga."""

    def show_main_menu(self) -> str:
        print("\n=== Dataset Wizard ===")
        print("1) Gerar sintético")
        print("2) Obter real (NCBI)")
        print("0) Sair")
        choice = _ask("Escolha", "1")
        if choice == "1":
            return "synthetic"
        if choice == "2":
            return "real"
        return "exit"

    # -------- Synthetic --------
    def collect_synthetic_params(self) -> Dict[str, Any]:  # noqa: D401
        # Importar valores padrão
        from src.application.services.dataset_generator import SYNTHETIC_DEFAULTS

        print("\n--- Parâmetros Sintético ---")
        methods = ["random", "noise", "clustered", "mutations"]
        print("Métodos:")
        for i, m in enumerate(methods, 1):
            print(f"  {i}) {m}")
        method_choice = _ask("Método", "1")
        try:
            method = methods[int(method_choice) - 1]
        except Exception:  # noqa: BLE001
            method = "random"

        n = _ask_int("Quantidade de sequências (n)", SYNTHETIC_DEFAULTS["n"], 1)
        length = _ask_int(
            "Tamanho de cada sequência (length)", SYNTHETIC_DEFAULTS["length"], 1
        )
        alphabet = _ask("Alfabeto", SYNTHETIC_DEFAULTS["alphabet"]).upper()
        seed_raw = _ask("Seed (opcional)", "")
        seed = int(seed_raw) if seed_raw else None

        params: Dict[str, Any] = {
            "method": method,
            "n": n,
            "length": length,
            "alphabet": alphabet,
        }
        if seed is not None:
            params["seed"] = seed

        if method == "noise":
            params["noise_rate"] = _ask_float(
                "Taxa de ruído (0-1)", SYNTHETIC_DEFAULTS["noise_rate"], 0.0, 1.0
            )
            center = _ask("Sequência central (vazio=auto)", "")
            if center:
                params["center_sequence"] = center.upper()
        elif method == "clustered":
            params["num_clusters"] = _ask_int(
                "Número de clusters", SYNTHETIC_DEFAULTS["num_clusters"], 1
            )
            params["cluster_noise"] = _ask_float(
                "Ruído intra cluster (0-1)", SYNTHETIC_DEFAULTS["cluster_distance"]
            )
        elif method == "mutations":
            params["mutation_rate"] = _ask_float(
                "Taxa de mutação (0-1)", SYNTHETIC_DEFAULTS["mutation_rate"]
            )
            base = _ask("Sequência base (vazio=auto)", "")
            if base:
                params["base_sequence"] = base.upper()
        return params

    # -------- Real (NCBI) --------
    def collect_real_params(self) -> Dict[str, Any]:
        from src.infrastructure.external.dataset_entrez import ENTREZ_DEFAULTS

        print("\n--- Dataset Real (NCBI) ---")
        print(
            "ℹ️ Para mais informações sobre a sintaxe de consultas Entrez, consulte: https://www.ncbi.nlm.nih.gov/books/NBK25499/"
        )
        # Apenas NCBI é suportado agora
        db = _ask("Banco (db)", ENTREZ_DEFAULTS["db"]) or ENTREZ_DEFAULTS["db"]
        query = _ask("Consulta (query)", ENTREZ_DEFAULTS["query"]) or "*"
        max_seq = _ask_optional_int(
            "Máx. sequências", ENTREZ_DEFAULTS["max_sequences"]
        )  # None permitido
        min_len = _ask_optional_int(
            "Comprimento mínimo", ENTREZ_DEFAULTS["min_length"]
        )  # None permitido
        max_len = _ask_optional_int(
            "Comprimento máximo", ENTREZ_DEFAULTS["max_length"]
        )  # None permitido

        # Validar relação entre min/max quando ambos informados
        while min_len is not None and max_len is not None and max_len < min_len:
            print("⚠️  'Comprimento máximo' deve ser >= 'Comprimento mínimo'.")
            max_len = _ask_optional_int("Comprimento máximo", None)
        return {
            "query": query,
            "db": db,
            "max_sequences": max_seq,
            "min_length": min_len,
            "max_length": max_len,
        }

    # -------- Util --------
    def generate_default_filename(self, kind: str, params: Dict[str, Any]) -> str:
        if kind == "synthetic":
            method = params.get("method", "random")
            n = params.get("n", 0)
            length = params.get("length", 0)
            return f"synthetic_{method}_n{n}_L{length}.fasta"
        if kind == "real":
            # Apenas NCBI suportado: usar query para compor nome
            query = params.get("query", "q")
            safe = re.sub(r"[^A-Za-z0-9_]+", "_", query)[:30].strip("_") or "query"
            return f"real_{safe}.fasta"
        return "dataset.fasta"

    def get_output_filename(self, default: str) -> str:
        print(f"\nNome de arquivo sugerido: {default}")
        name = _ask("Arquivo de saída", default)
        if not name.lower().endswith(".fasta"):
            name += ".fasta"
        return name
