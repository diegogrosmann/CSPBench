"""Dataset Wizard (CLI)

This module replaces the discontinued TUI module and provides a simple
interface based on ``input()`` for collecting parameters for synthetic
or real dataset generation. The API maintains the same methods used by
the ``DatasetGenerationOrchestrator`` to minimize impact on existing code:

    - show_main_menu() -> str
    - collect_synthetic_params() -> dict
    - collect_real_params() -> dict
    - generate_default_filename(kind, params) -> str
    - get_output_filename(default) -> str

Usage: instantiated internally by the orchestrator (datasetsave).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, Optional


def _ask(prompt: str, default: str | None = None) -> str:
    """Ask user for input with optional default value.

    Args:
        prompt: Question to display to the user
        default: Default value to use if user presses Enter without input

    Returns:
        User input string, or default value if provided and user input is empty
    """
    # Maintains simple behavior: if Enter and default is not None, return default; otherwise, empty string
    suffix = f" [{default}]" if default is not None else " [None]"
    value = input(f"{prompt}{suffix}: ").strip()
    if not value and default is not None:
        return default
    return value


def _ask_int(prompt: str, default: int, min_value: int | None = 1) -> int:
    """Ask user for integer input with validation.

    Args:
        prompt: Question to display to the user
        default: Default value to use if user input is empty
        min_value: Minimum allowed value for validation

    Returns:
        Valid integer value from user input or default
    """
    while True:
        raw = _ask(prompt, str(default))
        try:
            value = int(raw)
            if min_value is not None and value < min_value:
                print(f"⚠️  Value must be >= {min_value}.")
                continue
            return value
        except ValueError:
            print("⚠️  Please enter a valid integer number.")


def _ask_float(
    prompt: str,
    default: float,
    min_value: float | None = 0.0,
    max_value: float | None = 1.0,
) -> float:
    """Ask user for float input with validation.

    Args:
        prompt: Question to display to the user
        default: Default value to use if user input is empty
        min_value: Minimum allowed value for validation
        max_value: Maximum allowed value for validation

    Returns:
        Valid float value from user input or default
    """
    while True:
        raw = _ask(prompt, str(default))
        try:
            value = float(raw)
            if min_value is not None and value < min_value:
                print(f"⚠️  Value must be >= {min_value}.")
                continue
            if max_value is not None and value > max_value:
                print(f"⚠️  Value must be <= {max_value}.")
                continue
            return value
        except ValueError:
            print("⚠️  Please enter a valid numeric value.")


def _ask_optional_int(prompt: str, default: Optional[int]) -> Optional[int]:
    """Ask user for optional integer input.

    Args:
        prompt: Question to display to the user
        default: Default value to use if user input is empty (can be None)

    Returns:
        Valid integer value from user input, default value, or None
    """
    while True:
        suffix = f" [{default}]" if default is not None else " [None]"
        raw = input(f"{prompt}{suffix}: ").strip()
        if raw == "":
            return default  # can be None
        try:
            return int(raw)
        except ValueError:
            print("⚠️  Please enter a valid integer number or leave empty for None.")


@dataclass
class DatasetWizard:
    """CLI dataset wizard replacing old TUI implementation.

    This class provides a command-line interface for collecting dataset generation
    parameters from users. It supports both synthetic dataset generation with
    various methods (random, noise, clustered, mutations) and real dataset
    retrieval from NCBI databases.

    The wizard guides users through parameter collection with input validation
    and provides reasonable defaults for all parameters.
    """

    def show_main_menu(self) -> str:
        """Display main menu and get user choice for dataset type.

        Returns:
            User choice: "synthetic" for synthetic generation,
            "real" for NCBI retrieval, or "exit" to quit
        """
        print("\n=== Dataset Wizard ===")
        print("1) Generate synthetic")
        print("2) Get real (NCBI)")
        print("0) Exit")
        choice = _ask("Choose", "1")
        if choice == "1":
            return "synthetic"
        if choice == "2":
            return "real"
        return "exit"

    # -------- Synthetic --------
    def collect_synthetic_params(self) -> Dict[str, Any]:
        """Collect parameters for synthetic dataset generation.

        Guides the user through collecting all necessary parameters for synthetic
        dataset generation, including method selection, sequence parameters,
        and method-specific options.

        Returns:
            Dictionary containing all synthetic generation parameters
        """
        # Import default values
        from src.application.services.dataset_generator import SYNTHETIC_DEFAULTS

        print("\n--- Synthetic Parameters ---")
        methods = ["random", "noise", "clustered", "mutations"]
        print("Methods:")
        for i, m in enumerate(methods, 1):
            print(f"  {i}) {m}")
        method_choice = _ask("Method", "1")
        try:
            method = methods[int(method_choice) - 1]
        except Exception:  # noqa: BLE001
            method = "random"

        n = _ask_int("Number of sequences (n)", SYNTHETIC_DEFAULTS["n"], 1)
        length = _ask_int(
            "Length of each sequence (length)", SYNTHETIC_DEFAULTS["length"], 1
        )
        alphabet = _ask("Alphabet", SYNTHETIC_DEFAULTS["alphabet"]).upper()
        seed_raw = _ask("Seed (optional)", "")
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
                "Noise rate (0-1)", SYNTHETIC_DEFAULTS["noise_rate"], 0.0, 1.0
            )
            center = _ask("Central sequence (empty=auto)", "")
            if center:
                params["center_sequence"] = center.upper()
        elif method == "clustered":
            params["num_clusters"] = _ask_int(
                "Number of clusters", SYNTHETIC_DEFAULTS["num_clusters"], 1
            )
            params["cluster_noise"] = _ask_float(
                "Intra-cluster noise (0-1)", SYNTHETIC_DEFAULTS["cluster_distance"]
            )
        elif method == "mutations":
            params["mutation_rate"] = _ask_float(
                "Mutation rate (0-1)", SYNTHETIC_DEFAULTS["mutation_rate"]
            )
            base = _ask("Base sequence (empty=auto)", "")
            if base:
                params["base_sequence"] = base.upper()
        return params

    # -------- Real (NCBI) --------
    def collect_real_params(self) -> Dict[str, Any]:
        """Collect parameters for real dataset retrieval from NCBI.

        Guides the user through collecting parameters for retrieving real
        biological sequences from NCBI databases using Entrez queries.

        Returns:
            Dictionary containing all NCBI retrieval parameters
        """
        from src.infrastructure.external.dataset_entrez import ENTREZ_DEFAULTS

        print("\n--- Real Dataset (NCBI) ---")
        print(
            "ℹ️ For more information about Entrez query syntax, see: https://www.ncbi.nlm.nih.gov/books/NBK25499/"
        )
        # Only NCBI is supported now
        db = _ask("Database (db)", ENTREZ_DEFAULTS["db"]) or ENTREZ_DEFAULTS["db"]
        query = _ask("Query (query)", ENTREZ_DEFAULTS["query"]) or "*"
        max_seq = _ask_optional_int(
            "Max sequences", ENTREZ_DEFAULTS["max_sequences"]
        )  # None allowed
        min_len = _ask_optional_int(
            "Minimum length", ENTREZ_DEFAULTS["min_length"]
        )  # None allowed
        max_len = _ask_optional_int(
            "Maximum length", ENTREZ_DEFAULTS["max_length"]
        )  # None allowed

        # Validate relationship between min/max when both are provided
        while min_len is not None and max_len is not None and max_len < min_len:
            print("⚠️  'Maximum length' must be >= 'Minimum length'.")
            max_len = _ask_optional_int("Maximum length", None)
        return {
            "query": query,
            "db": db,
            "max_sequences": max_seq,
            "min_length": min_len,
            "max_length": max_len,
        }

    # -------- Util --------
    def generate_default_filename(self, kind: str, params: Dict[str, Any]) -> str:
        """Generate default filename based on dataset type and parameters.

        Creates a descriptive filename based on the dataset type and key parameters
        to help users identify the generated dataset files.

        Args:
            kind: Dataset type ("synthetic" or "real")
            params: Parameter dictionary containing generation settings

        Returns:
            Suggested filename with .fasta extension
        """
        if kind == "synthetic":
            method = params.get("method", "random")
            n = params.get("n", 0)
            length = params.get("length", 0)
            return f"synthetic_{method}_n{n}_L{length}.fasta"
        if kind == "real":
            # Only NCBI supported: use query to compose name
            query = params.get("query", "q")
            safe = re.sub(r"[^A-Za-z0-9_]+", "_", query)[:30].strip("_") or "query"
            return f"real_{safe}.fasta"
        return "dataset.fasta"

    def get_output_filename(self, default: str) -> str:
        """Get output filename from user with default suggestion.

        Prompts the user for an output filename and ensures it has the correct
        .fasta extension.

        Args:
            default: Default filename to suggest to the user

        Returns:
            User-specified filename with .fasta extension
        """
        print(f"\nSuggested filename: {default}")
        name = _ask("Output file", default)
        if not name.lower().endswith(".fasta"):
            name += ".fasta"
        return name
