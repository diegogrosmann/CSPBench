"""
dataset_utils.py
================

Utilitários para salvar e gerenciar datasets.

Funções:
    ensure_datasets_folder(): Garante existência da pasta de datasets.
    save_dataset_fasta(sequences, filename, description): Salva dataset em formato FASTA.
    ask_save_dataset(): Redirecionamento para o wizard da UI (deprecated).
"""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def ensure_datasets_folder():
    """Garante que a pasta saved_datasets/ existe."""
    datasets_path = Path("saved_datasets")
    datasets_path.mkdir(exist_ok=True)
    return datasets_path


def save_dataset_fasta(sequences: list[str], filename: str, description: str = ""):
    """Salva um dataset no formato FASTA na pasta saved_datasets/."""
    datasets_path = ensure_datasets_folder()
    filepath = datasets_path / filename

    with filepath.open("w", encoding="utf-8") as f:
        for i, seq in enumerate(sequences):
            header = f">seq_{i+1}"
            if description:
                header += f" {description}"
            f.write(f"{header}\n{seq}\n")

    logger.debug("Dataset salvo em %s", filepath)
    return filepath
