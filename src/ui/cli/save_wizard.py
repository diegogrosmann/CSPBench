"""
Wizard para salvar datasets.

Este módulo fornece uma interface para permitir que o usuário
escolha se deseja salvar datasets gerados.
"""

import logging
from pathlib import Path

from src.utils.config import safe_input

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


def ask_save_dataset(sequences: list[str], dataset_type: str, params: dict) -> bool:
    """
    Pergunta se o usuário deseja salvar o dataset e salva se confirmado.

    Args:
        sequences: Lista de sequências
        dataset_type: Tipo do dataset ('synthetic', 'entrez', etc.)
        params: Parâmetros usados para gerar o dataset

    Returns:
        True se o dataset foi salvo, False caso contrário
    """
    save_choice = safe_input(
        "\nDeseja salvar este dataset para uso futuro? [s/N]: "
    ).lower()

    if save_choice == "s":
        # Gerar nome do arquivo baseado no tipo e parâmetros
        filename, description = _generate_filename_and_description(dataset_type, params)

        try:
            saved_path = save_dataset_fasta(sequences, filename, description)
            print(f"Dataset salvo em: {saved_path}")
            return True
        except Exception as e:
            print(f"Erro ao salvar dataset: {e}")
            logger.error("Erro ao salvar dataset: %s", e)
            return False

    return False


def _generate_filename_and_description(
    dataset_type: str, params: dict
) -> tuple[str, str]:
    """
    Gera nome do arquivo e descrição baseados no tipo e parâmetros.

    Args:
        dataset_type: Tipo do dataset
        params: Parâmetros do dataset

    Returns:
        Tupla (filename, description)
    """
    if dataset_type == "synthetic":
        filename = (
            f"synthetic_n{params['n']}_L{params['L']}_noise{params['noise']}.fasta"
        )
        description = f"Synthetic dataset: n={params['n']}, L={params['L']}, alphabet={params['alphabet']}, noise={params['noise']}"
    elif dataset_type == "entrez":
        # Limpar caracteres especiais do termo para nome do arquivo
        clean_term = "".join(c for c in params["term"] if c.isalnum() or c in "_ -")[
            :50
        ]
        filename = f"entrez_{params['db']}_{clean_term}_n{params['n']}.fasta"
        description = (
            f"NCBI dataset: db={params['db']}, term={params['term']}, n={params['n']}"
        )
    else:
        filename = f"dataset_{dataset_type}.fasta"
        description = f"Dataset from {dataset_type}"

    return filename, description
