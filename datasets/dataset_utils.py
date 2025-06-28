"""
dataset_utils.py
================

Utilitários para salvar e gerenciar datasets.

Funções:
    ensure_datasets_folder(): Garante existência da pasta de datasets.
    save_dataset_fasta(sequences, filename, description): Salva dataset em formato FASTA.
    ask_save_dataset(sequences, dataset_type, params): Pergunta ao usuário e salva dataset se desejado.
"""

import os
from pathlib import Path
from typing import List
import logging
from utils.config import safe_input

logger = logging.getLogger(__name__)

def ensure_datasets_folder():
    """Garante que a pasta saved_datasets/ existe."""
    datasets_path = Path("saved_datasets")
    datasets_path.mkdir(exist_ok=True)
    return datasets_path

def save_dataset_fasta(sequences: List[str], filename: str, description: str = ""):
    """Salva um dataset no formato FASTA na pasta saved_datasets/."""
    datasets_path = ensure_datasets_folder()
    filepath = datasets_path / filename
    
    with filepath.open('w', encoding='utf-8') as f:
        for i, seq in enumerate(sequences):
            header = f">seq_{i+1}"
            if description:
                header += f" {description}"
            f.write(f"{header}\n{seq}\n")
    
    logger.debug(f"Dataset salvo em {filepath}")
    return filepath

def ask_save_dataset(sequences: List[str], dataset_type: str, params: dict) -> bool:
    """Pergunta se o usuário deseja salvar o dataset e salva se confirmado."""
    save_choice = safe_input(f"\nDeseja salvar este dataset para uso futuro? [s/N]: ").lower()
    
    if save_choice == 's':
        # Gerar nome do arquivo baseado no tipo e parâmetros
        if dataset_type == "synthetic":
            filename = f"synthetic_n{params['n']}_L{params['L']}_noise{params['noise']}.fasta"
            description = f"Synthetic dataset: n={params['n']}, L={params['L']}, alphabet={params['alphabet']}, noise={params['noise']}"
        elif dataset_type == "entrez":
            # Limpar caracteres especiais do termo para nome do arquivo
            clean_term = "".join(c for c in params['term'] if c.isalnum() or c in "_ -")[:50]
            filename = f"entrez_{params['db']}_{clean_term}_n{params['n']}.fasta"
            description = f"NCBI dataset: db={params['db']}, term={params['term']}, n={params['n']}"
        else:
            filename = f"dataset_{dataset_type}.fasta"
            description = f"Dataset from {dataset_type}"
        
        try:
            saved_path = save_dataset_fasta(sequences, filename, description)
            print(f"Dataset salvo em: {saved_path}")
            return True
        except Exception as e:
            print(f"Erro ao salvar dataset: {e}")
            logger.error(f"Erro ao salvar dataset: {e}")
            return False
    
    return False
