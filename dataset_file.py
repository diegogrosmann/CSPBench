# dataset_file.py
"""
dataset_file.py
===============

Lê um arquivo de texto ou FASTA e extrai as sequências (uma por linha
ou FASTA padrão). A interação é feita no main.
"""

from typing import List, Tuple, Dict, Any
from pathlib import Path
from Bio import SeqIO
import logging
from config import FILE_DEFAULTS
from dataset_utils import ensure_datasets_folder

logger = logging.getLogger(__name__)

def load_dataset() -> Tuple[List[str], Dict[str, Any]]:
    defaults = FILE_DEFAULTS
    
    # Mostrar arquivos disponíveis na pasta datasets
    datasets_path = ensure_datasets_folder()
    available_files = list(datasets_path.glob("*.fasta")) + list(datasets_path.glob("*.fa")) + list(datasets_path.glob("*.txt"))
    
    if available_files:
        print(f"\nArquivos disponíveis na pasta datasets/:")
        for i, file in enumerate(available_files, 1):
            print(f"  {i}. {file.name}")
        
        print(f"\nOpções:")
        print(f"  - Digite o número do arquivo para selecioná-lo")
        print(f"  - Digite o caminho completo do arquivo")
        print(f"  - Pressione Enter para usar o padrão [{defaults['filepath']}]")
        
        path_input = input(f"\nSua escolha: ").strip()
        
        # Verificar se é um número (seleção da lista)
        if path_input.isdigit():
            file_index = int(path_input) - 1
            if 0 <= file_index < len(available_files):
                path_str = str(available_files[file_index])
                print(f"Arquivo selecionado: {available_files[file_index].name}")
            else:
                print(f"Número inválido. Usando padrão: {defaults['filepath']}")
                path_str = defaults['filepath']
        elif path_input:
            # Caminho fornecido pelo usuário
            path_str = path_input
        else:
            # Usar padrão
            path_str = defaults['filepath']
    else:
        print(f"\nNenhum arquivo encontrado na pasta datasets/")
        path_input = input(f"Caminho do arquivo (.txt ou .fasta) [{defaults['filepath']}]: ").strip()
        path_str = path_input if path_input else defaults['filepath']
    
    path = Path(path_str)
    logger.debug(f"Tentando carregar dataset de '{path}'")
    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")

    if path.suffix.lower() in {".fa", ".fasta"}:
        logger.debug("Detectado FASTA; usando SeqIO.parse")
        records = list(SeqIO.parse(path, "fasta"))
        seqs = [str(rec.seq).upper() for rec in records]
    else:
        logger.debug("Detectado texto simples; lendo linha a linha")
        seqs = [line.strip().upper() for line in path.open() if line.strip()]

    print(f"Dataset carregado: n={len(seqs)}, L={len(seqs[0])}")
    params = {'filepath': path_str}
    return seqs, params
