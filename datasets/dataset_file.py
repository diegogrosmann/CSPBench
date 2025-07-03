"""
Leitura de arquivos de texto ou FASTA para extração de sequências.

Funções:
    load_dataset(): Lê um arquivo e retorna as sequências e parâmetros.
"""

# dataset_file.py
"""
dataset_file.py
===============

Lê um arquivo de texto ou FASTA e extrai as sequências (uma por linha
ou FASTA padrão). A interação é feita no main.
"""

from typing import List, Tuple, Dict, Any
from pathlib import Path
import logging
from utils.config import FILE_DEFAULTS, safe_input
from src.console_manager import console

from Bio import SeqIO

from datasets.dataset_utils import ensure_datasets_folder

logger = logging.getLogger(__name__)

def load_dataset(silent: bool = False) -> Tuple[List[str], Dict[str, Any]]:
    defaults = FILE_DEFAULTS
    datasets_path = ensure_datasets_folder()
    available_files = list(datasets_path.glob("*.fasta")) + list(datasets_path.glob("*.fa")) + list(datasets_path.glob("*.txt"))
    if silent:
        path_str = defaults['filepath']
    else:
        if available_files:
            print(f"\nArquivos disponíveis na pasta saved_datasets/:")
            for i, file in enumerate(available_files, 1):
                print(f"  {i}. {file.name}")
            print(f"\nOpções:")
            print(f"  - Digite o número do arquivo para selecioná-lo")
            print(f"  - Digite o caminho completo do arquivo")
            print(f"  - Pressione Enter para usar o padrão [{defaults['filepath']}]")
            path_input = input(f"\nSua escolha: ").strip()
            if path_input.isdigit():
                file_index = int(path_input) - 1
                if 0 <= file_index < len(available_files):
                    path_str = str(available_files[file_index])
                    print(f"Arquivo selecionado: {available_files[file_index].name}")
                else:
                    print(f"Número inválido. Usando padrão: {defaults['filepath']}")
                    path_str = defaults['filepath']
            elif path_input:
                path_str = path_input
            else:
                path_str = defaults['filepath']
        else:
            print(f"\nNenhum arquivo encontrado na pasta saved_datasets/")
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
    return seqs, {'filepath': path_str}

def _load_fasta(file_path: Path) -> List[str]:
    """Carrega sequências de arquivo FASTA."""
    sequences = []
    current_seq = ""
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Nova sequência
                    if current_seq:
                        clean_seq = ''.join(c.upper() for c in current_seq if c.isalpha())
                        if clean_seq:
                            sequences.append(clean_seq)
                    current_seq = ""
                else:
                    current_seq += line
            
            # Adicionar última sequência
            if current_seq:
                clean_seq = ''.join(c.upper() for c in current_seq if c.isalpha())
                if clean_seq:
                    sequences.append(clean_seq)
                    
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo FASTA {file_path}: {e}")
    
    return sequences

def _load_text(file_path: Path) -> List[str]:
    """Carrega sequências de arquivo texto simples (uma por linha)."""
    sequences = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):  # Ignora linhas vazias e comentários
                    clean_seq = ''.join(c.upper() for c in line if c.isalpha())
                    if clean_seq:
                        sequences.append(clean_seq)
                    elif line:  # Se tinha conteúdo mas não sobrou nada após limpeza
                        logger.warning(f"Linha {line_num} ignorada (sem caracteres válidos): {line[:50]}")
                        
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo texto {file_path}: {e}")
    
    return sequences

def load_dataset_with_params(params: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
    """
    Carrega dataset de arquivo com parâmetros específicos.
    
    Args:
        params: Dicionário com parâmetros (filepath)
        
    Returns:
        Tupla (sequências, parâmetros_usados)
    """
    # Merge com defaults
    merged_params = {**FILE_DEFAULTS}
    merged_params.update(params)
    
    filepath = merged_params['filepath']
    
    console.print(f"Carregando dataset de arquivo: {filepath}")
    
    file_path = Path(filepath)
    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {filepath}")
    
    sequences = []
    
    if file_path.suffix.lower() in ['.fasta', '.fa', '.fas']:
        sequences = _load_fasta(file_path)
    elif file_path.suffix.lower() == '.txt':
        sequences = _load_text(file_path)
    else:
        raise ValueError(f"Formato de arquivo não suportado: {file_path.suffix}")
    
    if not sequences:
        raise ValueError("Arquivo vazio ou sem sequências válidas")
    
    # Validação de comprimento uniforme
    L = len(sequences[0])
    valid_sequences = []
    for i, seq in enumerate(sequences):
        if len(seq) == L:
            valid_sequences.append(seq)
        else:
            logger.warning(f"Sequência {i+1} tem comprimento {len(seq)}, esperado {L} - ignorada")
    
    if not valid_sequences:
        raise ValueError("Nenhuma sequência com comprimento uniforme encontrada")
    
    if len(valid_sequences) < len(sequences):
        console.print(f"⚠️ {len(sequences) - len(valid_sequences)} sequências ignoradas por comprimento diferente")
    
    used_params = {
        'filepath': str(file_path.absolute()),
        'n': len(valid_sequences),
        'L': L,
        'format': file_path.suffix.lower(),
        'original_count': len(sequences)
    }
    
    return valid_sequences, used_params
