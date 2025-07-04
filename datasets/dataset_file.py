"""
Leitura de arquivos de texto ou FASTA para extração de sequências.

Funções:
    load_dataset(): Lê um arquivo e retorna as sequências e parâmetros.
"""

import logging
from pathlib import Path
from typing import Any

from csp_blfga.core.io import get_file_info, load_sequences, validate_sequences
from csp_blfga.ui.cli.console_manager import console
from csp_blfga.utils.config import FILE_DEFAULTS
from datasets.dataset_utils import ensure_datasets_folder

logger = logging.getLogger(__name__)


def load_dataset(silent: bool = False) -> tuple[list[str], dict[str, Any]]:
    defaults = FILE_DEFAULTS
    datasets_path = ensure_datasets_folder()
    available_files = (
        list(datasets_path.glob("*.fasta")) + list(datasets_path.glob("*.fa")) + list(datasets_path.glob("*.txt"))
    )
    if silent:
        path_str = defaults["filepath"]
    else:
        if available_files:
            print("\nArquivos disponíveis na pasta saved_datasets/:")
            for i, file in enumerate(available_files, 1):
                print(f"  {i}. {file.name}")
            print("\nOpções:")
            print("  - Digite o número do arquivo para selecioná-lo")
            print("  - Digite o caminho completo do arquivo")
            print(f"  - Pressione Enter para usar o padrão [{defaults['filepath']}]")
            path_input = input("\nSua escolha: ").strip()
            if path_input.isdigit():
                file_index = int(path_input) - 1
                if 0 <= file_index < len(available_files):
                    path_str = str(available_files[file_index])
                    print(f"Arquivo selecionado: {available_files[file_index].name}")
                else:
                    print(f"Número inválido. Usando padrão: {defaults['filepath']}")
                    path_str = defaults["filepath"]
            elif path_input:
                path_str = path_input
            else:
                path_str = defaults["filepath"]
        else:
            print("\nNenhum arquivo encontrado na pasta saved_datasets/")
            path_input = input(f"Caminho do arquivo (.txt ou .fasta) [{defaults['filepath']}]: ").strip()
            path_str = path_input if path_input else defaults["filepath"]
    path = Path(path_str)
    logger.debug(f"Tentando carregar dataset de '{path}'")

    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")

    # Usar o parser unificado
    try:
        seqs = load_sequences(path)
        # Converter para maiúsculas
        seqs = [seq.upper() for seq in seqs]
    except Exception as e:
        logger.error(f"Erro ao carregar dataset: {e}")
        raise

    return seqs, {"filepath": path_str}


def load_dataset_with_params(
    params: dict[str, Any],
) -> tuple[list[str], dict[str, Any]]:
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

    filepath = merged_params["filepath"]

    console.print(f"Carregando dataset de arquivo: {filepath}")

    file_path = Path(filepath)
    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {filepath}")

    # Usar parsers unificados
    try:
        sequences = load_sequences(file_path)
        # Converter para maiúsculas
        sequences = [seq.upper() for seq in sequences]
    except Exception as e:
        logger.error(f"Erro ao carregar dataset: {e}")
        raise

    if not sequences:
        raise ValueError("Arquivo vazio ou sem sequências válidas")

    # Validação usando o validador unificado
    validation = validate_sequences(sequences)

    if not validation["valid"]:
        # Log dos erros mas continua processamento
        for error in validation["errors"]:
            logger.warning(f"Validação: {error}")

    # Validação de comprimento uniforme
    if not validation["uniform_length"]:
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

        sequences = valid_sequences

    # Obter informações do arquivo
    file_info = get_file_info(file_path)

    used_params = {
        "filepath": str(file_path.absolute()),
        "n": len(sequences),
        "L": len(sequences[0]) if sequences else 0,
        "format": file_info["detected_format"],
        "original_count": validation["count"],
        "alphabet": validation["alphabet"],
    }

    return sequences, used_params
