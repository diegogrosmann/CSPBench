"""
Parsers unificados para I/O de datasets.

Módulo centralizado para parsing de diferentes formatos de arquivos de datasets,
fornecendo uma interface consistente para leitura de dados.

Formatos suportados:
- FASTA (.fasta, .fa, .fas)
- Texto simples (.txt)
- Arquivos customizados

Funções:
    load_fasta(file_path): Carrega sequências de um arquivo FASTA
    load_txt(file_path): Carrega sequências de um arquivo de texto
    load_sequences(file_path): Detecta formato e carrega automaticamente
    validate_sequences(sequences): Valida sequências carregadas
    get_alphabet_from_sequences(sequences): Extrai alfabeto das sequências
"""

import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def load_fasta(file_path: str | Path) -> list[str]:
    """
    Carrega sequências de um arquivo FASTA.

    Args:
        file_path: Caminho para o arquivo FASTA

    Returns:
        Lista de sequências (sem cabeçalhos)

    Raises:
        FileNotFoundError: Se o arquivo não existir
        ValueError: Se o formato FASTA estiver inválido
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

    sequences = []
    current_seq = ""
    line_num = 0

    logger.info(f"Carregando FASTA: {file_path}")

    try:
        with open(file_path, encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    # Cabeçalho FASTA - salvar sequência anterior se existir
                    if current_seq:
                        sequences.append(current_seq)
                        current_seq = ""
                else:
                    # Linha de sequência
                    current_seq += line

            # Adicionar última sequência
            if current_seq:
                sequences.append(current_seq)

    except Exception as e:
        raise ValueError(f"Erro ao processar FASTA na linha {line_num}: {e}")

    if not sequences:
        raise ValueError(f"Nenhuma sequência encontrada no arquivo: {file_path}")

    logger.info(f"Carregadas {len(sequences)} sequências do FASTA")
    return sequences


def load_txt(file_path: str | Path) -> list[str]:
    """
    Carrega sequências de um arquivo de texto simples.

    Formatos suportados:
    - Uma sequência por linha
    - Linhas vazias são ignoradas
    - Comentários iniciados com '#' são ignorados

    Args:
        file_path: Caminho para o arquivo de texto

    Returns:
        Lista de sequências

    Raises:
        FileNotFoundError: Se o arquivo não existir
        ValueError: Se nenhuma sequência válida for encontrada
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

    sequences = []
    line_num = 0

    logger.info(f"Carregando TXT: {file_path}")

    try:
        with open(file_path, encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # Ignorar linhas vazias e comentários
                if not line or line.startswith("#"):
                    continue

                sequences.append(line)

    except Exception as e:
        raise ValueError(f"Erro ao processar TXT na linha {line_num}: {e}")

    if not sequences:
        raise ValueError(f"Nenhuma sequência encontrada no arquivo: {file_path}")

    logger.info(f"Carregadas {len(sequences)} sequências do TXT")
    return sequences


def load_sequences(file_path: str | Path) -> list[str]:
    """
    Detecta o formato do arquivo e carrega as sequências automaticamente.

    Args:
        file_path: Caminho para o arquivo

    Returns:
        Lista de sequências

    Raises:
        FileNotFoundError: Se o arquivo não existir
        ValueError: Se o formato não for suportado ou inválido
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

    # Detectar formato pela extensão
    suffix = file_path.suffix.lower()

    if suffix in [".fasta", ".fa", ".fas"]:
        return load_fasta(file_path)
    elif suffix in [".txt"]:
        return load_txt(file_path)
    else:
        # Tentar detectar pelo conteúdo
        try:
            with open(file_path, encoding="utf-8") as f:
                first_line = f.readline().strip()

                if first_line.startswith(">"):
                    logger.info(f"Detectado formato FASTA em: {file_path}")
                    return load_fasta(file_path)
                else:
                    logger.info(f"Assumindo formato TXT em: {file_path}")
                    return load_txt(file_path)

        except Exception as e:
            raise ValueError(f"Não foi possível detectar formato do arquivo: {e}")


def validate_sequences(sequences: list[str]) -> dict[str, Any]:
    """
    Valida uma lista de sequências.

    Args:
        sequences: Lista de sequências para validar

    Returns:
        Dicionário com informações da validação:
        - valid: bool - Se todas as sequências são válidas
        - count: int - Número de sequências
        - lengths: list[int] - Comprimentos das sequências
        - uniform_length: bool - Se todas têm o mesmo comprimento
        - alphabet: str - Alfabeto detectado
        - errors: list[str] - Lista de erros encontrados
    """
    if not sequences:
        return {
            "valid": False,
            "count": 0,
            "lengths": [],
            "uniform_length": True,
            "alphabet": "",
            "errors": ["Lista de sequências vazia"],
        }

    errors = []
    lengths = []
    all_chars = set()

    for i, seq in enumerate(sequences):
        if not seq:
            errors.append(f"Sequência {i+1} está vazia")
            continue

        lengths.append(len(seq))
        all_chars.update(seq.upper())

    # Verificar se todas têm o mesmo comprimento
    uniform_length = len(set(lengths)) <= 1

    # Determinar alfabeto
    alphabet = "".join(sorted(all_chars))

    # Verificar se há caracteres inválidos comuns
    valid_chars = set("ACGTUNRYSWKMBDHV-")  # DNA/RNA + ambiguos + gap
    invalid_chars = all_chars - valid_chars

    if invalid_chars:
        errors.append(f"Caracteres inválidos encontrados: {invalid_chars}")

    if not uniform_length:
        errors.append(f"Sequências têm comprimentos diferentes: {set(lengths)}")

    return {
        "valid": len(errors) == 0,
        "count": len(sequences),
        "lengths": lengths,
        "uniform_length": uniform_length,
        "alphabet": alphabet,
        "errors": errors,
    }


def get_alphabet_from_sequences(sequences: list[str]) -> str:
    """
    Extrai o alfabeto único das sequências.

    Args:
        sequences: Lista de sequências

    Returns:
        String com os caracteres únicos encontrados, ordenados
    """
    if not sequences:
        return ""

    all_chars = set()
    for seq in sequences:
        all_chars.update(seq.upper())

    return "".join(sorted(all_chars))


def get_file_info(file_path: str | Path) -> dict[str, Any]:
    """
    Obtém informações sobre um arquivo de dataset.

    Args:
        file_path: Caminho para o arquivo

    Returns:
        Dicionário com informações do arquivo:
        - path: Path - Caminho do arquivo
        - exists: bool - Se o arquivo existe
        - size: int - Tamanho em bytes
        - extension: str - Extensão do arquivo
        - detected_format: str - Formato detectado
    """
    file_path = Path(file_path)

    info = {
        "path": file_path,
        "exists": file_path.exists(),
        "size": 0,
        "extension": file_path.suffix.lower(),
        "detected_format": "unknown",
    }

    if file_path.exists():
        info["size"] = file_path.stat().st_size

        # Detectar formato
        if info["extension"] in [".fasta", ".fa", ".fas"]:
            info["detected_format"] = "fasta"
        elif info["extension"] in [".txt"]:
            info["detected_format"] = "txt"
        else:
            # Tentar detectar pelo conteúdo
            try:
                with open(file_path, encoding="utf-8") as f:
                    first_line = f.readline().strip()
                    if first_line.startswith(">"):
                        info["detected_format"] = "fasta"
                    else:
                        info["detected_format"] = "txt"
            except Exception:
                info["detected_format"] = "unknown"

    return info
