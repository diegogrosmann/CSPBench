"""
Módulo de entrada/saída (I/O).

Contém funcionalidades para exportação de dados, formatação de resultados,
gerenciamento de arquivos e parsers unificados para datasets.
"""

from .parsers import (
    get_alphabet_from_sequences,
    get_file_info,
    load_fasta,
    load_sequences,
    load_txt,
    validate_sequences,
)

__all__ = [
    "load_fasta",
    "load_txt",
    "load_sequences",
    "validate_sequences",
    "get_alphabet_from_sequences",
    "get_file_info",
]
