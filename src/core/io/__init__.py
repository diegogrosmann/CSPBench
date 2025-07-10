"""
Módulo de Entrada/Saída (I/O) - CSPBench

Este módulo fornece um sistema unificado e robusto para entrada/saída de dados,
incluindo parsers especializados, validação automática e formatação de resultados.
É o ponto central para todas as operações de dados no CSPBench.

Funcionalidades Principais:
    - Parsers unificados para múltiplos formatos
    - Validação automática de dados
    - Detecção inteligente de formato
    - Normalização e processamento
    - Exportação em diversos formatos
    - Tratamento robusto de erros

Parsers Suportados:
    - FASTA: Formato padrão para sequências biológicas
    - Texto: Sequências simples, uma por linha
    - Detecção automática baseada em conteúdo
    - Validação de integridade e formato

Validação de Dados:
    - Verificação de formato e estrutura
    - Detecção automática de alfabeto
    - Validação de comprimento uniforme
    - Identificação de caracteres inválidos
    - Relatórios detalhados de problemas

Exemplo de Uso:
    ```python
    from src.core.io import (
        load_sequences,
        validate_sequences,
        get_file_info,
        get_alphabet_from_sequences
    )

    # Carregar e validar sequências
    sequences = load_sequences('data/sequences.fasta')
    validation = validate_sequences(sequences)

    # Analisar arquivo
    file_info = get_file_info('data/sequences.fasta')
    alphabet = get_alphabet_from_sequences(sequences)

    print(f"Formato: {file_info['detected_format']}")
    print(f"Válido: {validation['valid']}")
    print(f"Alfabeto: {alphabet}")
    ```

Integração:
    - Sistema de datasets unificado
    - Validação automática em pipelines
    - Suporte a múltiplos formatos
    - Logging detalhado de operações
    - Tratamento de erros padronizado

Autor: CSPBench Development Team
Data: 2024
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
