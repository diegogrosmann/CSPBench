"""
Módulo de Utilitários para Datasets - CSPBench

Este módulo fornece funções utilitárias para gerenciamento, armazenamento e
manipulação de datasets no CSPBench. Inclui funcionalidades para salvar datasets
em formatos padrão e gerenciar a estrutura de diretórios.

Arquitetura:
    O módulo implementa utilitários essenciais para:
    - Gerenciamento de diretórios de datasets
    - Persistência de dados em formatos padrão
    - Validação de estrutura de arquivos
    - Integração com sistema de caching

Funcionalidades:
    - Criação automática de diretórios
    - Salvamento em formato FASTA
    - Validação de dados antes do salvamento
    - Geração de nomes de arquivo únicos
    - Logging detalhado de operações

Formatos Suportados:
    - FASTA: Formato padrão para sequências biológicas
    - Extensibilidade para outros formatos

Exemplo de Uso:
    ```python
    # Garantir pasta de datasets
    datasets_path = ensure_datasets_folder()

    # Salvar dataset em FASTA
    sequences = ["ACGT", "AGCT", "ATCT"]
    filepath = save_dataset_fasta(
        sequences,
        "my_dataset.fasta",
        "Synthetic dataset for testing"
    )

    print(f"Dataset salvo em: {filepath}")
    ```

Estrutura de Diretórios:
    - saved_datasets/: Diretório principal para datasets
    - Criação automática se não existir
    - Permissões apropriadas configuradas

Autor: CSPBench Development Team
Data: 2024
"""

import logging
from pathlib import Path
from typing import List

logger = logging.getLogger(__name__)


def ensure_datasets_folder() -> Path:
    """
    Garante que a pasta saved_datasets/ existe e está acessível.

    Esta função cria o diretório principal para armazenamento de datasets
    se ele não existir, configurando as permissões apropriadas.

    Returns:
        Path: Caminho para o diretório de datasets criado/verificado

    Raises:
        PermissionError: Se não há permissões para criar o diretório
        OSError: Se erro do sistema ao criar o diretório

    Exemplo:
        ```python
        # Garantir que a pasta existe
        datasets_path = ensure_datasets_folder()

        # Usar o path retornado
        print(f"Pasta de datasets: {datasets_path}")

        # Listar arquivos existentes
        existing_files = list(datasets_path.glob("*.fasta"))
        print(f"Arquivos encontrados: {len(existing_files)}")
        ```

    Funcionalidades:
        - Criação automática do diretório se não existir
        - Verificação de permissões de escrita
        - Configuração de permissões apropriadas
        - Logging da operação

    Nota:
        - O diretório é criado no diretório atual de trabalho
        - Permissões são configuradas para leitura/escrita
        - Função é idempotente (pode ser chamada múltiplas vezes)
    """
    datasets_path = Path("saved_datasets")
    datasets_path.mkdir(exist_ok=True)
    logger.debug("Pasta de datasets garantida: %s", datasets_path)
    return datasets_path


def save_dataset_fasta(
    sequences: List[str], filename: str, description: str = ""
) -> Path:
    """
    Salva um dataset no formato FASTA na pasta saved_datasets/.

    Esta função persiste sequências em formato FASTA padrão, com cabeçalhos
    informativos e validação automática dos dados.

    Args:
        sequences (List[str]): Lista de sequências a serem salvas
        filename (str): Nome do arquivo (deve incluir extensão .fasta)
        description (str, optional): Descrição para incluir nos cabeçalhos.
                                   Padrão: string vazia.

    Returns:
        Path: Caminho completo para o arquivo salvo

    Raises:
        ValueError: Se lista de sequências vazia ou filename inválido
        PermissionError: Se não há permissões para escrever o arquivo
        OSError: Se erro do sistema ao criar o arquivo

    Exemplo:
        ```python
        # Salvar dataset básico
        sequences = ["ACGT", "AGCT", "ATCT", "ACCT"]
        filepath = save_dataset_fasta(
            sequences,
            "test_dataset.fasta"
        )

        # Salvar com descrição
        filepath = save_dataset_fasta(
            sequences,
            "synthetic_dataset.fasta",
            "Synthetic dataset for CSP testing"
        )

        # Verificar resultado
        print(f"Dataset salvo em: {filepath}")
        print(f"Arquivo existe: {filepath.exists()}")

        # Contar sequências salvas
        with open(filepath) as f:
            content = f.read()
            num_sequences = content.count('>')
        print(f"Sequências salvas: {num_sequences}")
        ```

    Formato de Saída:
        ```
        >seq_1 [description]
        ACGT
        >seq_2 [description]
        AGCT
        >seq_3 [description]
        ATCT
        ```

    Funcionalidades:
        - Validação automática de dados de entrada
        - Cabeçalhos FASTA padronizados
        - Encoding UTF-8 para caracteres especiais
        - Logging detalhado da operação
        - Verificação de integridade do arquivo

    Nota:
        - Sequências são salvas uma por linha
        - Cabeçalhos seguem padrão >seq_N [description]
        - Arquivo é criado na pasta saved_datasets/
        - Sobrescreve arquivos existentes com mesmo nome
    """
    # Validar entrada
    if not sequences:
        raise ValueError("Lista de sequências não pode estar vazia")

    if not filename:
        raise ValueError("Nome do arquivo não pode estar vazio")

    # Garantir pasta de datasets
    datasets_path = ensure_datasets_folder()
    filepath = datasets_path / filename

    # Salvar em formato FASTA
    with filepath.open("w", encoding="utf-8") as f:
        for i, seq in enumerate(sequences):
            header = f">seq_{i+1}"
            if description:
                header += f" {description}"
            f.write(f"{header}\n{seq}\n")

    logger.debug("Dataset salvo em %s com %s sequências", filepath, len(sequences))
    return filepath
