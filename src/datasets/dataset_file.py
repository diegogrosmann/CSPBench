"""
Módulo de Carregamento de Datasets de Arquivos - CSPBench

Este módulo fornece funcionalidades para carregar datasets de arquivos em formato FASTA
ou texto simples, com validação automática, detecção de formato e tratamento de erros.
É parte da infraestrutura de dados do CSPBench para problemas de Shortest Common
Superstring (SCS).

Arquitetura:
    O módulo implementa um sistema de carregamento de dados com:
    - Detecção automática de formato (FASTA, texto)
    - Validação de sequências e uniformidade de comprimento
    - Interface interativa para seleção de arquivos
    - Tratamento robusto de erros
    - Integração com sistema de configuração

Formatos Suportados:
    - FASTA (.fasta, .fa): Formato padrão de sequências biológicas
    - Texto simples (.txt): Uma sequência por linha
    - Detecção automática baseada em conteúdo

Funcionalidades:
    - Carregamento interativo com menu de seleção
    - Carregamento programático com parâmetros
    - Validação de formato e conteúdo
    - Normalização de sequências (maiúsculas)
    - Filtragem por comprimento uniforme
    - Relatórios de validação detalhados

Exemplo de Uso:
    ```python
    # Carregamento interativo
    sequences, params = load_dataset()

    # Carregamento programático
    sequences, params = load_dataset_with_params({
        'filepath': 'data/sequences.fasta'
    })

    print(f"Carregadas {len(sequences)} sequências")
    print(f"Comprimento: {params['L']}")
    print(f"Formato: {params['format']}")
    ```

Autor: CSPBench Development Team
Data: 2024
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Tuple

from src.core.io import get_file_info, load_sequences, validate_sequences
from src.datasets.dataset_utils import ensure_datasets_folder
from src.ui.cli.console_manager import console
from src.utils.config import FILE_DEFAULTS

logger = logging.getLogger(__name__)


def load_dataset(silent: bool = False) -> Tuple[List[str], Dict[str, Any]]:
    """
    Carrega dataset de arquivo com interface interativa ou silenciosa.

    Esta função fornece uma interface amigável para carregar sequências de arquivos,
    com detecção automática de arquivos disponíveis e validação completa.

    Args:
        silent (bool): Se True, usa configuração padrão sem interação.
                      Se False, apresenta menu interativo de seleção.

    Returns:
        Tuple[List[str], Dict[str, Any]]: Tupla contendo:
            - Lista de sequências normalizadas (maiúsculas)
            - Dicionário com parâmetros utilizados

    Raises:
        FileNotFoundError: Se o arquivo especificado não existe
        ValueError: Se o arquivo está vazio ou sequências inválidas
        Exception: Para outros erros de carregamento

    Exemplo:
        ```python
        # Carregamento interativo
        sequences, params = load_dataset()

        # Carregamento silencioso
        sequences, params = load_dataset(silent=True)

        print(f"Carregadas {len(sequences)} sequências")
        print(f"Arquivo: {params['filepath']}")
        ```

    Nota:
        - Busca automaticamente arquivos em saved_datasets/
        - Suporta extensões: .fasta, .fa, .txt
        - Normaliza sequências para maiúsculas
        - Fornece feedback visual do processo
    """
    defaults = FILE_DEFAULTS
    datasets_path = ensure_datasets_folder()

    # Buscar arquivos disponíveis nos formatos suportados
    available_files = (
        list(datasets_path.glob("*.fasta"))
        + list(datasets_path.glob("*.fa"))
        + list(datasets_path.glob("*.txt"))
    )

    if silent:
        # Modo silencioso: usar configuração padrão
        path_str = defaults["filepath"]
    else:
        # Modo interativo: apresentar opções ao usuário
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
                # Seleção por número
                file_index = int(path_input) - 1
                if 0 <= file_index < len(available_files):
                    path_str = str(available_files[file_index])
                    print(f"Arquivo selecionado: {available_files[file_index].name}")
                else:
                    print(f"Número inválido. Usando padrão: {defaults['filepath']}")
                    path_str = defaults["filepath"]
            elif path_input:
                # Caminho personalizado
                path_str = path_input
            else:
                # Configuração padrão
                path_str = defaults["filepath"]
        else:
            # Nenhum arquivo encontrado: solicitar caminho
            print("\nNenhum arquivo encontrado na pasta saved_datasets/")
            path_input = input(
                f"Caminho do arquivo (.txt ou .fasta) [{defaults['filepath']}]: "
            ).strip()
            path_str = path_input if path_input else defaults["filepath"]

    path = Path(path_str)
    logger.debug("Tentando carregar dataset de '%s'", path)

    # Verificar se o arquivo existe
    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")

    # Carregar sequências usando o parser unificado
    try:
        seqs = load_sequences(path)
        # Normalizar para maiúsculas
        seqs = [seq.upper() for seq in seqs]
    except Exception as e:
        logger.error("Erro ao carregar dataset: %s", e)
        raise

    return seqs, {"filepath": path_str}


def load_dataset_with_params(
    params: Dict[str, Any],
) -> Tuple[List[str], Dict[str, Any]]:
    """
    Carrega dataset de arquivo com parâmetros específicos.

    Esta função permite carregamento programático de datasets com configuração
    personalizada, incluindo validação completa e tratamento de erros.

    Args:
        params (Dict[str, Any]): Dicionário com parâmetros de configuração.
                                Campos suportados:
                                - filepath (str): Caminho para o arquivo
                                - Outros campos são mesclados com defaults

    Returns:
        Tuple[List[str], Dict[str, Any]]: Tupla contendo:
            - Lista de sequências validadas e normalizadas
            - Dicionário com parâmetros utilizados e metadados:
                - filepath: Caminho absoluto do arquivo
                - n: Número de sequências carregadas
                - L: Comprimento das sequências
                - format: Formato detectado do arquivo
                - original_count: Número original de sequências
                - alphabet: Alfabeto detectado

    Raises:
        FileNotFoundError: Se o arquivo especificado não existe
        ValueError: Se arquivo vazio ou sequências inválidas
        Exception: Para outros erros de carregamento/validação

    Exemplo:
        ```python
        # Carregamento básico
        sequences, params = load_dataset_with_params({
            'filepath': 'data/sequences.fasta'
        })

        # Verificar resultados
        print(f"Carregadas: {params['n']} sequências")
        print(f"Comprimento: {params['L']}")
        print(f"Formato: {params['format']}")
        print(f"Alfabeto: {params['alphabet']}")

        # Usar sequências
        for i, seq in enumerate(sequences[:5]):
            print(f"Seq {i+1}: {seq}")
        ```

    Funcionalidades:
        - Mesclagem com configurações padrão
        - Validação de formato e conteúdo
        - Filtragem por comprimento uniforme
        - Detecção automática de formato
        - Relatórios detalhados de validação
        - Tratamento robusto de erros

    Nota:
        - Sequências são normalizadas para maiúsculas
        - Apenas sequências de comprimento uniforme são mantidas
        - Erros de validação são registrados mas não impedem carregamento
        - Metadados completos são retornados para análise
    """
    # Mesclar parâmetros com configurações padrão
    merged_params = {**FILE_DEFAULTS}
    merged_params.update(params)

    filepath = merged_params["filepath"]

    console.print(f"Carregando dataset de arquivo: {filepath}")

    file_path = Path(filepath)
    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {filepath}")

    # Carregar sequências usando parsers unificados
    try:
        sequences = load_sequences(file_path)
        # Normalizar para maiúsculas
        sequences = [seq.upper() for seq in sequences]
    except Exception as e:
        logger.error("Erro ao carregar dataset: %s", e)
        raise

    if not sequences:
        raise ValueError("Arquivo vazio ou sem sequências válidas")

    # Validação usando o validador unificado
    validation = validate_sequences(sequences)

    if not validation["valid"]:
        # Registrar erros mas continuar processamento
        for error in validation["errors"]:
            logger.warning("Validação: %s", error)

    # Validação de comprimento uniforme
    if not validation["uniform_length"]:
        L = len(sequences[0])
        valid_sequences = []

        for i, seq in enumerate(sequences):
            if len(seq) == L:
                valid_sequences.append(seq)
            else:
                logger.warning(
                    f"Sequência {i+1} tem comprimento {len(seq)}, esperado {L} - ignorada"
                )

        if not valid_sequences:
            raise ValueError("Nenhuma sequência com comprimento uniforme encontrada")

        if len(valid_sequences) < len(sequences):
            console.print(
                f"⚠️ {len(sequences) - len(valid_sequences)} sequências ignoradas por comprimento diferente"
            )

        sequences = valid_sequences

    # Aplicar limite de sequências se especificado
    n_sequences = params.get("n_sequences", -1)
    if n_sequences > 0 and len(sequences) > n_sequences:
        sequences = sequences[:n_sequences]
        logger.info(f"Limitando a {n_sequences} sequências conforme solicitado")

    # Obter informações detalhadas do arquivo
    file_info = get_file_info(file_path)

    # Preparar parâmetros de retorno com metadados completos
    used_params = {
        "filepath": str(file_path.absolute()),
        "n": len(sequences),
        "n_sequences": len(sequences),  # Alias para compatibilidade
        "L": len(sequences[0]) if sequences else 0,
        "format": file_info["detected_format"],
        "original_count": validation["count"],
        "alphabet": validation["alphabet"],
    }

    return sequences, used_params
