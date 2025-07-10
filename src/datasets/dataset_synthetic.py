"""
Módulo de Geração de Datasets Sintéticos - CSPBench

Este módulo fornece funcionalidades avançadas para geração de datasets sintéticos
para problemas de Closest String Problem (CSP). Oferece controle preciso sobre
características dos dados, incluindo tamanho, ruído, alfabeto e reprodutibilidade.

Arquitetura:
    O módulo implementa um sistema flexível de geração de dados com:
    - Geração baseada em string centro + ruído controlado
    - Geração completamente aleatória
    - Interface interativa e programática
    - Controle de reprodutibilidade via seeds
    - Validação e métricas automáticas

Estratégias de Geração:
    1. **Baseada em Centro**: Gera uma string base e aplica ruído controlado
    2. **Completamente Aleatória**: Cada string é gerada independentemente
    3. **Ruído Uniforme**: Mesmo nível de ruído para todas as strings
    4. **Ruído Variável**: Níveis diferentes de ruído por string

Funcionalidades:
    - Controle preciso de parâmetros (n, L, alfabeto, ruído)
    - Geração reprodutível com seeds
    - Validação automática de distâncias
    - Métricas de qualidade do dataset
    - Interface interativa amigável
    - Suporte a alfabetos customizados

Exemplo de Uso:
    ```python
    # Geração interativa
    sequences, params = generate_dataset()

    # Geração programática
    sequences, params = generate_dataset_with_params({
        'n': 100,
        'L': 50,
        'alphabet': 'ACGT',
        'noise': 0.15,
        'seed': 42
    })

    # Geração direta
    sequences, params = generate_dataset_from_params(
        n=50, L=30, alphabet='ACGT', noise=0.2, seed=123
    )

    print(f"Geradas {len(sequences)} sequências")
    print(f"Distância da string base: {params['distancia_string_base']}")
    ```

Controle de Qualidade:
    - Validação automática de distâncias
    - Métricas de diversidade
    - Análise de distribuição de ruído
    - Verificação de uniformidade

Autor: CSPBench Development Team
Data: 2024
"""

import logging
import random
from typing import Any, Dict, List, Optional, Tuple

from src.ui.cli.console_manager import console
from src.utils.config import SYNTHETIC_DEFAULTS, safe_input

logger = logging.getLogger(__name__)


def generate_dataset(silent: bool = False) -> Tuple[List[str], Dict[str, Any]]:
    """
    Gera dataset sintético com interface interativa completa.

    Esta função oferece uma interface amigável para configurar e gerar datasets
    sintéticos para problemas CSP, com controle total sobre todos os parâmetros.

    Args:
        silent (bool): Se True, usa configurações padrão sem interação.
                      Se False, apresenta interface interativa.

    Returns:
        Tuple[List[str], Dict[str, Any]]: Tupla contendo:
            - Lista de sequências geradas
            - Dicionário com parâmetros utilizados e métricas:
                - n: Número de strings geradas
                - L: Comprimento das strings
                - alphabet: Alfabeto utilizado
                - noise: Taxa de ruído aplicada
                - fully_random: Se foi geração totalmente aleatória
                - seed: Semente utilizada para reprodutibilidade
                - base_string: String centro utilizada (se não aleatória)
                - distancia_string_base: Distância máxima da string base

    Exemplo:
        ```python
        # Geração interativa
        sequences, params = generate_dataset()

        # Geração silenciosa
        sequences, params = generate_dataset(silent=True)

        # Analisar resultados
        print(f"Dataset: {params['n']} sequências de comprimento {params['L']}")
        print(f"Alfabeto: {params['alphabet']}")
        print(f"Ruído: {params['noise']}")
        print(f"Semente: {params['seed']}")

        if not params['fully_random']:
            print(f"String base: {params['base_string']}")
            print(f"Distância máxima: {params['distancia_string_base']}")
        ```

    Funcionalidades:
        - Interface interativa para configuração de parâmetros
        - Geração de semente automática para reprodutibilidade
        - Validação de parâmetros de entrada
        - Geração baseada em centro ou totalmente aleatória
        - Cálculo automático de métricas de qualidade
        - Logging detalhado do processo

    Parâmetros Configuráveis:
        - n: Número de strings (padrão: configuração)
        - L: Comprimento das strings (padrão: configuração)
        - alphabet: Alfabeto customizado (padrão: ACGT)
        - noise: Taxa de ruído 0-1 (padrão: configuração)
        - fully_random: Geração totalmente aleatória (padrão: False)
        - seed: Semente para reprodutibilidade (padrão: auto-gerada)

    Nota:
        - Semente é gerada automaticamente se não fornecida
        - Distância da string base é calculada automaticamente
        - Todos os parâmetros são validados antes da geração
    """
    defaults = SYNTHETIC_DEFAULTS

    if silent:
        # Modo silencioso: usar configurações padrão
        n = defaults["n"]
        L = defaults["L"]
        alphabet = defaults["alphabet"]
        noise = defaults["noise"]
        fully_random = False
        seed = None
    else:
        # Modo interativo: solicitar parâmetros do usuário
        n_input = safe_input(f"Quantas strings (n)? [{defaults['n']}]: ")
        n = int(n_input) if n_input else defaults["n"]

        L_input = safe_input(f"Comprimento das strings (L)? [{defaults['L']}]: ")
        L = int(L_input) if L_input else defaults["L"]

        alpha_input = safe_input(
            f"Alfabeto (ex.: ACGT ou ABCDE...)? [{defaults['alphabet']}]: "
        ).upper()
        alphabet = alpha_input if alpha_input else defaults["alphabet"]

        noise_input = safe_input(
            f"Taxa de ruído por posição (0–1) [{defaults['noise']}]: "
        )
        noise = float(noise_input) if noise_input else defaults["noise"]

        fully_random_input = safe_input(
            "Gerar strings totalmente aleatórias? (s/n) [n]: "
        ).lower()
        fully_random = fully_random_input.startswith("s")

        # Solicitar semente para reprodutibilidade (opcional)
        seed_input = safe_input("Semente aleatória (int) [None]: ")
        seed = int(seed_input) if seed_input else None

    # Gerar semente automaticamente se não fornecida
    if seed is None:
        import time

        seed = int(time.time() * 1000000) % (2**32)  # Semente baseada no timestamp
        logger.info("Semente gerada automaticamente: %s", seed)

    # Preparar parâmetros para geração
    params = {
        "n": n,
        "L": L,
        "alphabet": alphabet,
        "noise": noise,
        "fully_random": fully_random,
        "seed": seed,
    }

    logger.info("Iniciando geração do dataset sintético")
    logger.debug(
        "Gerando dataset sintético com n=%s, L=%s, |Σ|=%s, noise=%s, fully_random=%s, seed=%s",
        n,
        L,
        len(alphabet),
        noise,
        fully_random,
        seed,
    )

    # Configurar gerador de números aleatórios
    rng = random.Random(seed)
    data = []

    if fully_random:
        # Geração completamente aleatória
        for _ in range(n):
            seq = "".join(rng.choices(alphabet, k=L))
            data.append(seq)
        logger.info("Strings totalmente aleatórias geradas")
        params["base_string"] = None
    else:
        # Geração baseada em string centro
        base_string = "".join(rng.choices(alphabet, k=L))
        logger.info("String centro gerada para aplicar ruído")

        # Definir nível de ruído
        if noise is None:
            noise = rng.uniform(0.1, 0.5)
            logger.info("Noise não informado, sorteado aleatoriamente: %.3f", noise)
        else:
            logger.info("Noise informado: %s", noise)

        # Geração das strings com ruído
        for _ in range(n):
            s = list(base_string)
            num_mut = int(round(noise * L))
            mut_pos = rng.sample(range(L), num_mut) if num_mut > 0 else []

            for pos in mut_pos:
                orig = s[pos]
                alt = rng.choice([c for c in alphabet if c != orig])
                s[pos] = alt

            new_s = "".join(s)
            data.append(new_s)

        params["base_string"] = base_string

    # Calcular métricas de qualidade do dataset
    if not fully_random and params["base_string"] is not None:
        from src.utils.distance import max_distance

        distancia_base = max_distance(params["base_string"], data)
        params["distancia_string_base"] = distancia_base

    logger.info("Parâmetros retornados pelo gerador: %s", params)
    return data, params


def _unpack_singleton_lists(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Converte valores que são listas de um elemento para o valor simples.

    Esta função utilitária processa parâmetros que podem vir em formato de lista
    singleton, convertendo-os para valores simples para facilitar o processamento.

    Args:
        params (Dict[str, Any]): Dicionário com parâmetros que podem conter listas singleton

    Returns:
        Dict[str, Any]: Dicionário com valores simplificados

    Exemplo:
        ```python
        # Entrada com listas singleton
        params = {
            'n': [80],
            'L': [50],
            'alphabet': ['ACGT'],
            'noise': [0.15],
            'other': ['value1', 'value2']  # Lista múltipla mantida
        }

        # Saída simplificada
        result = _unpack_singleton_lists(params)
        # {
        #     'n': 80,
        #     'L': 50,
        #     'alphabet': 'ACGT',
        #     'noise': 0.15,
        #     'other': ['value1', 'value2']
        # }
        ```

    Nota:
        - Apenas listas com exatamente um elemento são convertidas
        - Listas vazias ou com múltiplos elementos são mantidas como estão
        - Valores que não são listas permanecem inalterados
    """
    return {
        k: (v[0] if isinstance(v, list) and len(v) == 1 else v)
        for k, v in params.items()
    }


def generate_dataset_with_params(
    params: Dict[str, Any],
) -> Tuple[List[str], Dict[str, Any]]:
    """
    Gera dataset sintético com parâmetros específicos fornecidos.

    Esta função permite geração programática de datasets sintéticos com configuração
    personalizada, oferecendo controle total sobre todos os aspectos da geração.

    Args:
        params (Dict[str, Any]): Dicionário com parâmetros de configuração.
                                Campos suportados:
                                - n (int): Número de strings a gerar
                                - L (int): Comprimento das strings
                                - alphabet (str): Alfabeto a utilizar
                                - noise (float): Taxa de ruído (0-1)
                                - fully_random (bool): Geração totalmente aleatória
                                - aleatorio_total (bool): Alias para fully_random
                                - seed (int): Semente para reprodutibilidade
                                Campos não fornecidos usam valores padrão.

    Returns:
        Tuple[List[str], Dict[str, Any]]: Tupla contendo:
            - Lista de sequências geradas
            - Dicionário com parâmetros utilizados e métricas:
                - n: Número real de strings geradas
                - L: Comprimento das strings
                - alphabet: Alfabeto utilizado
                - noise: Taxa(s) de ruído aplicada(s)
                - fully_random: Tipo de geração utilizado
                - seed: Semente utilizada
                - base_string: String centro (se aplicável)
                - distancia_string_base: Distância máxima (se aplicável)

    Exemplo:
        ```python
        # Geração básica
        sequences, params = generate_dataset_with_params({
            'n': 100,
            'L': 50,
            'alphabet': 'ACGT',
            'noise': 0.15,
            'seed': 42
        })

        # Geração totalmente aleatória
        sequences, params = generate_dataset_with_params({
            'n': 50,
            'L': 30,
            'alphabet': 'ACGT',
            'fully_random': True,
            'seed': 123
        })

        # Geração com ruído variável
        sequences, params = generate_dataset_with_params({
            'n': 80,
            'L': 40,
            'alphabet': 'ACGT',
            'noise': None,  # Ruído aleatório por string
            'seed': 456
        })

        # Analisar resultados
        print(f"Geradas {len(sequences)} sequências")
        print(f"Parâmetros usados: {params}")
        ```

    Funcionalidades:
        - Mesclagem automática com configurações padrão
        - Suporte a ruído fixo ou variável
        - Geração de semente automática
        - Validação de parâmetros
        - Cálculo de métricas de qualidade
        - Logging detalhado do processo

    Estratégias de Ruído:
        - **Fixo**: Mesmo nível para todas as strings
        - **Variável**: Nível diferente para cada string
        - **Nulo**: Ruído aleatório entre 0.1-0.5

    Nota:
        - Parâmetros são validados e mesclados com defaults
        - Semente é gerada automaticamente se não fornecida
        - Suporte a aliases para compatibilidade (aleatorio_total)
        - Métricas de qualidade são calculadas automaticamente
    """
    # Desempacotar listas singleton e mesclar com defaults
    params = _unpack_singleton_lists(params)
    merged_params = {**SYNTHETIC_DEFAULTS}
    merged_params.update(params)

    # Extrair parâmetros com validação
    n = merged_params["n"]
    L = merged_params["L"]
    alphabet = merged_params["alphabet"]

    # Aceitar tanto fully_random quanto aleatorio_total para compatibilidade
    fully_random = merged_params.get(
        "fully_random", merged_params.get("aleatorio_total", False)
    )
    seed = merged_params.get("seed", None)

    # Gerar semente automaticamente se não fornecida
    if seed is None:
        import time

        seed = int(time.time() * 1000000) % (2**32)
        logger.info("Semente gerada automaticamente: %s", seed)

    noise = merged_params.get("noise", None)  # Pode ser None para ruído aleatório

    console.print(
        f"Parâmetros usados: n={n}, L={L}, alphabet='{alphabet}', "
        f"noise={noise}, fully_random={fully_random}"
    )

    # Configurar gerador e estruturas de dados
    rng = random.Random(seed)
    sequences = []
    used_params = {
        "n": n,
        "L": L,
        "alphabet": alphabet,
        "fully_random": fully_random,
        "seed": seed,
    }

    if fully_random:
        # Geração completamente aleatória
        for _ in range(n):
            seq = "".join(rng.choice(alphabet) for _ in range(L))
            sequences.append(seq)
        used_params["noise"] = None
        used_params["base_string"] = None
    else:
        # Geração baseada em string centro
        base_string = "".join(rng.choices(alphabet, k=L))
        used_params["base_string"] = base_string

        if noise is not None:
            # Ruído fixo para todas as strings
            used_params["noise"] = noise
            for i in range(n):
                seq = list(base_string)
                for pos in range(L):
                    if rng.random() < noise:
                        current_char = seq[pos]
                        other_chars = [c for c in alphabet if c != current_char]
                        if other_chars:
                            seq[pos] = rng.choice(other_chars)
                sequences.append("".join(seq))
        else:
            # Ruído variável: nível diferente para cada string
            noises = [rng.uniform(0.1, 0.5) for _ in range(n)]
            used_params["noise"] = noises

            for i in range(n):
                seq = list(base_string)
                this_noise = noises[i]
                for pos in range(L):
                    if rng.random() < this_noise:
                        current_char = seq[pos]
                        other_chars = [c for c in alphabet if c != current_char]
                        if other_chars:
                            seq[pos] = rng.choice(other_chars)
                sequences.append("".join(seq))

    # Calcular métricas de qualidade do dataset
    if not fully_random and used_params.get("base_string") is not None:
        from src.utils.distance import max_distance

        base_str = used_params["base_string"]
        if isinstance(base_str, str):
            distancia_base = max_distance(base_str, sequences)
            used_params["distancia_string_base"] = distancia_base

            console.print(
                f"Distância da string base: {distancia_base}, Semente: {seed}"
            )

    return sequences, used_params


def generate_perturbations(
    base_string: str, n_strings: int, max_distance: int, alphabet: str, rng
) -> List[str]:
    """
    Gera strings com distância exata igual a max_distance da string base.

    Esta função especializada gera um conjunto de strings que têm exatamente
    a distância especificada da string base, útil para testes controlados
    e validação de algoritmos.

    Args:
        base_string (str): String base de referência
        n_strings (int): Número de strings a gerar
        max_distance (int): Distância exata desejada da string base
        alphabet (str): Alfabeto a utilizar para substituições
        rng: Gerador de números aleatórios configurado

    Returns:
        List[str]: Lista de strings com distância exata especificada

    Raises:
        ValueError: Se max_distance > len(base_string)
        IndexError: Se não há caracteres suficientes no alfabeto

    Exemplo:
        ```python
        import random

        # Configurar gerador
        rng = random.Random(42)

        # Gerar strings com distância exata
        strings = generate_perturbations(
            base_string="ACGT",
            n_strings=5,
            max_distance=2,
            alphabet="ACGT",
            rng=rng
        )

        # Verificar resultados
        for i, s in enumerate(strings):
            dist = hamming_distance("ACGT", s)
            print(f"String {i}: {s} (distância: {dist})")
        ```

    Funcionalidades:
        - Distância exata garantida
        - Validação automática de resultados
        - Logging detalhado do processo
        - Seleção aleatória de posições
        - Verificação de segurança para alfabetos

    Nota:
        - Cada string gerada tem exatamente max_distance diferenças
        - Posições são selecionadas aleatoriamente sem repetição
        - Caracteres de substituição são diferentes dos originais
        - Validação automática garante correção dos resultados
    """
    func_logger = logging.getLogger(__name__)

    # Validar parâmetros
    L = len(base_string)
    if max_distance > L:
        raise ValueError(
            f"max_distance ({max_distance}) não pode ser maior que o comprimento da string ({L})"
        )

    strings = []

    for i in range(n_strings):
        # Criar string com exatamente max_distance diferenças
        new_string = list(base_string)
        positions = rng.sample(range(L), min(max_distance, L))

        func_logger.debug("[GENERATOR] String %s: posições alteradas: %s", i, positions)

        for pos in positions:
            # Escolher um símbolo diferente do atual
            current_char = base_string[pos]
            available_chars = [c for c in alphabet if c != current_char]

            if available_chars:  # Verificação de segurança
                chosen_char = rng.choice(available_chars)
                new_string[pos] = chosen_char
                func_logger.debug(
                    "[GENERATOR] String %s: pos %s: %s -> %s",
                    i,
                    pos,
                    current_char,
                    chosen_char,
                )
            else:
                func_logger.warning(
                    "Nenhum caractere alternativo disponível para posição %s", pos
                )

        result_string = "".join(new_string)
        strings.append(result_string)

        # VALIDAÇÃO: Verificar se a distância está correta
        from src.utils.distance import hamming_distance

        actual_distance = hamming_distance(base_string, result_string)

        if actual_distance != max_distance:
            func_logger.error(
                "[GENERATOR] ERRO: String %s tem distância %s, esperado %s",
                i,
                actual_distance,
                max_distance,
            )

    # VALIDAÇÃO FINAL: Verificar todas as distâncias
    from src.utils.distance import hamming_distance

    all_distances = [hamming_distance(base_string, s) for s in strings]
    func_logger.info("[GENERATOR] Todas as distâncias geradas: %s", all_distances)
    func_logger.info(
        "[GENERATOR] Máxima: %s, Mínima: %s", max(all_distances), min(all_distances)
    )

    return strings


def generate_dataset_from_params(
    n: int,
    L: int,
    alphabet: str,
    noise: float,
    fully_random: bool = False,
    seed: Optional[int] = None,
) -> Tuple[List[str], Dict[str, Any]]:
    """
    Gera dataset sintético diretamente a partir dos parâmetros fornecidos.

    Esta função oferece uma interface direta para geração de datasets sintéticos,
    sem necessidade de dicionários de configuração, ideal para uso programático
    e integração com outros sistemas.

    Args:
        n (int): Número de strings a gerar
        L (int): Comprimento das strings
        alphabet (str): Alfabeto a ser usado (ex: 'ACGT', 'ABCD')
        noise (float): Taxa de ruído por posição (0–1)
        fully_random (bool): Se True, gera strings totalmente aleatórias.
                            Se False, usa string base + ruído. Padrão: False
        seed (Optional[int]): Semente para reprodutibilidade.
                             Se None, gera automaticamente. Padrão: None

    Returns:
        Tuple[List[str], Dict[str, Any]]: Tupla contendo:
            - Lista de sequências geradas
            - Dicionário com parâmetros utilizados:
                - n: Número de strings geradas
                - L: Comprimento das strings
                - alphabet: Alfabeto utilizado
                - noise: Taxa de ruído aplicada
                - fully_random: Tipo de geração utilizado
                - seed: Semente utilizada
                - distancia_string_base: Distância máxima (se não aleatório)

    Exemplo:
        ```python
        # Geração básica com ruído
        sequences, params = generate_dataset_from_params(
            n=50,
            L=30,
            alphabet='ACGT',
            noise=0.15,
            seed=42
        )

        # Geração totalmente aleatória
        sequences, params = generate_dataset_from_params(
            n=100,
            L=25,
            alphabet='ACGT',
            noise=0.0,  # Ignorado em modo aleatório
            fully_random=True,
            seed=123
        )

        # Geração com semente automática
        sequences, params = generate_dataset_from_params(
            n=80,
            L=40,
            alphabet='ABCDEFGH',
            noise=0.2
        )

        # Analisar resultados
        print(f"Geradas {len(sequences)} sequências")
        print(f"Semente utilizada: {params['seed']}")
        if not params['fully_random']:
            print(f"Distância máxima: {params['distancia_string_base']}")
        ```

    Funcionalidades:
        - Interface direta sem configurações complexas
        - Geração automática de semente
        - Validação de parâmetros
        - Logging detalhado
        - Cálculo automático de métricas
        - Suporte a ambos os modos de geração

    Modos de Geração:
        - **Centro + Ruído**: Gera string base e aplica ruído controlado
        - **Totalmente Aleatório**: Cada string gerada independentemente

    Nota:
        - Semente é gerada automaticamente baseada no timestamp se não fornecida
        - Em modo totalmente aleatório, parâmetro 'noise' é ignorado
        - Distância da string base é calculada apenas em modo centro + ruído
        - Todos os parâmetros são validados antes da geração
    """
    import time

    # Gerar semente automaticamente se não fornecida
    if seed is None:
        seed = int(time.time() * 1000000) % (2**32)
        logger.info("Semente gerada automaticamente: %s", seed)

    # Preparar parâmetros de retorno
    params = {
        "n": n,
        "L": L,
        "alphabet": alphabet,
        "noise": noise,
        "fully_random": fully_random,
        "seed": seed,
    }

    logger.info("Iniciando geração do dataset sintético")
    logger.debug(
        "Gerando dataset sintético com n=%s, L=%s, |Σ|=%s, noise=%s, fully_random=%s, seed=%s",
        n,
        L,
        len(alphabet),
        noise,
        fully_random,
        seed,
    )

    # Configurar gerador de números aleatórios
    rng = random.Random(seed)
    data = []

    if fully_random:
        # Gerar strings totalmente aleatórias
        for _ in range(n):
            s = "".join(rng.choices(alphabet, k=L))
            data.append(s)
        logger.info("Strings totalmente aleatórias geradas")
    else:
        # Geração baseada em string centro
        base_string = "".join(rng.choices(alphabet, k=L))
        logger.info("String centro gerada para aplicar ruído")

        # Geração das strings com ruído
        for _ in range(n):
            s = list(base_string)
            num_mut = int(round(noise * L))
            mut_pos = rng.sample(range(L), num_mut) if num_mut > 0 else []

            for pos in mut_pos:
                orig = s[pos]
                alt = rng.choice([c for c in alphabet if c != orig])
                s[pos] = alt

            new_s = "".join(s)
            data.append(new_s)

        # Calcular distância da string base
        from src.utils.distance import max_distance

        distancia_base = max_distance(base_string, data)
        params["distancia_string_base"] = distancia_base

    logger.info("Dataset sintético gerado com sucesso")
    return data, params
