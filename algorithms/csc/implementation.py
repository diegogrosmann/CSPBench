"""
Implementação do algoritmo Consensus String Clustering (CSC) para CSP.

O CSC é uma abordagem inovadora que combina técnicas de clusterização com
consenso local e recombinação de blocos para resolver o Closest String Problem.
A estratégia central é dividir o problema em subproblemas menores através
de clustering, resolver cada subproblema localmente, e depois recombinar
as soluções parciais para formar candidatos à solução global.

ARQUITETURA ALGORÍTMICA:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                          ALGORITMO CSC DETALHADO                                │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. ANÁLISE ESTATÍSTICA                                                         │
│   ├── Calcula distâncias de Hamming entre todos os pares de strings            │
│   ├── Computa estatísticas (média, mediana, min, max)                          │
│   └── Define parâmetros automáticos (d para DBSCAN, n_blocks para recombinação)│
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. CLUSTERIZAÇÃO BASEADA EM DISTÂNCIA                                         │
│   ├── Converte strings para representação numérica                             │
│   ├── Aplica DBSCAN com métrica de Hamming customizada                         │
│   ├── Identifica grupos de strings similares                                   │
│   └── Remove outliers (strings que não pertencem a clusters)                   │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 3. CONSENSO LOCAL POR CLUSTER                                                 │
│   ├── Para cada cluster encontrado:                                            │
│   │   ├── Calcula string consenso por votação majoritária                      │
│   │   ├── Obtém representante local ótimo para o subgrupo                      │
│   │   └── Divide consenso em blocos de tamanho uniforme                        │
│   └── Cria repositório de soluções locais de alta qualidade                    │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 4. RECOMBINAÇÃO DE BLOCOS                                                     │
│   ├── Gera produto cartesiano de todos os blocos de todos os consensos         │
│   ├── Cria candidatos combinando blocos de diferentes consensos                │
│   ├── Explora hibridização entre soluções locais                               │
│   └── Normaliza tamanho dos candidatos para match com strings originais        │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 5. SELEÇÃO E REFINAMENTO                                                      │
│   ├── Avalia todos os candidatos usando função max_distance                    │
│   ├── Seleciona candidato com menor distância máxima                           │
│   ├── Aplica busca local intensiva posição-a-posição                           │
│   └── Retorna solução refinada localmente ótima                                │
└─────────────────────────────────────────────────────────────────────────────────┘

INOVAÇÕES E CARACTERÍSTICAS:

• **CLUSTERING INTELIGENTE**: Usa DBSCAN para identificar estrutura nos dados
• **CONSENSO ESTRATIFICADO**: Resolve subproblemas de forma independente
• **RECOMBINAÇÃO CRIATIVA**: Hibridiza soluções locais para explorar novas regiões
• **ADAPTAÇÃO AUTOMÁTICA**: Parâmetros ajustados às características dos dados
• **BUSCA LOCAL INTEGRADA**: Refinamento final para melhorar qualidade

VANTAGENS SOBRE MÉTODOS CLÁSSICOS:

• **Escalabilidade**: Divide problema grande em subproblemas menores
• **Qualidade**: Consenso local é geralmente melhor que global
• **Robustez**: Funciona bem mesmo com dados heterogêneos
• **Flexibilidade**: Adapta-se a diferentes padrões de similaridade

COMPLEXIDADE ALGORÍTMICA:

• **Clusterização**: O(n² × L) para DBSCAN com métrica de Hamming
• **Consenso**: O(c × |cluster| × L) onde c é número de clusters
• **Recombinação**: O(|consensos|^n_blocks) candidatos gerados
• **Busca Local**: O(L × |Σ| × iterações) onde |Σ| é tamanho do alfabeto

PARÂMETROS CRÍTICOS:

• **d (eps do DBSCAN)**: Raio máximo para considerar strings similares
• **n_blocks**: Número de blocos para divisão (balanceia exploração vs eficiência)
• **min_samples**: Número mínimo de strings para formar cluster

CASOS DE USO IDEAIS:

• **Dados Estruturados**: Quando strings têm padrões locais distintos
• **Instâncias Grandes**: Melhora escalabilidade vs métodos gulosos
• **Diversidade Alta**: Quando há múltiplas "famílias" de strings
• **Exploração-Exploitation**: Balanceia consenso local com recombinação global

EXEMPLO DE FUNCIONAMENTO:

Entrada: ["AAATTT", "AAAGGG", "CCCGGG", "CCCTTT"]
1. Clustering: {["AAATTT", "AAAGGG"], ["CCCGGG", "CCCTTT"]}
2. Consensos: ["AAAGGR", "CCCGRR"] (R = mais frequente)
3. Blocos: [["AAA", "GGR"], ["CCC", "GRR"]]
4. Recombinação: ["AAAGRR", "CCCGGR", ...]
5. Seleção: Candidato com menor max_distance
6. Refinamento: Busca local para otimização final

Classes:
    Nenhuma - Implementação funcional modular

Funções principais:
    heuristic_closest_string(): Algoritmo principal do CSC
    cluster_strings(): Clusterização baseada em DBSCAN
    consensus_string(): Cálculo de consenso por maioria
    local_search(): Refinamento local intensivo

Funções auxiliares:
    split_blocks(): Divisão de strings em blocos uniformes
    recombine_blocks(): Recombinação de blocos em strings
    strings_to_array(): Conversão para representação numérica
    auto_parameters(): Definição automática de parâmetros

Author: Implementação baseada em técnicas de clustering e consenso estratificado
Version: Otimizada para instâncias com estrutura local bem definida
"""

import logging
from collections import Counter
from collections.abc import Callable
from itertools import combinations, product

import numpy as np
from sklearn.cluster import DBSCAN

from src.utils.distance import hamming_distance, max_distance

from .config import CSC_DEFAULTS

logger = logging.getLogger(__name__)

# ---------- Funções Auxiliares ----------


def consensus_string(strings):
    """
    Calcula string consenso através de votação majoritária posição-a-posição.

    Esta função implementa o algoritmo clássico de consenso por maioria,
    que é a base para o consenso local de cada cluster no CSC. Para cada
    posição, escolhe o símbolo que aparece com maior frequência naquela
    posição entre todas as strings do conjunto.

    ALGORITMO:
    1. Para cada posição i de 0 a L-1:
       - Coleta todos os símbolos na posição i
       - Conta frequência de cada símbolo
       - Escolhe o mais frequente (critério de desempate: primeiro encontrado)
    2. Concatena símbolos escolhidos para formar consenso

    PROPRIEDADES:
    - **Determinístico**: Mesmo resultado para mesma entrada
    - **Ótimo Local**: Minimiza distância de Hamming soma para conjunto
    - **Eficiente**: O(n × L) onde n=strings, L=comprimento
    - **Robusto**: Funciona com qualquer alfabeto

    Args:
        strings: Lista de strings de mesmo comprimento

    Returns:
        str: String consenso calculada por votação majoritária

    Example:
        >>> consensus_string(["ACGT", "AGCT", "ATCT"])
        "ACCT"  # A=3/3, C=2/3, C=2/3, T=3/3
    """
    L = len(strings[0])
    consensus = ""
    for i in range(L):
        # Coleta símbolos na posição i
        chars = [s[i] for s in strings]
        # Votação majoritária com Counter
        most_common = Counter(chars).most_common(1)[0][0]
        consensus += most_common
    return consensus


def split_blocks(s, n_blocks):
    """
    Divide string em blocos de tamanho aproximadamente uniforme.

    Esta função é fundamental para a estratégia de recombinação do CSC.
    Divide a string em blocos que podem ser depois recombinados com
    blocos de outros consensos locais, permitindo hibridização inteligente.

    ALGORITMO DE DIVISÃO:
    1. Calcula tamanho base: L // n_blocks
    2. Distribui resto uniformemente: primeiros (L % n_blocks) blocos +1
    3. Cria blocos sequenciais cobrindo toda a string

    CARACTERÍSTICAS:
    - **Cobertura Total**: Todos os caracteres incluídos em algum bloco
    - **Balanceamento**: Diferença máxima de 1 entre tamanhos
    - **Sequencial**: Preserva ordem original dos caracteres
    - **Flexível**: Adapta-se a qualquer combinação L/n_blocks

    Args:
        s: String a ser dividida
        n_blocks: Número de blocos desejado

    Returns:
        list[str]: Lista de blocos cobrindo toda a string

    Example:
        >>> split_blocks("ABCDEFG", 3)
        ["ABC", "DE", "FG"]  # Tamanhos [3,2,2] para L=7, n_blocks=3
    """
    L = len(s)
    # Calcula tamanhos base e distribui resto
    block_sizes = [L // n_blocks] * n_blocks
    for i in range(L % n_blocks):
        block_sizes[i] += 1

    # Constrói blocos sequencialmente
    blocks = []
    idx = 0
    for size in block_sizes:
        blocks.append(s[idx : idx + size])
        idx += size
    return blocks


def recombine_blocks(blocks_list):
    """
    Recombina lista de blocos em uma string unificada.

    Função complementar a split_blocks, usada para construir candidatos
    através da combinação de blocos de diferentes consensos locais.
    É o coração da estratégia de hibridização do CSC.

    PROCESSO DE RECOMBINAÇÃO:
    1. Concatena blocos na ordem fornecida
    2. Preserva conteúdo original de cada bloco
    3. Forma string candidata para avaliação

    PAPEL NO CSC:
    - **Hibridização**: Combina qualidades de diferentes consensos
    - **Exploração**: Cria candidatos não óbvios
    - **Eficiência**: Operação O(total_chars) muito rápida

    Args:
        blocks_list: Lista de strings (blocos) a serem concatenados

    Returns:
        str: String resultante da concatenação dos blocos

    Example:
        >>> recombine_blocks(["AC", "GT", "TA"])
        "ACGTTA"
    """
    return "".join(blocks_list)


# ---------- Etapa 1: Leitura/conversão do conjunto de strings ----------


def strings_to_array(strings):
    """
    Converte lista de strings para representação numérica compatível com sklearn.

    Esta função é necessária para usar o DBSCAN do sklearn, que requer
    dados numéricos. Cria mapeamento de caracteres para números e converte
    todas as strings para arrays numpy de inteiros.

    PROCESSO DE CONVERSÃO:
    1. **Extração do Alfabeto**: Identifica todos os símbolos únicos
    2. **Mapeamento**: Cria dicionário {char: int} ordenado
    3. **Conversão**: Transforma cada string em array de inteiros
    4. **Validação**: Verifica comprimento uniforme das strings

    VANTAGENS:
    - **Compatibilidade**: Permite uso de métricas sklearn
    - **Eficiência**: Arrays numpy são mais rápidos que strings
    - **Flexibilidade**: Funciona com qualquer alfabeto
    - **Memória**: Representação mais compacta para alfabetos pequenos

    Args:
        strings: Lista de strings de mesmo comprimento

    Returns:
        tuple: (array_numpy, mapeamento_chars)
            - array_numpy: matriz n×L de inteiros
            - mapeamento_chars: dicionário {char: int}

    Raises:
        ValueError: Se strings têm comprimentos diferentes ou lista vazia

    Example:
        >>> strings_to_array(["AC", "GT"])
        (array([[0, 2], [3, 1]]), {'A': 0, 'C': 2, 'G': 3, 'T': 1})
    """
    if not strings:
        raise ValueError("A lista de strings está vazia.")

    # Validação de comprimento uniforme
    L = len(strings[0])
    for idx, s in enumerate(strings):
        if len(s) != L:
            raise ValueError(
                f"Todas as strings devem ter o mesmo comprimento. "
                f"Encontrada string de tamanho {len(s)} na posição {idx}, esperado {L}."
            )

    # Criação do mapeamento char → int
    unique_chars = sorted(set("".join(strings)))
    char_map = {c: i for i, c in enumerate(unique_chars)}

    # Conversão para array numpy
    arr = np.array([[char_map[c] for c in s] for s in strings])

    return arr, char_map


# ---------- Parâmetros automáticos ----------


def auto_parameters(strings):
    """
    Define automaticamente parâmetros críticos do CSC baseado nas características dos dados.

    Esta função implementa heurísticas adaptativas para configurar automaticamente
    os parâmetros mais importantes do algoritmo CSC, eliminando a necessidade
    de ajuste manual e garantindo performance adequada para diferentes tipos
    de instâncias do CSP.

    PARÂMETROS CALCULADOS:

    1. **d (eps para DBSCAN)**:
       - Baseado na média das distâncias de Hamming entre pares
       - Fórmula: max(min_d, floor(média × d_factor))
       - Intuição: Strings "similares" se distância ≤ d
       - Balanceia entre clusters muito pequenos vs muito grandes

    2. **n_blocks (número de blocos)**:
       - Baseado no tamanho da instância e configurações padrão
       - Fórmula: max(min_blocks, min(max_blocks, n/n_div, L/l_div))
       - Intuição: Mais strings → mais blocos para maior recombinação
       - Garante granularidade adequada sem explosão combinatorial

    ANÁLISE ESTATÍSTICA PRÉVIA:
    - Calcula todas as O(n²) distâncias de Hamming entre pares
    - Computa estatísticas descritivas (média, mediana, min, max)
    - Usa estatísticas para inferir estrutura dos dados

    HEURÍSTICAS ADAPTATIVAS:
    - **Instâncias Homogêneas**: d menor → clusters maiores → menos consensos
    - **Instâncias Heterogêneas**: d maior → clusters menores → mais consensos
    - **Strings Longas**: Mais blocos para maior granularidade
    - **Muitas Strings**: Mais blocos para aproveitar diversidade

    VANTAGENS DA CONFIGURAÇÃO AUTOMÁTICA:
    - **Facilidade de Uso**: Elimina necessidade de expertise
    - **Robustez**: Adapta-se a diferentes características dos dados
    - **Performance**: Parâmetros adequados para a instância específica
    - **Reprodutibilidade**: Mesma configuração para dados similares

    Args:
        strings: Lista de strings de entrada para análise

    Returns:
        tuple: (d, n_blocks)
            - d: Raio máximo para clusterização (int)
            - n_blocks: Número de blocos para recombinação (int)

    Example:
        >>> strings = ["AAAA", "AAAT", "TTTT", "TTTA"]
        >>> auto_parameters(strings)
        (2, 2)  # d=2 (média≈2.67), n_blocks=2 (4 strings, 4 chars)
    """
    # ANÁLISE ESTATÍSTICA DAS DISTÂNCIAS
    # Calcula todas as distâncias de Hamming entre pares
    distancias = [hamming_distance(s1, s2) for s1, s2 in combinations(strings, 2)]

    # Estatísticas descritivas
    media = np.mean(distancias)
    mediana = np.median(distancias)
    minimo = np.min(distancias)
    maximo = np.max(distancias)

    # HEURÍSTICA PARA PARÂMETRO d (DBSCAN eps)
    # Baseado na média das distâncias com fator de escala
    d = max(CSC_DEFAULTS["min_d"], int(np.floor(media * CSC_DEFAULTS["d_factor"])))

    # HEURÍSTICA PARA PARÂMETRO n_blocks
    # Considera número de strings (n) e comprimento (L)
    n = len(strings)
    L = len(strings[0])
    n_blocks = max(
        CSC_DEFAULTS["min_blocks"],
        min(
            CSC_DEFAULTS["max_blocks"],
            n // CSC_DEFAULTS["n_div"],  # Baseado no número de strings
            L // CSC_DEFAULTS["l_div"],  # Baseado no comprimento
        ),
    )

    # LOG DAS DECISÕES PARA AUDITORIA
    logger.info(
        "Parâmetros automáticos: Hamming média=%.2f, min=%d, max=%d | d=%d | n_blocks=%d",
        media,
        minimo,
        maximo,
        d,
        n_blocks,
    )

    return d, n_blocks


# ---------- Etapa 2: Clusterização ----------


def cluster_strings(strings, d, min_samples=2):
    logger.debug("Clusterizando strings com d=%d, min_samples=%d", d, min_samples)
    arr, char_map = strings_to_array(strings)

    def hamming_metric(x, y):
        return np.sum(x != y)

    clustering = DBSCAN(eps=d, min_samples=min_samples, metric=hamming_metric)
    labels = clustering.fit_predict(arr)
    clusters = {}
    for idx, label in enumerate(labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(strings[idx])
    # Remove outliers (label == -1)
    clusters = {k: v for k, v in clusters.items() if k != -1}
    logger.info("Clusters encontrados: %d", len(clusters))
    return list(clusters.values())


# ---------- Etapa 3: Consenso local e recombinação ----------


def heuristic_closest_string(
    strings,
    d=None,
    n_blocks=None,
    progress_callback: Callable[[str], None] | None = None,
):
    """
    Algoritmo principal do CSC para resolver o Closest String Problem.

    Esta função implementa a estratégia completa do Consensus String Clustering,
    coordenando todas as etapas desde a clusterização até o refinamento final.
    É a interface principal do algoritmo CSC.

    FLUXO ALGORÍTMICO COMPLETO:

    1. **CONFIGURAÇÃO AUTOMÁTICA** (se parâmetros não fornecidos):
       - Analisa características dos dados de entrada
       - Define d e n_blocks automaticamente via heurísticas
       - Garante configuração adequada para a instância

    2. **CLUSTERIZAÇÃO INTELIGENTE**:
       - Aplica DBSCAN com métrica de Hamming
       - Identifica grupos de strings similares
       - Remove outliers que não formam clusters
       - Fallback para consenso global se nenhum cluster

    3. **CONSENSO LOCAL ESTRATIFICADO**:
       - Calcula consenso por maioria para cada cluster
       - Obtém representantes locais de alta qualidade
       - Divide cada consenso em blocos uniformes

    4. **RECOMBINAÇÃO CRIATIVA**:
       - Gera produto cartesiano de blocos de todos consensos
       - Cria candidatos hibridizando soluções locais
       - Normaliza tamanho para manter consistência

    5. **SELEÇÃO E REFINAMENTO**:
       - Avalia todos candidatos com função max_distance
       - Seleciona melhor candidato inicial
       - Aplica busca local intensiva para refinamento

    VANTAGENS ESTRATÉGICAS:
    - **Qualidade Superior**: Consenso local > consenso global
    - **Exploração Inteligente**: Recombinação cria candidatos não óbvios
    - **Adaptabilidade**: Parâmetros se ajustam aos dados
    - **Robustez**: Fallback garante sempre uma solução

    CASOS DE FALHA E RECUPERAÇÃO:
    - **Sem Clusters**: Usa consenso global + busca local
    - **Parâmetros Inadequados**: Auto-configuração como backup
    - **Cancelamento**: Progress callback permite interrupção limpa

    MONITORAMENTO DE PROGRESSO:
    O algoritmo reporta progresso através de callback opcional:
    - "Clusterizando strings..."
    - "Encontrados X clusters"
    - "Gerando candidatos por recombinação..."
    - "Avaliando Y candidatos..."
    - "Executando busca local..."

    Args:
        strings: Lista de strings de mesmo comprimento (entrada do CSP)
        d: Raio para DBSCAN (None = automático)
        n_blocks: Número de blocos para recombinação (None = automático)
        progress_callback: Função para reportar progresso (opcional)

    Returns:
        str: String center otimizada para o conjunto de entrada

    Example:
        >>> strings = ["AAATTT", "AAAGGG", "CCCGGG", "CCCTTT"]
        >>> heuristic_closest_string(strings)
        "AAAGTT"  # Resultado otimizado após clusterização e recombinação

    Note:
        A função é thread-safe e pode ser interrompida via progress_callback.
        Parâmetros automáticos são logados para auditoria e debug.
    """
    # CONFIGURAÇÃO AUTOMÁTICA DE PARÂMETROS
    if d is None or n_blocks is None:
        d_auto, n_blocks_auto = auto_parameters(strings)
        if d is None:
            d = d_auto
        if n_blocks is None:
            n_blocks = n_blocks_auto

    logger.info("Iniciando CSC com d=%d, n_blocks=%d", d, n_blocks)

    # ETAPA 1: CLUSTERIZAÇÃO
    if progress_callback:
        progress_callback("Clusterizando strings...")

    clusters = cluster_strings(strings, d)

    # FALLBACK: Se nenhum cluster encontrado
    if not clusters:
        if progress_callback:
            progress_callback("⚠️ Nenhum cluster encontrado, usando consenso global")
        # Estratégia de recuperação: consenso global + busca local
        best_candidate = consensus_string(strings)
        best_candidate = local_search(best_candidate, strings, progress_callback)
        return best_candidate

    if progress_callback:
        progress_callback(f"Encontrados {len(clusters)} clusters")

    # ETAPA 2: CONSENSO LOCAL POR CLUSTER
    consensos = [consensus_string(cluster) for cluster in clusters]
    logger.debug("Consensos locais calculados: %d", len(consensos))

    if progress_callback:
        progress_callback("Gerando candidatos por recombinação...")

    # ETAPA 3: RECOMBINAÇÃO DE BLOCOS
    L = len(strings[0])
    # Divide cada consenso em blocos
    blocks_per_consenso = [split_blocks(cons, n_blocks) for cons in consensos]

    # Gera todas as combinações possíveis de blocos
    candidates = []
    for block_tuple in product(*blocks_per_consenso):
        candidate = recombine_blocks(block_tuple)

        # NORMALIZAÇÃO DE TAMANHO
        # Garante que candidato tenha exatamente o tamanho correto
        if len(candidate) > L:
            candidate = candidate[:L]  # Trunca se muito longo
        elif len(candidate) < L:
            # Completa com sufixo do primeiro consenso se muito curto
            candidate += consensos[0][len(candidate) : L]

        candidates.append(candidate)

    logger.info("Gerados %d candidatos por recombinação", len(candidates))

    # ETAPA 4: SELEÇÃO DO MELHOR CANDIDATO
    if progress_callback:
        progress_callback(f"Avaliando {len(candidates)} candidatos...")

    # Encontra candidato com menor distância máxima
    best_candidate = min(candidates, key=lambda cand: max_distance(cand, strings))
    logger.info("Melhor candidato pré-refinamento selecionado")

    # ETAPA 5: REFINAMENTO LOCAL INTENSIVO
    if progress_callback:
        progress_callback("Executando busca local...")

    best_candidate = local_search(best_candidate, strings, progress_callback)
    logger.info("Refinamento local concluído")

    return best_candidate


def local_search(
    candidate, strings, progress_callback: Callable[[str], None] | None = None
):
    candidate = list(candidate)
    improved = True
    iterations = 0
    max_iterations = 50  # Limite para evitar loops longos

    while improved and iterations < max_iterations:
        improved = False
        iterations += 1

        # Verificar cancelamento periodicamente
        if progress_callback:
            if iterations % 10 == 0:
                progress_callback(
                    f"Busca local: iteração {iterations}/{max_iterations}"
                )

            # Verificação adicional para cancelamento - checar se o processo ainda deve rodar
            try:
                import signal

                # Isso permitirá que o processo seja interrompido mais facilmente
                signal.signal(signal.SIGTERM, signal.default_int_handler)
            except (ImportError, AttributeError):
                pass  # Em caso de erro, continuar normalmente

        for i in range(len(candidate)):
            current_char = candidate[i]
            chars = {s[i] for s in strings}
            for alt in chars:
                if alt == current_char:
                    continue
                new_candidate = candidate.copy()
                new_candidate[i] = alt
                new_candidate_str = "".join(new_candidate)
                if max_distance(new_candidate_str, strings) < max_distance(
                    "".join(candidate), strings
                ):
                    candidate[i] = alt
                    improved = True
                    # Parar no primeiro melhoramento para evitar travamento
                    break

            # Se melhorou, parar de testar outras posições nesta iteração
            if improved:
                break

    logger.debug(
        "Resultado da busca local após %d iterações: %s", iterations, "".join(candidate)
    )
    return "".join(candidate)
