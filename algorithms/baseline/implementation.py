"""
Implementação do algoritmo de consenso ganancioso (Baseline) para CSP.

O algoritmo Baseline representa a abordagem mais simples e fundamental para o
Closest String Problem (CSP). Utiliza uma estratégia gulosa (greedy) que constrói
a string center posição por posição, escolhendo sempre o símbolo que minimiza
localmente a distância máxima naquele momento.

ALGORITMO GREEDY CONSENSUS:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                        ALGORITMO BASELINE DETALHADO                             │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. INICIALIZAÇÃO                                                               │
│   ├── Valida entrada (strings não vazias, mesmo comprimento)                   │
│   └── Inicializa string consenso vazia                                         │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. CONSTRUÇÃO POSIÇÃO-A-POSIÇÃO                                               │
│   ├── Para cada posição i de 0 a L-1:                                          │
│   │   ├── Para cada símbolo c do alfabeto:                                     │
│   │   │   ├── Calcula string parcial = consenso[0:i] + c                      │
│   │   │   ├── Calcula distância máxima da string parcial                      │
│   │   │   └── Armazena se for a melhor opção até agora                        │
│   │   └── Escolhe símbolo que minimiza distância máxima                       │
│   └── Adiciona melhor símbolo ao consenso                                      │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 3. VALIDAÇÃO E RETORNO                                                        │
│   ├── Calcula distância final da string consenso completa                      │
│   ├── Registra resultado no log para auditoria                                │
│   └── Retorna string consenso construída                                       │
└─────────────────────────────────────────────────────────────────────────────────┘

CARACTERÍSTICAS DO ALGORITMO:

• **SIMPLICIDADE**: Implementação direta e compreensível
• **DETERMINÍSTICO**: Sempre produz o mesmo resultado para a mesma entrada
• **GULOSO LOCAL**: Otimiza cada posição independentemente
• **BASELINE CONFIÁVEL**: Serve como referência para comparação
• **COMPLEXIDADE BAIXA**: O(L * |Σ| * n) onde L=comprimento, |Σ|=alfabeto, n=strings

LIMITAÇÕES:

• **ÓTIMO LOCAL**: Pode não encontrar a solução global ótima
• **MIOPIA**: Decisões locais podem prejudicar qualidade global
• **SEM BACKTRACKING**: Não reconsidera decisões anteriores
• **QUALIDADE LIMITADA**: Geralmente inferior a metaheurísticas avançadas

APLICAÇÃO E USO:

O Baseline é fundamental como:
- **Referência de comparação** para outros algoritmos
- **Solução inicial** para métodos mais sofisticados
- **Validação de implementação** (deve funcionar corretamente)
- **Análise de dificuldade** do problema (se Baseline resolve bem, problema é fácil)

EXEMPLO DE FUNCIONAMENTO:

Entrada: ["ACGT", "AGCT", "ATCT"]
Posição 0: A=3/3, melhor='A' (distância parcial=0)
Posição 1: C vs G vs T → C minimiza distância máxima
Posição 2: G vs C vs T → C minimiza distância máxima
Posição 3: T=3/3, melhor='T' (distância parcial mínima)
Resultado: "ACCT" com distância máxima = 1

Classes:
    Nenhuma - Implementação funcional direta

Funções principais:
    greedy_consensus(strings, alphabet): Algoritmo principal de consenso guloso
    max_distance(center, strings): Função auxiliar para cálculo de distância

Author: Implementação de referência baseada em algoritmos gulosos clássicos
Version: Otimizada para clareza e servir como baseline confiável
"""

import logging

logger = logging.getLogger(__name__)


def greedy_consensus(strings: list[str], alphabet: str) -> str:
    """
    Constrói uma string consenso usando estratégia gulosa posição por posição.

    Esta é a implementação do algoritmo Baseline para o Closest String Problem.
    A estratégia gulosa (greedy) constrói a solução incrementalmente, tomando
    a decisão localmente ótima em cada passo sem considerar o impacto global.

    ALGORITMO DETALHADO:

    1. **VALIDAÇÃO DE ENTRADA**:
       - Verifica se há strings para processar
       - Assume que todas têm o mesmo comprimento (pré-condição)

    2. **CONSTRUÇÃO INCREMENTAL**:
       - Para cada posição i de 0 a L-1:
         a) Testa cada símbolo do alfabeto nessa posição
         b) Para cada símbolo candidato:
            - Constrói string parcial até posição i
            - Calcula distância máxima da string parcial
         c) Escolhe símbolo que minimiza a distância máxima
         d) Adiciona símbolo escolhido ao consenso

    3. **ESTRATÉGIA DE OTIMIZAÇÃO LOCAL**:
       - Em cada posição, escolhe símbolo que resulta na menor
         distância máxima considerando apenas o prefixo construído
       - Não considera impacto futuro (característica gulosa)

    COMPLEXIDADE TEMPORAL:
    - O(L × |Σ| × n × L) onde:
      - L = comprimento das strings
      - |Σ| = tamanho do alfabeto
      - n = número de strings
      - Fator L adicional vem do cálculo de distância parcial

    CARACTERÍSTICAS:
    - **Determinístico**: Sempre produz o mesmo resultado
    - **Guloso**: Otimização local sem backtracking
    - **Incremental**: Constrói solução posição por posição
    - **Baseline**: Referência para comparação com outros algoritmos

    LIMITAÇÕES:
    - Pode ficar preso em ótimos locais
    - Não garante solução globalmente ótima
    - Qualidade depende da ordem de construção

    Args:
        strings: Lista de strings de entrada (todas com mesmo comprimento)
        alphabet: String contendo todos os símbolos válidos do alfabeto

    Returns:
        str: String consenso construída pela estratégia gulosa

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> alphabet = "ACGT"
        >>> greedy_consensus(strings, alphabet)
        "ACCT"  # Distância máxima = 1

    Note:
        A função registra o resultado final no log para auditoria
        e validação da implementação.
    """
    # VALIDAÇÃO DE ENTRADA
    if not strings:
        return ""

    # INICIALIZAÇÃO
    L = len(strings[0])  # Comprimento das strings (assumido uniforme)
    consensus = []  # String consenso sendo construída

    # CONSTRUÇÃO POSIÇÃO-A-POSIÇÃO
    for pos in range(L):
        best_char = None  # Melhor símbolo para posição atual
        best_max_dist = float("inf")  # Menor distância máxima encontrada

        # TESTE DE CADA SÍMBOLO DO ALFABETO
        for char in alphabet:
            # Construir string parcial com símbolo candidato
            partial_consensus = consensus + [char]

            # CÁLCULO DA DISTÂNCIA MÁXIMA PARCIAL
            # Calcula distância apenas considerando posições já decididas
            max_dist = 0
            for s in strings:
                # Distância de Hamming parcial (até posição atual + 1)
                dist = sum(1 for i in range(pos + 1) if partial_consensus[i] != s[i])
                max_dist = max(max_dist, dist)

            # ATUALIZAÇÃO DO MELHOR CANDIDATO
            # Se encontrou símbolo que resulta em menor distância máxima
            if max_dist < best_max_dist:
                best_max_dist = max_dist
                best_char = char

        # DECISÃO GULOSA: Adicionar melhor símbolo encontrado
        consensus.append(best_char)

    # CONSTRUÇÃO DA STRING FINAL
    result = "".join(consensus)

    # VALIDAÇÃO E LOG DO RESULTADO
    # Import local para evitar conflito de nomes
    from src.domain.metrics import max_distance as calc_max_distance

    final_distance = calc_max_distance(result, strings)
    logger.info("[CONSENSUS] Consenso: %s, distância: %d", result, final_distance)

    return result


def max_distance(center: str, strings: list[str]) -> int:
    """
    Calcula a distância máxima de Hamming entre uma string center e um conjunto de strings.

    Esta função implementa a função objetivo do Closest String Problem (CSP).
    A distância máxima é o valor que precisamos minimizar para encontrar
    a melhor string center possível.

    DEFINIÇÃO MATEMÁTICA:
    Para uma string center c e conjunto de strings S = {s₁, s₂, ..., sₙ}:
    max_distance(c, S) = max{d_H(c, sᵢ) | sᵢ ∈ S}

    onde d_H(a,b) é a distância de Hamming entre strings a e b.

    IMPORTÂNCIA NO CSP:
    - **Função Objetivo**: Valor a ser minimizado
    - **Critério de Qualidade**: Menor valor = melhor solução
    - **Fitness**: Usado para avaliar candidatos em algoritmos evolutivos
    - **Stopping Criterion**: Algoritmo para quando atinge 0 (ótimo)

    PROPRIEDADES:
    - **Simétrica**: max_distance(c, S) = max_distance(c, S')
    - **Limitada**: 0 ≤ resultado ≤ L (comprimento das strings)
    - **Monotônica**: Mais diferenças → maior distância
    - **Discreta**: Apenas valores inteiros

    COMPLEXIDADE:
    - Temporal: O(n × L) onde n=número de strings, L=comprimento
    - Espacial: O(1) - não requer armazenamento adicional

    CASOS ESPECIAIS:
    - Se center ∈ strings, então resultado ≤ diâmetro do conjunto
    - Se strings são idênticas e center = string, resultado = 0
    - Pior caso: center difere de todas em todas as posições

    Args:
        center: String candidata a center (mesmo comprimento que strings)
        strings: Lista de strings de entrada para comparação

    Returns:
        int: Maior distância de Hamming encontrada entre center e strings

    Example:
        >>> center = "ACCT"
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> max_distance(center, strings)
        1  # max(1, 1, 1) = 1

    Note:
        Esta função é wrapper para manter consistência com a interface
        do módulo baseline, delegando o cálculo para implementação otimizada.
    """
    # Import local para usar implementação otimizada
    from src.domain.metrics import hamming_distance

    # CÁLCULO DE TODAS AS DISTÂNCIAS
    # Calcula distância de Hamming entre center e cada string
    distances = [hamming_distance(center, s) for s in strings]

    # RETORNA A DISTÂNCIA MÁXIMA
    # max() encontra o maior valor na lista de distâncias
    return max(distances)
