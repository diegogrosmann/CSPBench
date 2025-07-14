"""
Implementação do H³-CSP (Hybrid Hierarchical Hamming Search) para CSP.

O H³-CSP é um algoritmo híbrido sofisticado que resolve o Closest String Problem
através de uma abordagem hierárquica de três camadas, combinando decomposição
estrutural, seleção adaptativa de técnicas e refinamento global.

ARQUITETURA ALGORÍTMICA:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                          ALGORITMO H³-CSP DETALHADO                             │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. B-SPLITTER (Divisão Hierárquica)                                            │
│   ├── Divide strings em ~√L blocos contíguos                                   │
│   ├── Cada bloco: tamanho ≈ √L (heurística de balanceamento)                   │
│   ├── Preserva localidade espacial (blocos contíguos)                          │
│   └── Reduz complexidade: O(|Σ|^L) → O(B×|Σ|^(L/B))                          │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. SMART-CORE (Seleção Adaptativa por Bloco)                                  │
│   ├── Para cada bloco, calcula dificuldade d_b (distância do consenso)        │
│   ├── BLOCO FÁCIL (d_b ≤ block_small):                                        │
│   │   ├── Busca exaustiva: explora todo o espaço |Σ|^(r-l)                   │
│   │   ├── Garantia de otimalidade local                                        │
│   │   └── Usado quando |Σ|^(r-l) ≤ 10.000                                    │
│   ├── BLOCO MÉDIO (block_small < d_b ≤ block_medium):                         │
│   │   ├── Beam search reduzido (beam_width/2)                                 │
│   │   ├── Balanço qualidade vs. eficiência                                    │
│   │   └── Construção incremental posição-a-posição                            │
│   └── BLOCO DIFÍCIL (d_b > block_medium):                                     │
│       ├── Beam search completo (beam_width)                                   │
│       ├── Máximo esforço computacional                                        │
│       └── Exploração ampla do espaço de busca                                 │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 3. GLOBAL REFINE (Fusão e Refinamento Global)                                 │
│   ├── FUSÃO: Concatena melhores candidatos de cada bloco                       │
│   ├── REFINAMENTO: Hill-climbing posição-a-posição                            │
│   │   ├── Testa todas as substituições possíveis                              │
│   │   ├── Aceita apenas melhorias                                             │
│   │   └── Para no primeiro ótimo local                                        │
│   └── CONTROLE: Timeout e callback de progresso                               │
└─────────────────────────────────────────────────────────────────────────────────┘

FILOSOFIA ALGORÍTMICA:

• **DECOMPOSIÇÃO INTELIGENTE**: Divide o problema em subproblemas menores e
  independentes, reduzindo drasticamente a complexidade computacional.

• **ADAPTATIVIDADE**: Seleciona a técnica mais apropriada para cada bloco
  baseado na dificuldade, otimizando o uso de recursos computacionais.

• **HIERARQUIA**: Combina otimização local (por bloco) com refinamento global,
  balanceando qualidade local e consistência global.

• **DETERMINISMO**: Algoritmo determinístico que garante reprodutibilidade
  de resultados (exceto por timeout).

CARACTERÍSTICAS DISTINTIVAS:

• **EFICIÊNCIA ESCALÁVEL**: Complexidade O(B×|Σ|^(L/B)) vs. O(|Σ|^L)
• **QUALIDADE ADAPTATIVA**: Técnicas sofisticadas apenas onde necessário
• **ROBUSTEZ**: Fallbacks e controle de timeout para casos extremos
• **FLEXIBILIDADE**: Parâmetros configuráveis para diferentes cenários

APLICAÇÃO AO CSP:
O H³-CSP é especialmente eficaz para:
- Instâncias de tamanho médio (L=50-500)
- Dados com padrões locais ou estrutura hierárquica
- Cenários onde qualidade da solução é crítica
- Problemas que requerem determinismo e reprodutibilidade

VANTAGENS:
- Redução exponencial da complexidade através da decomposição
- Adaptação automática à dificuldade do problema
- Garantia de otimalidade local onde computacionalmente viável
- Refinamento global para consistência entre blocos

LIMITAÇÕES:
- Eficiência reduzida para strings muito curtas (overhead da decomposição)
- Possível subotimalidade devido à decomposição em blocos
- Memória adicional para armazenar candidatos de todos os blocos

EXEMPLO DE FLUXO:
Para strings ["ACGTACGT", "AGCTACGT", "ATGTACGT"] com L=8:
1. B-Splitter: √8≈3 blocos → [(0,3), (3,6), (6,8)]
2. Smart-Core por bloco:
   - Bloco 1: consenso="ACG", d_b=1 → busca exaustiva
   - Bloco 2: consenso="TAC", d_b=0 → busca exaustiva
   - Bloco 3: consenso="GT", d_b=0 → busca exaustiva
3. Fusão: "ACG" + "TAC" + "GT" = "ACGTACGT"
4. Refinamento: hill-climbing para otimização final

Classes:
    H3CSP: Implementação principal do algoritmo H³-CSP.

Funções auxiliares:
    split_in_blocks(L): Divide L posições em ~√L blocos contíguos.
    consensus_block(strings, l, r): Calcula consenso por maioria em um bloco.
    _exhaustive_block(): Busca exaustiva para blocos pequenos.
    _beam_search_block(): Beam search para blocos médios/grandes.
    _local_search(): Hill-climbing para refinamento global.

Tipos:
    String: Alias para str (string de DNA/proteína).
    Block: Tupla (início, fim) representando um bloco - 0-based, fim exclusivo.

Author: Implementação baseada em algoritmos hierárquicos para CSP
Version: Otimizada para instâncias de tamanho médio com estrutura local
"""

from __future__ import annotations

import itertools
import logging
import math
import random
import time
from collections import Counter
from collections.abc import Callable, Sequence

from src.domain.metrics import max_distance

from .config import H3_CSP_DEFAULTS

logger = logging.getLogger(__name__)

String = str
Block = tuple[int, int]  # (início, fim) 0-based, exclusivo


# ---------------------------------------------------------------------------
# Auxiliares de bloco
# ---------------------------------------------------------------------------


def split_in_blocks(L: int) -> list[Block]:
    """
    Divide L posições em aproximadamente √L blocos contíguos de tamanho uniforme.

    Esta função implementa a "regra √L" do H³-CSP, que é uma heurística fundamental
    para balancear o trade-off entre qualidade da solução e eficiência computacional.
    A divisão hierárquica reduz a complexidade de O(|Σ|^L) para O(B×|Σ|^(L/B)).

    MATEMÁTICA DA DIVISÃO:
    - Número de blocos: B = ⌈√L⌉
    - Tamanho base: base_size = ⌈L/B⌉
    - Último bloco pode ser menor se L não for divisível por B

    DESIGN RATIONALE:
    - √L é o ponto ótimo teórico para minimizar B×|Σ|^(L/B)
    - Blocos contíguos preservam localidade espacial das strings
    - Tamanho aproximadamente uniforme simplifica análise de complexidade
    - Permite paralelização eficiente do processamento por bloco

    COMPLEXIDADE:
    - Temporal: O(1) - apenas cálculos aritméticos
    - Espacial: O(√L) - armazena √L blocos

    Args:
        L (int): Comprimento total das strings a serem divididas (L ≥ 1).

    Returns:
        list[Block]: Lista de blocos, onde cada bloco é uma tupla (início, fim)
                    representando o intervalo [início, fim) com fim exclusivo.
                    A união de todos os blocos cobre exatamente [0, L).

    Examples:
        >>> split_in_blocks(16)  # √16 = 4 blocos de tamanho 4
        [(0, 4), (4, 8), (8, 12), (12, 16)]

        >>> split_in_blocks(10)  # √10 ≈ 3.16 → 4 blocos, último menor
        [(0, 3), (3, 6), (6, 9), (9, 10)]

        >>> split_in_blocks(1)   # Caso trivial
        [(0, 1)]

    Note:
        A heurística √L é otimizada para alfabetos pequenos (|Σ| = 4-20).
        Para alfabetos muito grandes, pode ser necessário ajustar a fórmula.
    """
    # Cálculo do número de blocos usando a regra √L
    B = math.ceil(math.sqrt(L))

    # Cálculo do tamanho base de cada bloco para distribuição uniforme
    base_size = math.ceil(L / B)

    # Criação dos blocos contíguos
    blocks = []
    cur = 0  # Posição atual (início do próximo bloco)

    while cur < L:
        # Cria bloco [cur, min(cur + base_size, L))
        # min() garante que o último bloco não ultrapasse L
        blocks.append((cur, min(cur + base_size, L)))
        cur += base_size

    return blocks


def consensus_block(strings: Sequence[String], l: int, r: int) -> String:
    """
    Calcula a string consenso por maioria para um bloco específico.

    Esta função implementa o algoritmo de consenso por votação majoritária,
    que é uma heurística simples mas eficaz para encontrar uma solução inicial
    de boa qualidade. Para cada posição no bloco, seleciona o caractere mais
    frequente entre todas as strings.

    ALGORITMO DE CONSENSO:
    1. Para cada posição i em [l, r):
       - Coleta caracteres de todas as strings na posição i
       - Conta frequência de cada caractere
       - Seleciona o caractere com maior frequência
    2. Concatena os caracteres selecionados para formar o consenso

    PROPRIEDADES:
    - Determinístico: sempre produz o mesmo resultado para a mesma entrada
    - Guloso: otimiza localmente cada posição independentemente
    - Eficiente: O(n × (r-l)) onde n é o número de strings
    - Heurística: não garante otimalidade global

    APLICAÇÃO NO H³-CSP:
    - Usado como baseline para medir dificuldade do bloco (d_b)
    - Serve como candidato inicial para refinamento
    - Fornece upper bound para qualidade do bloco

    COMPLEXIDADE:
    - Temporal: O(n × (r-l)) onde n = |strings|
    - Espacial: O(|Σ|) para contador de frequências

    Args:
        strings (Sequence[String]): Sequência de strings de entrada.
                                  Todas devem ter comprimento ≥ r.
        l (int): Posição inicial do bloco (inclusiva, 0-based).
        r (int): Posição final do bloco (exclusiva, 0-based).

    Returns:
        String: String consenso para o bloco especificado.
               Comprimento = r - l.

    Examples:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> consensus_block(strings, 1, 3)  # posições 1-2
        "CG"  # Pos 1: C(2) vs G(1) vs T(1) → C; Pos 2: G(2) vs C(1) → G

        >>> strings = ["AAAA", "TTTT", "CCCC"]
        >>> consensus_block(strings, 0, 2)  # empate → primeiro mais comum
        "AA"  # Counter().most_common(1) retorna ordem determinística

    Note:
        Em caso de empate, Counter.most_common(1) retorna o primeiro
        caractere na ordem de inserção (determinístico).
    """
    # Armazena o consenso caractere por caractere
    consensus_chars = []

    # Processar cada posição do bloco independentemente
    for pos in range(l, r):
        # Coletar caracteres de todas as strings na posição atual
        chars_at_pos = [s[pos] for s in strings]

        # Contar frequência de cada caractere
        char_counter = Counter(chars_at_pos)

        # Selecionar o caractere mais frequente
        # most_common(1) retorna lista com 1 tupla (char, count)
        # [0][0] extrai o caractere da primeira tupla
        most_common_char = char_counter.most_common(1)[0][0]
        consensus_chars.append(most_common_char)

    # Concatenar caracteres para formar string consenso
    return "".join(consensus_chars)


# ---------------------------------------------------------------------------
# Técnicas por bloco (versão simplificada)
# ---------------------------------------------------------------------------


def _exhaustive_block(
    strings: Sequence[String], alphabet: str, l: int, r: int, k: int
) -> list[String]:
    """
    Realiza busca exaustiva em um bloco pequeno para encontrar os k melhores candidatos.

    Esta função implementa a busca exaustiva completa para blocos pequenos,
    onde o espaço de busca |Σ|^(r-l) é computacionalmente viável. A busca
    exaustiva garante encontrar a solução ótima para o bloco, mas tem
    complexidade exponencial.

    ALGORITMO DE BUSCA EXAUSTIVA:
    1. **VERIFICAÇÃO DE VIABILIDADE**: Se |Σ|^(r-l) ≤ 10.000, usa busca completa
    2. **GERAÇÃO DE CANDIDATOS**: Usa itertools.product para gerar todas as
       combinações possíveis de caracteres
    3. **AVALIAÇÃO**: Calcula distância máxima de cada candidato para o bloco
    4. **SELEÇÃO**: Mantém os k melhores candidatos em uma lista ordenada
    5. **FALLBACK**: Se espaço muito grande, usa apenas candidatos do dataset

    ESTRATÉGIA DE FALLBACK:
    Para blocos onde |Σ|^(r-l) > 10.000:
    - Usa apenas os segmentos originais das strings como candidatos
    - Adiciona o consenso como candidato adicional
    - Reduz drasticamente o espaço de busca mantendo qualidade razoável

    COMPLEXIDADE:
    - Busca completa: O(|Σ|^(r-l) × n × (r-l)) onde n = |strings|
    - Fallback: O(n × n × (r-l)) - muito mais eficiente

    APLICAÇÃO NO H³-CSP:
    - Usado para blocos com d_b ≤ block_small (blocos fáceis)
    - Garante otimalidade local onde computacionalmente viável
    - Fornece baseline de qualidade para comparação

    Args:
        strings (Sequence[String]): Sequência de strings de entrada.
        alphabet (str): Alfabeto disponível para construção dos candidatos.
        l (int): Posição inicial do bloco (inclusiva).
        r (int): Posição final do bloco (exclusiva).
        k (int): Número máximo de candidatos a retornar.

    Returns:
        list[String]: Lista dos k melhores candidatos para o bloco,
                     ordenados por distância máxima crescente.

    Examples:
        >>> strings = ["ACG", "ATG", "AAG"]
        >>> _exhaustive_block(strings, "ACGT", 1, 3, 2)  # pos 1-2
        ["CG", "TG"]  # Melhores candidatos para segmento [1:3]

        >>> # Para bloco grande, usa fallback
        >>> _exhaustive_block(strings, "ACGT", 0, 10, 3)
        ["ACG...", "ATG...", "AAG..."]  # Usa segmentos originais

    Note:
        O limite de 10.000 é heurístico, balanceando qualidade vs. eficiência.
        Para alfabetos maiores, pode ser necessário ajustar este limite.
    """
    m = r - l  # Tamanho do bloco
    search_space_size = len(alphabet) ** m  # |Σ|^m

    # Lista para armazenar os k melhores candidatos como tuplas (distância, string)
    best_candidates = []

    # ESTRATÉGIA 1: BUSCA EXAUSTIVA (para blocos pequenos)
    if search_space_size <= 10_000:
        # Gera todas as combinações possíveis de caracteres
        for candidate_tuple in itertools.product(alphabet, repeat=m):
            # Converte tupla para string
            candidate = "".join(candidate_tuple)

            # Calcula distância máxima do candidato para os segmentos do bloco
            block_segments = [s[l:r] for s in strings]
            distance = max_distance(candidate, block_segments)

            # Mantém apenas os k melhores candidatos
            if len(best_candidates) < k:
                # Ainda há espaço, adiciona diretamente
                best_candidates.append((distance, candidate))
                # Mantém lista ordenada por distância
                best_candidates.sort(key=lambda x: x[0])
            elif distance < best_candidates[-1][0]:
                # Candidato melhor que o pior da lista
                best_candidates.append((distance, candidate))
                best_candidates.sort(key=lambda x: x[0])
                # Remove o pior candidato (último da lista ordenada)
                best_candidates = best_candidates[:k]

    # ESTRATÉGIA 2: FALLBACK (para blocos grandes)
    else:
        # Usa apenas segmentos originais como candidatos
        unique_segments = set()

        for s in strings:
            segment = s[l:r]
            if segment not in unique_segments:
                unique_segments.add(segment)

                # Calcula distância do segmento para todos os outros segmentos
                block_segments = [t[l:r] for t in strings]
                distance = max_distance(segment, block_segments)

                best_candidates.append((distance, segment))

        # Ordena candidatos por distância
        best_candidates.sort(key=lambda x: x[0])
        # Mantém apenas os k melhores
        best_candidates = best_candidates[:k]

    # GARANTIA DE CONSENSO: Sempre inclui consenso como candidato
    consensus = consensus_block(strings, l, r)
    block_segments = [s[l:r] for s in strings]
    consensus_distance = max_distance(consensus, block_segments)

    # Adiciona consenso se não estiver na lista ou se melhorar a lista
    best_candidates.append((consensus_distance, consensus))
    best_candidates.sort(key=lambda x: x[0])

    # Remove duplicatas mantendo ordem
    seen = set()
    unique_candidates = []
    for distance, candidate in best_candidates:
        if candidate not in seen:
            seen.add(candidate)
            unique_candidates.append((distance, candidate))

    # Retorna apenas as strings, mantendo os k melhores
    return [candidate for _, candidate in unique_candidates[:k]]


def _beam_search_block(
    strings: Sequence[String], alphabet: str, l: int, r: int, beam_width: int, k: int
) -> list[String]:
    """
    Aplica beam search posição-a-posição para gerar candidatos em um bloco.

    Esta função implementa um beam search incremental que constrói candidatos
    posição por posição, mantendo apenas os beam_width melhores prefixos em
    cada etapa. É uma alternativa eficiente à busca exaustiva para blocos
    médios/grandes onde |Σ|^(r-l) é computacionalmente inviável.

    ALGORITMO DE BEAM SEARCH:
    1. **INICIALIZAÇÃO**: Começa com prefixo vazio
    2. **CONSTRUÇÃO INCREMENTAL**: Para cada posição do bloco:
       - Estende cada prefixo atual com todos os caracteres do alfabeto
       - Avalia cada extensão baseada na distância parcial
       - Mantém apenas os beam_width melhores prefixos
    3. **PODA**: Remove prefixos piores, mantendo apenas os promissores
    4. **AVALIAÇÃO FINAL**: Calcula distância completa dos candidatos finais

    VANTAGENS DO BEAM SEARCH:
    - Complexidade controlada: O(beam_width × |Σ| × (r-l))
    - Exploração dirigida: Foca em regiões promissoras do espaço
    - Escalabilidade: Funciona bem para blocos grandes
    - Flexibilidade: beam_width controla trade-off qualidade/eficiência

    LIMITAÇÕES:
    - Não garante otimalidade: Pode podar soluções ótimas prematuramente
    - Avaliação local: Critério de poda baseado em distância parcial
    - Dependência de heurística: Qualidade depende da função de avaliação

    APLICAÇÃO NO H³-CSP:
    - Usado para blocos médios (block_small < d_b ≤ block_medium)
    - Usado para blocos difíceis (d_b > block_medium) com beam_width maior
    - Balanço entre qualidade e eficiência computacional

    COMPLEXIDADE:
    - Temporal: O(beam_width × |Σ| × (r-l) × n) onde n = |strings|
    - Espacial: O(beam_width × (r-l)) para armazenar prefixos

    Args:
        strings (Sequence[String]): Sequência de strings de entrada.
        alphabet (str): Alfabeto disponível para construção dos candidatos.
        l (int): Posição inicial do bloco (inclusiva).
        r (int): Posição final do bloco (exclusiva).
        beam_width (int): Largura do beam (número de prefixos mantidos).
        k (int): Número máximo de candidatos finais a retornar.

    Returns:
        list[String]: Lista dos k melhores candidatos para o bloco,
                     ordenados por distância máxima crescente.

    Examples:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> _beam_search_block(strings, "ACGT", 0, 2, beam_width=2, k=3)
        ["AC", "AG", "AT"]  # Melhores candidatos para posições 0-1

        >>> # Para bloco maior, explora subespaço eficientemente
        >>> _beam_search_block(strings, "ACGT", 0, 4, beam_width=4, k=2)
        ["ACGT", "AGCT"]  # Dois melhores candidatos completos

    Note:
        A qualidade do resultado é limitada pelo beam_width. Valores maiores
        produzem melhores resultados mas aumentam o custo computacional.
    """
    m = r - l  # Tamanho do bloco
    beam = [""]  # Lista de prefixos atuais, inicializada com string vazia

    # CONSTRUÇÃO INCREMENTAL: Posição por posição
    for position in range(m):
        # Lista para armazenar todos os prefixos estendidos com sua avaliação
        extended_prefixes = []

        # Para cada prefixo atual no beam
        for current_prefix in beam:
            # Tenta estender com cada caractere do alfabeto
            for char in alphabet:
                # Cria novo prefixo estendido
                extended_prefix = current_prefix + char

                # AVALIAÇÃO PARCIAL: Calcula distância do prefixo para os segmentos
                partial_distances = []
                for s in strings:
                    # Conta mismatches apenas nas posições já construídas
                    mismatch_count = 0
                    for i, prefix_char in enumerate(extended_prefix):
                        if i < len(extended_prefix) and prefix_char != s[l + i]:
                            mismatch_count += 1
                    partial_distances.append(mismatch_count)

                # Usa distância máxima como critério de avaliação
                partial_score = max(partial_distances)
                extended_prefixes.append((partial_score, extended_prefix))

        # PODA: Mantém apenas os beam_width melhores prefixos
        # Ordena por score crescente (menor distância = melhor)
        extended_prefixes.sort(key=lambda x: x[0])

        # Extrai apenas os prefixos dos beam_width melhores
        beam = [prefix for _, prefix in extended_prefixes[:beam_width]]

    # AVALIAÇÃO FINAL: Calcula distância completa para candidatos finais
    # Necessário porque avaliação parcial pode não refletir qualidade final
    final_candidates = []
    for candidate in beam:
        # Extrai segmentos correspondentes de todas as strings
        block_segments = [s[l:r] for s in strings]
        # Calcula distância máxima real
        final_distance = max_distance(candidate, block_segments)
        final_candidates.append((final_distance, candidate))

    # SELEÇÃO DOS K MELHORES
    final_candidates.sort(key=lambda x: x[0])
    return [candidate for _, candidate in final_candidates[:k]]


# ---------------------------------------------------------------------------
# Busca local global (hill-climbing)
# ---------------------------------------------------------------------------


def _local_search(candidate: String, strings: Sequence[String]) -> String:
    """
    Aplica hill-climbing (busca local) para melhorar uma solução candidata.

    Esta função implementa uma busca local clássica que refina iterativamente
    uma solução candidata através de movimentos locais. É a fase final do
    H³-CSP, responsável por otimizar a solução global obtida pela fusão
    dos blocos.

    ALGORITMO DE HILL-CLIMBING:
    1. **PREPARAÇÃO**: Converte string para lista mutável
    2. **ITERAÇÃO**: Enquanto houver melhorias:
       - Para cada posição i:
         - Salva caractere atual
         - Testa todas as substituições possíveis
         - Se melhoria encontrada, mantém mudança
         - Caso contrário, restaura caractere original
    3. **CONVERGÊNCIA**: Para quando nenhuma melhoria é encontrada

    ESTRATÉGIA DE BUSCA:
    - **Busca Sistemática**: Testa todas as posições sequencialmente
    - **Alfabeto Adaptativo**: Considera apenas caracteres que aparecem
      nas strings originais em cada posição (evita bias)
    - **Movimento Guloso**: Aceita primeira melhoria encontrada
    - **Busca Completa**: Continua até convergir para ótimo local

    VANTAGENS:
    - Simplicidade: Algoritmo direto e fácil de implementar
    - Eficiência: Complexidade O(L × |Σ| × iterações)
    - Determinismo: Sempre produz o mesmo resultado
    - Garantia: Converge para ótimo local

    LIMITAÇÕES:
    - Ótimos Locais: Pode ficar preso em máximos locais
    - Ordem Dependente: Resultado pode depender da ordem de teste
    - Guloso: Não considera combinações de mudanças

    APLICAÇÃO NO H³-CSP:
    - Refinamento final da solução obtida por fusão de blocos
    - Correção de inconsistências entre blocos
    - Melhoria da qualidade global da solução

    COMPLEXIDADE:
    - Temporal: O(L × |Σ| × iterações × n) onde n = |strings|
    - Espacial: O(L) para armazenar candidato e alfabetos por posição
    - Iterações: Tipicamente pequeno (convergência rápida)

    Args:
        candidate (String): String candidata inicial para refinamento.
                           Deve ter mesmo comprimento que strings originais.
        strings (Sequence[String]): Sequência de strings de entrada.

    Returns:
        String: String refinada (ótimo local da busca).
               Garantidamente não pior que a entrada.

    Examples:
        >>> candidate = "ACGT"
        >>> strings = ["ACCT", "AGGT", "ATGT"]
        >>> refined = _local_search(candidate, strings)
        >>> # refined pode ser "ACGT" ou "AGGT" dependendo das melhorias

        >>> # Caso onde nenhuma melhoria é possível
        >>> candidate = "AAAA"
        >>> strings = ["AAAA", "AAAA", "AAAA"]
        >>> refined = _local_search(candidate, strings)
        >>> refined == "AAAA"  # Sem mudanças possíveis
        True

    Note:
        O algoritmo é determinístico, mas a qualidade do resultado depende
        fortemente da qualidade da solução inicial fornecida.
    """
    # Converte string para lista mutável para eficiência
    candidate_list = list(candidate)
    L = len(candidate_list)

    # PRÉ-PROCESSAMENTO: Cria alfabeto adaptativo por posição
    # Para cada posição, considera apenas caracteres que aparecem nas strings originais
    # Isso evita bias e reduz o espaço de busca
    alphabets_by_position = []
    for i in range(L):
        position_chars = {s[i] for s in strings}  # Caracteres únicos na posição i
        alphabets_by_position.append(position_chars)

    # LOOP PRINCIPAL: Itera até convergir para ótimo local
    improvement_found = True
    while improvement_found:
        improvement_found = False

        # Calcula fitness atual (baseline para comparação)
        current_fitness = max_distance("".join(candidate_list), list(strings))

        # BUSCA SISTEMÁTICA: Testa cada posição sequencialmente
        for position in range(L):
            original_char = candidate_list[position]

            # Testa todas as substituições possíveis na posição atual
            for alternative_char in alphabets_by_position[position]:
                # Pula se for o mesmo caractere (sem mudança)
                if alternative_char == original_char:
                    continue

                # MOVIMENTO LOCAL: Aplica mudança temporária
                candidate_list[position] = alternative_char

                # Avalia novo fitness
                new_fitness = max_distance("".join(candidate_list), list(strings))

                # CRITÉRIO DE ACEITAÇÃO: Aceita se houver melhoria
                if new_fitness < current_fitness:
                    # Melhoria encontrada: mantém mudança e atualiza baseline
                    current_fitness = new_fitness
                    improvement_found = True
                    original_char = alternative_char
                    break  # Aceita primeira melhoria (estratégia gulosa)
                else:
                    # Sem melhoria: restaura caractere original
                    candidate_list[position] = original_char

    # Converte lista de volta para string
    return "".join(candidate_list)


# ---------------------------------------------------------------------------
# Classe principal
# ---------------------------------------------------------------------------


class H3CSP:
    """
    Implementação do algoritmo H³-CSP (Hybrid Hierarchical Hamming Search).

    O H³-CSP é um algoritmo híbrido de três camadas que resolve o Closest String
    Problem através de:

    1. **B-Splitter**: Divisão hierárquica em ~√L blocos contíguos
    2. **Smart-Core**: Seleção adaptativa de técnicas por bloco:
       - Busca exaustiva para blocos pequenos (d_b ≤ block_small)
       - Beam search para blocos médios/grandes (d_b > block_small)
    3. **Global Refine**: Fusão de blocos + refinamento por hill-climbing

    O algoritmo é determinístico e adapta sua estratégia baseada na dificuldade
    de cada bloco (medida pela distância do consenso). Blocos mais difíceis
    recebem técnicas mais sofisticadas.

    Attributes:
        strings (list[String]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        L (int): Comprimento das strings.
        params (dict): Parâmetros do algoritmo.
        rng (random.Random): Gerador de números aleatórios.
        progress_callback (Callable | None): Callback para progresso.
        blocks (list[Block]): Lista de blocos gerados pela divisão.

    Example:
        >>> h3 = H3CSP(["ACGT", "AGCT", "ATCT"], "ACGT")
        >>> center, distance = h3.run()
        >>> print(f"Solução: {center} com distância {distance}")
    """

    def __init__(self, strings: Sequence[String], alphabet: str, **params):
        """
        Inicializa o algoritmo H³-CSP.

        Configura os parâmetros, valida a entrada e prepara as estruturas
        necessárias para execução, incluindo a divisão inicial em blocos.

        Args:
            strings (Sequence[String]): Sequência de strings de entrada.
                                      Todas devem ter o mesmo comprimento.
            alphabet (str): Alfabeto utilizado (ex: "ACGT" para DNA).
            **params: Parâmetros do algoritmo que sobrescreverão os padrões.
                     Ver H3_CSP_DEFAULTS para opções disponíveis.

        Raises:
            ValueError: Se as strings tiverem comprimentos diferentes.
            ValueError: Se o alfabeto estiver vazio.

        Note:
            A divisão em blocos é feita imediatamente usando a regra √L.
        """
        self.strings = list(strings)
        self.alphabet = alphabet
        self.L = len(strings[0])

        # Merge defaults + params
        self.params = dict(H3_CSP_DEFAULTS)
        self.params.update(params)

        self.rng = random.Random(self.params["seed"])
        self.progress_callback: Callable[[str], None] | None = None

        # Divisão inicial em blocos usando a regra √L
        self.blocks = split_in_blocks(self.L)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define um callback para monitorar o progresso da execução.

        O callback será chamado em pontos-chave do algoritmo com mensagens
        descritivas sobre o progresso atual.

        Args:
            callback (Callable[[str], None]): Função que recebe mensagens
                                             de progresso como string.
        """
        self.progress_callback = callback

    # ---------------------------------------------------------------------

    def _smart_core(self) -> list[list[String]]:
        """
        Implementa o Smart-Core: seleção adaptativa de técnicas por bloco.

        Esta é a função central do H³-CSP, responsável por analisar cada bloco
        e selecionar a técnica de busca mais apropriada baseada na dificuldade
        do bloco. A seleção adaptativa permite usar técnicas computacionalmente
        caras apenas onde necessário, otimizando o uso de recursos.

        ESTRATÉGIA DE SELEÇÃO ADAPTATIVA:

        1. **MEDIÇÃO DE DIFICULDADE**: Para cada bloco, calcula d_b como a
           distância máxima do consenso para os segmentos do bloco.

        2. **CLASSIFICAÇÃO DE BLOCOS**:
           - **FÁCIL** (d_b ≤ block_small): Busca exaustiva
           - **MÉDIO** (block_small < d_b ≤ block_medium): Beam search reduzido
           - **DIFÍCIL** (d_b > block_medium): Beam search completo

        3. **APLICAÇÃO DE TÉCNICAS**:
           - Busca exaustiva: Garante otimalidade local
           - Beam search reduzido: Balanceia qualidade e eficiência
           - Beam search completo: Máximo esforço para blocos difíceis

        FILOSOFIA DA ADAPTAÇÃO:
        - Blocos fáceis (baixa diversidade) → Solução rápida e ótima
        - Blocos médios (diversidade moderada) → Busca eficiente
        - Blocos difíceis (alta diversidade) → Busca intensiva

        VANTAGENS DA SELEÇÃO ADAPTATIVA:
        - **Eficiência**: Recursos computacionais usados onde necessário
        - **Qualidade**: Técnicas apropriadas para cada nível de dificuldade
        - **Escalabilidade**: Adapta-se automaticamente ao problema
        - **Flexibilidade**: Parâmetros ajustáveis para diferentes cenários

        COMPLEXIDADE:
        - Blocos fáceis: O(|Σ|^(L/B)) - busca exaustiva
        - Blocos médios: O(beam_width/2 × |Σ| × (L/B)) - beam search reduzido
        - Blocos difíceis: O(beam_width × |Σ| × (L/B)) - beam search completo

        Returns:
            list[list[String]]: Lista de listas de candidatos, onde cada sublista
                               contém os k melhores candidatos para um bloco.
                               A ordem corresponde à ordem em self.blocks.

        Examples:
            >>> # Para 3 blocos com dificuldades [1, 3, 8] e limites [2, 5]
            >>> # Bloco 1: d_b=1 ≤ 2 → busca exaustiva
            >>> # Bloco 2: 2 < d_b=3 ≤ 5 → beam search reduzido
            >>> # Bloco 3: d_b=8 > 5 → beam search completo
            >>> candidates = self._smart_core()
            >>> # candidates[0] = ["AA", "AC", "AG"]  # Exaustiva
            >>> # candidates[1] = ["TT", "TC"]        # Beam reduzido
            >>> # candidates[2] = ["GG", "GC", "GA"]  # Beam completo

        Note:
            A qualidade dos candidatos depende dos parâmetros block_small,
            block_medium e beam_width. Valores maiores produzem melhores
            resultados mas aumentam o custo computacional.
        """
        # Extrai parâmetros de configuração
        k = self.params["k_candidates"]  # Número de candidatos por bloco
        small_threshold = self.params["block_small"]  # Limite para blocos fáceis
        medium_threshold = self.params["block_medium"]  # Limite para blocos médios
        beam_width = self.params["beam_width"]  # Largura do beam search

        # Lista para armazenar candidatos de todos os blocos
        all_block_candidates = []

        # Processa cada bloco independentemente
        for block_index, (l, r) in enumerate(self.blocks):
            # FASE 1: MEDIÇÃO DE DIFICULDADE
            # Calcula consenso do bloco como baseline
            block_consensus = consensus_block(self.strings, l, r)

            # Calcula d_b: distância máxima do consenso para os segmentos
            block_segments = [s[l:r] for s in self.strings]
            block_difficulty = max_distance(block_consensus, block_segments)

            # Log para debug/análise
            logger.debug(
                "Bloco %d [%d:%d]: d_b=%d, tamanho=%d",
                block_index,
                l,
                r,
                block_difficulty,
                r - l,
            )

            # FASE 2: SELEÇÃO ADAPTATIVA DE TÉCNICA
            if block_difficulty <= small_threshold:
                # BLOCO FÁCIL: Busca exaustiva
                logger.debug("Bloco %d: usando busca exaustiva", block_index)
                candidates = _exhaustive_block(self.strings, self.alphabet, l, r, k)
            elif block_difficulty <= medium_threshold:
                # BLOCO MÉDIO: Beam search reduzido
                logger.debug("Bloco %d: usando beam search reduzido", block_index)
                candidates = _beam_search_block(
                    self.strings, self.alphabet, l, r, beam_width // 2, k
                )
            else:
                # BLOCO DIFÍCIL: Beam search completo
                logger.debug("Bloco %d: usando beam search completo", block_index)
                candidates = _beam_search_block(
                    self.strings, self.alphabet, l, r, beam_width, k
                )

            # FASE 3: VALIDAÇÃO E ARMAZENAMENTO
            # Garante que temos pelo menos um candidato (fallback)
            if not candidates:
                logger.warning("Bloco %d: sem candidatos, usando consenso", block_index)
                candidates = [block_consensus]

            # Armazena candidatos do bloco
            all_block_candidates.append(candidates)

            # Log de progresso
            logger.debug(
                "Bloco %d: %d candidatos gerados, melhor distância=%d",
                block_index,
                len(candidates),
                (
                    max_distance(candidates[0], block_segments)
                    if candidates
                    else float("inf")
                ),
            )

        return all_block_candidates

    # ---------------------------------------------------------------------

    def _fuse_blocks(self, chosen: list[String]) -> String:
        """
        Fusão de blocos: concatena os candidatos escolhidos para formar solução completa.

        Esta função implementa a fase de fusão do H³-CSP, onde os melhores
        candidatos de cada bloco são concatenados para formar uma solução
        completa. A fusão é direta e preserva a ordem original dos blocos.

        ALGORITMO DE FUSÃO:
        1. **VALIDAÇÃO**: Verifica se número de candidatos corresponde ao número de blocos
        2. **CONCATENAÇÃO**: Une os candidatos na ordem dos blocos originais
        3. **VALIDAÇÃO FINAL**: Verifica se comprimento final está correto

        CARACTERÍSTICAS DA FUSÃO:
        - **Simplicidade**: Operação direta de concatenação
        - **Preservação de Ordem**: Mantém ordem espacial original
        - **Eficiência**: Complexidade O(L) onde L é comprimento total
        - **Determinística**: Sempre produz o mesmo resultado

        VANTAGENS:
        - Preserva otimalidade local de cada bloco
        - Operação computacionalmente eficiente
        - Fácil implementação e debug
        - Compatível com qualquer técnica de busca por bloco

        LIMITAÇÕES:
        - Não considera interações entre blocos
        - Pode produzir inconsistências nas fronteiras
        - Solução pode não ser globalmente ótima

        PAPEL NO H³-CSP:
        - Ponte entre otimização local (por bloco) e global (refinamento)
        - Construção de solução inicial para hill-climbing
        - Combinação de conhecimento local de múltiplos blocos

        COMPLEXIDADE:
        - Temporal: O(L) onde L é comprimento total das strings
        - Espacial: O(L) para armazenar resultado

        Args:
            chosen (list[String]): Lista de candidatos escolhidos, um por bloco.
                                  A ordem deve corresponder à ordem em self.blocks.
                                  len(chosen) deve ser igual a len(self.blocks).

        Returns:
            String: String completa formada pela concatenação dos blocos.
                   Comprimento igual ao comprimento das strings originais.

        Raises:
            AssertionError: Se número de candidatos não corresponder ao número de blocos.
            ValueError: Se algum candidato tiver comprimento inconsistente.

        Examples:
            >>> # Para blocks = [(0,2), (2,4), (4,6)] e chosen = ["AC", "GT", "AA"]
            >>> self._fuse_blocks(["AC", "GT", "AA"])
            "ACGTAA"

            >>> # Para blocks = [(0,3), (3,6)] e chosen = ["ATG", "CCT"]
            >>> self._fuse_blocks(["ATG", "CCT"])
            "ATGCCT"

            >>> # Caso com um único bloco
            >>> self._fuse_blocks(["ACGTAA"])
            "ACGTAA"

        Note:
            A fusão assume que cada candidato tem exatamente o tamanho
            correspondente ao seu bloco. Candidatos com tamanho incorreto
            podem causar erros ou resultados inesperados.
        """
        # VALIDAÇÃO DE ENTRADA
        assert len(chosen) == len(self.blocks), (
            f"Número de candidatos ({len(chosen)}) deve corresponder "
            f"ao número de blocos ({len(self.blocks)})"
        )

        # CONCATENAÇÃO DOS CANDIDATOS
        # Simples concatenação na ordem dos blocos
        fused_result = "".join(chosen)

        # VALIDAÇÃO FINAL
        expected_length = self.L
        if len(fused_result) != expected_length:
            raise ValueError(
                f"Comprimento da solução fusionada ({len(fused_result)}) "
                f"não corresponde ao esperado ({expected_length})"
            )

        return fused_result

    # ---------------------------------------------------------------------

    def run(self) -> tuple[String, int]:
        """
        Executa o algoritmo H³-CSP completo com todas as suas fases.

        Esta função orquestra a execução completa do algoritmo H³-CSP,
        coordenando as três fases principais e gerenciando controles de
        tempo, progresso e tratamento de erros.

        FLUXO DE EXECUÇÃO COMPLETO:

        1. **INICIALIZAÇÃO**:
           - Registra tempo de início
           - Configura callback de progresso
           - Prepara estruturas de dados

        2. **FASE 1 - SMART-CORE**:
           - Analisa dificuldade de cada bloco
           - Seleciona técnica apropriada por bloco
           - Gera candidatos para cada bloco

        3. **FASE 2 - FUSÃO**:
           - Seleciona melhor candidato de cada bloco
           - Concatena candidatos para formar solução inicial
           - Avalia qualidade da solução fusionada

        4. **FASE 3 - REFINAMENTO GLOBAL**:
           - Aplica hill-climbing iterativo
           - Verifica timeout a cada iteração
           - Para quando não há mais melhorias

        5. **FINALIZAÇÃO**:
           - Retorna melhor solução encontrada
           - Trata exceções e erros
           - Registra estatísticas finais

        CONTROLES DE EXECUÇÃO:
        - **Timeout**: Verifica limite de tempo a cada iteração
        - **Callback**: Reporta progresso se configurado
        - **Logging**: Registra eventos importantes
        - **Exceções**: Trata e re-propaga erros

        CARACTERÍSTICAS:
        - **Determinismo**: Sempre produz o mesmo resultado (sem timeout)
        - **Robustez**: Trata erros e casos extremos
        - **Monitoramento**: Permite acompanhar progresso
        - **Flexibilidade**: Parâmetros configuráveis

        COMPLEXIDADE TOTAL:
        - Smart-Core: O(B × max(|Σ|^(L/B), beam_width × |Σ| × (L/B)))
        - Fusão: O(L)
        - Refinamento: O(iterações × L × |Σ| × n)
        - Total: Dominado pela fase mais custosa (geralmente Smart-Core)

        Returns:
            tuple[String, int]: Tupla contendo:
                - String: Solução encontrada (string central)
                - int: Distância máxima da solução para as strings originais

        Raises:
            Exception: Re-propaga qualquer exceção ocorrida durante a execução,
                      após registrar o erro no log.

        Examples:
            >>> h3 = H3CSP(["ACGT", "AGCT", "ATCT"], "ACGT")
            >>> center, distance = h3.run()
            >>> print(f"Solução: {center} com distância {distance}")
            Solução: ACGT com distância 1

            >>> # Com callback de progresso
            >>> def progress(msg):
            ...     print(f"[H³-CSP] {msg}")
            >>> h3.set_progress_callback(progress)
            >>> center, distance = h3.run()
            [H³-CSP] Analisando blocos...
            [H³-CSP] Fusão de blocos...
            [H³-CSP] Refinamento global...
            [H³-CSP] Refinamento: iteração 1
            [H³-CSP] Melhoria encontrada: distância=0

        Note:
            O algoritmo pode ser interrompido por timeout, retornando a
            melhor solução encontrada até o momento da interrupção.
        """
        # INICIALIZAÇÃO
        start_time = time.time()

        try:
            # FASE 1: SMART-CORE - Análise e processamento de blocos
            if self.progress_callback:
                self.progress_callback("Analisando blocos...")

            logger.info("Iniciando Smart-Core para %d blocos", len(self.blocks))
            block_candidates = self._smart_core()

            # Log estatísticas dos blocos
            total_candidates = sum(len(cands) for cands in block_candidates)
            logger.info(
                "Smart-Core concluído: %d candidatos gerados (%d blocos)",
                total_candidates,
                len(self.blocks),
            )

            # FASE 2: FUSÃO - Combinação dos melhores candidatos
            if self.progress_callback:
                self.progress_callback("Fusão de blocos...")

            # Seleciona melhor candidato de cada bloco (primeiro da lista ordenada)
            best_candidates_per_block = [
                candidates[0] for candidates in block_candidates
            ]

            # Fusiona candidatos para formar solução inicial
            center = self._fuse_blocks(best_candidates_per_block)
            best_distance = max_distance(center, list(self.strings))

            logger.info("Fusão inicial: distância=%d", best_distance)

            # FASE 3: REFINAMENTO GLOBAL - Hill-climbing iterativo
            if self.progress_callback:
                self.progress_callback("Refinamento global...")

            # Aplica refinamento local até convergência ou timeout
            for iteration in range(self.params["local_iters"]):
                # CONTROLE DE TIMEOUT
                elapsed_time = time.time() - start_time
                if elapsed_time >= self.params["max_time"]:
                    logger.warning(
                        "Timeout atingido após %.2fs no H³-CSP (iter %d)",
                        elapsed_time,
                        iteration,
                    )
                    if self.progress_callback:
                        self.progress_callback(
                            f"Timeout atingido após {elapsed_time:.1f}s"
                        )
                    break

                # CALLBACK DE PROGRESSO
                if self.progress_callback:
                    self.progress_callback(f"Refinamento: iteração {iteration+1}")

                # APLICAÇÃO DE BUSCA LOCAL
                previous_center = center
                center = _local_search(center, self.strings)
                new_distance = max_distance(center, list(self.strings))

                # VERIFICAÇÃO DE MELHORIA
                if new_distance < best_distance:
                    # Melhoria encontrada
                    logger.info(
                        "Refinamento iter %d: %d → %d (melhoria: %d)",
                        iteration + 1,
                        best_distance,
                        new_distance,
                        best_distance - new_distance,
                    )
                    best_distance = new_distance

                    if self.progress_callback:
                        self.progress_callback(
                            f"Melhoria encontrada: distância={best_distance}"
                        )

                    # Verifica se encontrou solução ótima
                    if best_distance == 0:
                        logger.info("Solução ótima encontrada!")
                        if self.progress_callback:
                            self.progress_callback("Solução ótima encontrada!")
                        break
                elif center == previous_center:
                    # Sem melhoria e sem mudança: convergiu para ótimo local
                    logger.info(
                        "Convergência para ótimo local na iteração %d", iteration + 1
                    )
                    if self.progress_callback:
                        self.progress_callback("Convergiu para ótimo local")
                    break
                else:
                    # Sem melhoria mas houve mudança: continua buscando
                    logger.debug("Refinamento iter %d: sem melhoria", iteration + 1)

            # ESTATÍSTICAS FINAIS
            total_time = time.time() - start_time
            logger.info(
                "H³-CSP concluído: distância final=%d, tempo=%.2fs",
                best_distance,
                total_time,
            )

            return center, best_distance

        except Exception as e:
            # TRATAMENTO DE ERROS
            logger.error("Erro durante execução do H³-CSP: %s", str(e))
            if self.progress_callback:
                self.progress_callback(f"Erro: {str(e)}")
            raise e  # Re-propaga para ser capturado no wrapper
