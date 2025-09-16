"""
Genetic Operations for BLF-GA (Block Letter Frequency Genetic Algorithm)

This module implements a complete suite of genetic operators specialized for the
Closest String Problem (CSP) using genetic algorithms. BLF-GA is an evolutionary algorithm that
evolves a population of candidate strings to find the consensus string that minimizes the
maximum distance to a set of input strings.

MODULE ARCHITECTURE:
===================
1. Diversity Measurement:
    - Calculation of average Hamming distance between individuals
    - Population convergence monitoring

2. Mutation Operators:
    - Multi-point mutation: Randomly changes multiple positions
    - Inversion mutation: Inverts string segments
    - Transposition mutation: Moves segments to other positions

3. Crossover Operators:
    - One-point crossover: Exchanges suffixes between parents
    - Uniform crossover: Combines genes position by position
    - Block crossover: Exchanges structural segments

4. Local Refinement:
    - Greedy search: Position-by-position optimization
    - Swap refinement: Permutes positions to improve fitness
    - Insertion refinement: Moves small segments
    - 2-opt refinement: Inverts segments for optimization

ALGORITHMIC STRATEGIES:
======================
- Each individual is represented as a string of the same length as input strings
- Fitness is calculated as the maximum distance from individual to all reference strings
- Operators are designed to preserve string length
- Local refinement is used to escape local optima and accelerate convergence

USE CASES:
=========
- Consensus search in biological sequences
- String optimization in alignment problems
- Genetic variability analysis
- Distance-based clustering problems
"""

import random

import numpy as np

String = str
Population = list[String]


# =============================================================================
# POPULATION DIVERSITY MEASUREMENT
# =============================================================================


def mean_hamming_distance(pop: Population, distance_func=None) -> float:
    """
    Calculate the average Hamming distance between all pairs of individuals in the population.

    Hamming distance is a fundamental metric for measuring genetic diversity
    in string populations. It counts the number of positions where two individuals differ.
    A population with high diversity (high average distance) indicates greater exploration
    of the search space, while low diversity may indicate premature convergence.

    ALGORITHMIC STRATEGY:
    - Uses standardized distance calculation method when provided
    - Falls back to direct calculation only if no distance function provided
    - Considers only unique pairs to avoid duplication
    - Complexity: O(n²·m) where n is population size and m is string length

    Args:
        pop (Population): List of individuals (strings) in the population
        distance_func: Distance calculation function (should be HammingDistanceCalculator.distance)

    Returns:
        float: Average Hamming distance between all pairs.
               Returns 0.0 if population has fewer than 2 individuals.

    Examples:
        >>> from src.domain.distance import HammingDistanceCalculator
        >>> pop = ["ACGT", "ATGT", "GCGT"]
        >>> mean_hamming_distance(pop, HammingDistanceCalculator.distance)
        1.3333333333333333

        >>> pop = ["AAAA", "AAAA"]  # Homogeneous population
        >>> mean_hamming_distance(pop, HammingDistanceCalculator.distance)
        0.0

    Note:
        - Assumes all strings have the same length
        - Useful for monitoring genetic algorithm convergence
        - Values close to 0 indicate converged population
        - Uses standardized distance calculation when available
    """
    if len(pop) < 2:
        return 0.0

    # Use standardized distance calculation if provided
    if distance_func is not None:
        total_distance = 0.0
        total_pairs = 0

        for i in range(len(pop)):
            for j in range(i + 1, len(pop)):
                total_distance += distance_func(pop[i], pop[j])
                total_pairs += 1

        return total_distance / total_pairs if total_pairs > 0 else 0.0

    # Fallback to direct calculation (for backward compatibility)
    arr = np.array([list(ind) for ind in pop])
    dists = np.sum(arr[:, None, :] != arr[None, :, :], axis=2)
    iu = np.triu_indices(len(pop), 1)
    return np.mean(dists[iu])


# =============================================================================
# MUTATION OPERATORS
# =============================================================================


def mutate_multi(ind: str, alphabet: str, rng: random.Random, n: int = 2) -> str:
    """
    Realiza mutação multi-ponto alterando até n posições aleatórias do indivíduo.

    Este operador implementa mutação pontual múltipla, onde cada posição selecionada
    é alterada para um símbolo diferente do alfabeto. É uma estratégia de diversificação
    que introduz variabilidade genética controlada na população.

    ESTRATÉGIA ALGORÍTMICA:
    - Seleciona n posições aleatórias (com possível repetição)
    - Para cada posição, escolhe um símbolo diferente do atual
    - Preserva o comprimento da string
    - Evita mutações "neutras" (mesmo símbolo)

    Args:
        ind (str): String original (indivíduo a ser mutado)
        alphabet (str): Conjunto de símbolos válidos para mutação
        rng (random.Random): Instância do gerador de números aleatórios
        n (int, optional): Número de mutações a aplicar. Defaults to 2.

    Returns:
        str: String mutada com até n posições alteradas

    Examples:
        >>> rng = random.Random(42)
        >>> mutate_multi("ACGT", "ACGT", rng, n=1)
        'ACTT'  # Posição 2 alterada de G para T

        >>> mutate_multi("AAAA", "ACGT", rng, n=2)
        'ACGA'  # Duas posições alteradas

    Note:
        - Se o alfabeto contém apenas um símbolo, não haverá mutação
        - Posições podem ser selecionadas múltiplas vezes
        - Útil para escapar de ótimos locais
    """
    chars = list(ind)
    L = len(chars)

    # Aplica n mutações
    for _ in range(n):
        # Seleciona posição aleatória
        pos = rng.randint(0, L - 1)
        old = chars[pos]

        # Escolhe símbolo diferente do atual
        choices = [c for c in alphabet if c != old]
        if choices:
            chars[pos] = rng.choice(choices)

    return "".join(chars)


def mutate_inversion(ind: str, rng: random.Random) -> str:
    """
    Realiza mutação por inversão de segmento.

    Seleciona dois pontos aleatórios na string e inverte (reverte) o segmento entre eles.
    Este operador é inspirado em inversões cromossômicas da biologia evolutiva e é
    eficaz para reorganizar estruturas locais mantendo o conteúdo global.

    ESTRATÉGIA ALGORÍTMICA:
    - Seleciona dois pontos de corte aleatórios
    - Inverte o segmento entre os pontos (inclusivo)
    - Preserva o comprimento e conteúdo da string
    - Reorganiza apenas a ordem dos símbolos

    Args:
        ind (str): String original a ser mutada
        rng (random.Random): Instância do gerador de números aleatórios

    Returns:
        str: String com segmento invertido

    Examples:
        >>> rng = random.Random(42)
        >>> mutate_inversion("ABCDEF", rng)
        'AEDCBF'  # Segmento BCDE invertido para EDCB

        >>> mutate_inversion("ACGT", rng)
        'AGCT'  # Segmento CG invertido para GC

    Note:
        - Útil para escapar de ótimos locais através de rearranjos
        - Mantém a composição de símbolos inalterada
        - Pode produzir grandes mudanças estruturais com uma única operação
    """
    chars = list(ind)
    L = len(chars)

    # Seleciona dois pontos e ordena
    a, b = sorted(rng.sample(range(L), 2))

    # Inverte o segmento entre os pontos (inclusivo)
    chars[a : b + 1] = chars[a : b + 1][::-1]

    return "".join(chars)


def mutate_transposition(ind: str, rng: random.Random) -> str:
    """
    Realiza mutação por transposição de segmento.

    Seleciona um segmento aleatório da string e o move para uma nova posição.
    Este operador simula transposições genéticas e é eficaz para reorganizar
    estruturas mantendo segmentos intactos.

    ESTRATÉGIA ALGORÍTMICA:
    - Seleciona dois pontos para definir o segmento
    - Remove o segmento da posição original
    - Insere o segmento em uma nova posição aleatória
    - Preserva o comprimento e conteúdo da string

    Args:
        ind (str): String original a ser mutada
        rng (random.Random): Instância do gerador de números aleatórios

    Returns:
        str: String com segmento transposto

    Examples:
        >>> rng = random.Random(42)
        >>> mutate_transposition("ABCDEF", rng)
        'ADEFBC'  # Segmento BC movido para o final

        >>> mutate_transposition("ACGT", rng)
        'CAGT'  # Segmento C movido para posição 1

    Note:
        - Útil para reorganizar estruturas locais
        - Mantém a composição de símbolos inalterada
        - Pode produzir mudanças estruturais significativas
        - Eficaz para problemas onde a ordem dos elementos importa
    """
    chars = list(ind)
    L = len(chars)

    # Seleciona dois pontos para definir o segmento
    a, b = sorted(rng.sample(range(L), 2))

    # Extrai o segmento
    seg = chars[a : b + 1]

    # Remove o segmento da posição original
    del chars[a : b + 1]

    # Insere o segmento em nova posição aleatória
    pos = rng.randint(0, len(chars))
    chars[pos:pos] = seg

    return "".join(chars)


# =============================================================================
# OPERADORES DE CROSSOVER (RECOMBINAÇÃO)
# =============================================================================


def crossover_one_point(
    p1: String, p2: String, rng: random.Random
) -> tuple[String, String]:
    """
    Realiza crossover de um ponto entre dois indivíduos pais.

    O crossover de um ponto é o operador de recombinação mais simples e clássico.
    Corta ambos os pais em um ponto aleatório e troca os sufixos, criando dois filhos
    que combinam características de ambos os pais.

    ESTRATÉGIA ALGORÍTMICA:
    - Seleciona um ponto de corte aleatório (evita extremos)
    - Combina prefixo do pai 1 com sufixo do pai 2
    - Combina prefixo do pai 2 com sufixo do pai 1
    - Preserva o comprimento das strings

    Args:
        p1 (String): Primeiro pai
        p2 (String): Segundo pai
        rng (random.Random): Instância do gerador de números aleatórios

    Returns:
        tuple[String, String]: Dois filhos resultantes do crossover

    Examples:
        >>> rng = random.Random(42)
        >>> crossover_one_point("ABCD", "EFGH", rng)
        ('ABGH', 'EFCD')  # Corte na posição 2

        >>> crossover_one_point("ACGT", "TGCA", rng)
        ('ACCA', 'TGGT')  # Corte na posição 2

    Note:
        - Preserva blocos contíguos de genes
        - Baixa disruptividade (mantém estruturas locais)
        - Pode não ser eficaz se genes importantes estão separados
        - Ponto de corte evita posições extremas (1 a L-1)
    """
    L = len(p1)

    # Seleciona ponto de corte (evita extremos)
    point = rng.randint(1, L - 1)

    # Cria filhos trocando sufixos
    c1 = p1[:point] + p2[point:]
    c2 = p2[:point] + p1[point:]

    return c1, c2


def crossover_uniform(
    p1: String, p2: String, rng: random.Random
) -> tuple[String, String]:
    """
    Realiza crossover uniforme entre dois indivíduos pais.

    O crossover uniforme examina cada posição independentemente e decide
    aleatoriamente (50% de chance) de qual pai herdar o gene. Este operador
    é mais disruptivo que o crossover de um ponto, mas permite maior
    recombinação entre genes distantes.

    ESTRATÉGIA ALGORÍTMICA:
    - Para cada posição, gera um número aleatório
    - Se < 0.5, filho 1 herda do pai 1, filho 2 herda do pai 2
    - Se ≥ 0.5, filho 1 herda do pai 2, filho 2 herda do pai 1
    - Preserva o comprimento das strings

    Args:
        p1 (String): Primeiro pai
        p2 (String): Segundo pai
        rng (random.Random): Instância do gerador de números aleatórios

    Returns:
        tuple[String, String]: Dois filhos resultantes do crossover

    Examples:
        >>> rng = random.Random(42)
        >>> crossover_uniform("ABCD", "EFGH", rng)
        ('ABGH', 'EFCD')  # Herança posição por posição

        >>> crossover_uniform("ACGT", "TGCA", rng)
        ('ACCA', 'TGGT')  # Cada posição decidida independentemente

    Note:
        - Alta disruptividade (pode quebrar estruturas locais)
        - Excelente para recombinar genes distantes
        - Maior exploração do espaço de busca
        - Pode ser menos eficaz para problemas com dependências locais
    """
    L = len(p1)
    c1 = []
    c2 = []

    # Decide herança para cada posição
    for i in range(L):
        if rng.random() < 0.5:
            # Herança direta
            c1.append(p1[i])
            c2.append(p2[i])
        else:
            # Herança cruzada
            c1.append(p2[i])
            c2.append(p1[i])

    return "".join(c1), "".join(c2)


def crossover_blend_blocks(
    p1: String, p2: String, blocks: list[tuple[int, int]], rng: random.Random
) -> tuple[String, String]:
    """
    Realiza crossover por blocos estruturais entre dois indivíduos pais.

    Este operador respeita a estrutura em blocos do problema, trocando segmentos
    completos entre os pais. É útil quando o problema tem estrutura natural em
    blocos ou quando queremos preservar certas regiões funcionais.

    ESTRATÉGIA ALGORÍTMICA:
    - Para cada bloco definido, decide aleatoriamente (50% chance) se trocar
    - Troca segmentos completos entre os pais
    - Preserva a estrutura interna dos blocos
    - Permite recombinação respeitando fronteiras naturais

    Args:
        p1 (String): Primeiro pai
        p2 (String): Segundo pai
        blocks (list[tuple[int, int]]): Lista de blocos definidos como (início, fim)
        rng (random.Random): Instância do gerador de números aleatórios

    Returns:
        tuple[String, String]: Dois filhos resultantes do crossover

    Examples:
        >>> rng = random.Random(42)
        >>> blocks = [(0, 2), (2, 4)]  # Dois blocos
        >>> crossover_blend_blocks("ABCD", "EFGH", blocks, rng)
        ('EFCD', 'ABGH')  # Primeiro bloco trocado

        >>> crossover_blend_blocks("ACGT", "TGCA", [(0, 2), (2, 4)], rng)
        ('TGGT', 'ACCA')  # Blocos trocados independentemente

    Note:
        - Preserva estrutura funcional dos blocos
        - Útil para problemas com regiões conservadas
        - Permite controle fino sobre recombinação
        - Blocos devem ser não-sobrepostos e cobrir toda a string
    """
    c1 = list(p1)
    c2 = list(p2)

    # Para cada bloco, decide se trocar
    for l, r in blocks:
        if rng.random() < 0.5:
            # Troca segmentos do bloco
            c1[l:r] = p2[l:r]
            c2[l:r] = p1[l:r]

    return "".join(c1), "".join(c2)


# =============================================================================
# OPERADORES DE REFINAMENTO LOCAL
# =============================================================================


def refine_greedy(ind: String, strings: list[String], max_distance_func=None) -> String:
    """
    Aplica refinamento local greedy posição por posição.

    Este refinador testa cada posição da string, experimentando todos os
    símbolos do alfabeto e mantendo aquele que produz melhor fitness.
    É um método de busca local que garante melhoria monotônica.

    ESTRATÉGIA ALGORÍTMICA:
    - Itera sobre cada posição da string
    - Para cada posição, testa todos os símbolos do alfabeto extraído das strings
    - Mantém a alteração que resulta no melhor fitness (menor distância máxima)
    - Repete até que nenhuma posição possa ser melhorada
    - Usa greedy: aceita a primeira melhoria encontrada

    Args:
        ind (String): Indivíduo a ser refinado
        strings (list[String]): Lista de strings de referência para cálculo do fitness
        max_distance_func: Function to calculate max distance (should accept individual and strings)

    Returns:
        String: Indivíduo refinado através de busca greedy posição por posição

    Examples:
        >>> strings = ["ACGT", "ATGT", "AGGT"]
        >>> refine_greedy("AAAA", strings, some_max_distance_func)
        'ACGT'  # Refinado para minimizar distância máxima

        >>> refine_greedy("ACGT", ["ACGT", "ACGT"], some_max_distance_func)
        'ACGT'  # Já é ótimo, não há mudança

    Note:
        - Garante melhoria monotônica do fitness
        - Pode ficar preso em ótimos locais
        - Complexidade O(k * n * |Σ|) onde k é número de iterações, n é comprimento, |Σ| é tamanho do alfabeto
        - Muito eficaz para refinamento final de soluções
        - Uses standardized distance calculation when available
    """
    if max_distance_func is None:
        # Fallback for backward compatibility
        from src.domain.distance import HammingDistanceCalculator

        distance_calc = HammingDistanceCalculator(strings)
        max_distance_func = distance_calc.max_distance

    # Extrai alfabeto das strings de referência
    alphabet = set()
    for s in strings:
        alphabet.update(s)
    alphabet = sorted(alphabet)

    best_ind = ind
    best_fitness = max_distance_func(best_ind)
    improved = True

    # Continua até não haver melhoria
    while improved:
        improved = False
        current = list(best_ind)

        # Tenta melhorar cada posição
        for pos in range(len(current)):
            original_char = current[pos]
            best_char = original_char

            # Testa cada símbolo do alfabeto
            for char in alphabet:
                if char != original_char:
                    current[pos] = char
                    test_ind = "".join(current)
                    test_fitness = max_distance_func(test_ind)

                    # Se encontrou melhoria, atualiza
                    if test_fitness < best_fitness:
                        best_fitness = test_fitness
                        best_char = char
                        improved = True

            # Aplica a melhor mudança encontrada
            current[pos] = best_char

        best_ind = "".join(current)

    return best_ind


def refine_swap(ind: String, strings: list[String], max_distance_func=None) -> String:
    """
    Aplica refinamento local por troca de posições.

    Este operador tenta trocar cada par de posições da string para verificar
    se a troca melhora o fitness. É uma busca local que preserva a composição
    de símbolos mas reorganiza sua distribuição.

    ESTRATÉGIA ALGORÍTMICA:
    - Enumera todos os pares de posições possíveis
    - Para cada par, tenta a troca e avalia o fitness
    - Mantém a troca se melhora o fitness
    - Desfaz a troca se não há melhoria
    - Usa abordagem gulosa (primeira melhoria encontrada)

    Args:
        ind (String): Indivíduo a ser refinado
        strings (list[String]): Lista de strings de referência para cálculo do fitness

    Returns:
        String: Indivíduo refinado através de trocas de posições

    Examples:
        >>> strings = ["ACGT", "ATCG", "AGCT"]
        >>> refine_swap("TGCA", strings)
        'ACGT'  # Posições trocadas para melhorar alinhamento

        >>> refine_swap("ACGT", ["ACGT"])
        'ACGT'  # Já é ótimo, não há mudança

    Note:
        - Preserva a composição de símbolos da string
        - Complexidade O(n²) onde n é o comprimento da string
        - Útil quando a ordem dos símbolos importa mais que sua identidade
        - Pode ser combinado com outros operadores de refinamento
    """
    if max_distance_func is None:
        # Fallback for backward compatibility
        from src.domain.distance import HammingDistanceCalculator

        distance_calc = HammingDistanceCalculator(strings)
        max_distance_func = distance_calc.max_distance

    best_ind = ind
    best_fitness = max_distance_func(best_ind)
    current = list(best_ind)
    L = len(current)

    # Tenta todas as trocas de pares
    for i in range(L):
        for j in range(i + 1, L):
            # Troca posições i e j
            current[i], current[j] = current[j], current[i]
            test_ind = "".join(current)
            test_fitness = max_distance_func(test_ind)

            # Se melhorou, mantém a troca
            if test_fitness < best_fitness:
                best_fitness = test_fitness
                best_ind = test_ind
            else:
                # Desfaz a troca se não melhorou
                current[i], current[j] = current[j], current[i]

    return best_ind


def refine_insertion(
    ind: String, strings: list[String], max_distance_func=None
) -> String:
    """
    Aplica refinamento local por inserção/remoção de segmentos.

    Este operador tenta mover segmentos pequenos (1-3 caracteres) para
    diferentes posições da string. É eficaz para reorganizar estruturas
    locais e encontrar melhores arranjos de símbolos.

    ESTRATÉGIA ALGORÍTMICA:
    - Testa segmentos de tamanho 1, 2 e 3
    - Para cada segmento, tenta todas as posições de inserção
    - Mantém o comprimento original da string
    - Avalia fitness após cada movimento
    - Mantém a melhor configuração encontrada

    Args:
        ind (String): Indivíduo a ser refinado
        strings (list[String]): Lista de strings de referência para cálculo do fitness

    Returns:
        String: Indivíduo refinado através de movimentação de segmentos

    Examples:
        >>> strings = ["ABCDEF", "ABCFED", "ABCFDE"]
        >>> refine_insertion("ABCDEF", strings)
        'ABCFED'  # Segmento "EF" movido para melhor posição

        >>> refine_insertion("ACGT", ["ACGT"])
        'ACGT'  # Já é ótimo, não há mudança

    Note:
        - Preserva a composição de símbolos da string
        - Complexidade O(n³) onde n é o comprimento da string
        - Útil para reorganizar estruturas locais
        - Trabalha com segmentos pequenos para evitar disrupção excessiva
        - Uses standardized distance calculation when available
    """
    if max_distance_func is None:
        # Fallback for backward compatibility
        from src.domain.distance import HammingDistanceCalculator

        distance_calc = HammingDistanceCalculator(strings)
        max_distance_func = distance_calc.max_distance

    best_ind = ind
    best_fitness = max_distance_func(best_ind)
    current = list(best_ind)
    L = len(current)

    # Testa mover segmentos de tamanho 1-3
    for seg_len in range(1, min(4, L)):
        for start in range(L - seg_len + 1):
            # Extrai o segmento
            segment = current[start : start + seg_len]

            # Remove segmento da posição original
            temp = current[:start] + current[start + seg_len :]

            # Tenta inserir em cada posição possível
            for insert_pos in range(len(temp) + 1):
                new_current = temp[:insert_pos] + segment + temp[insert_pos:]

                # Verifica se mantém comprimento original
                if len(new_current) == L:
                    test_ind = "".join(new_current)
                    test_fitness = max_distance_func(test_ind)

                    # Se melhorou, atualiza
                    if test_fitness < best_fitness:
                        best_fitness = test_fitness
                        best_ind = test_ind
                        current = new_current

    return best_ind


def refine_2opt(ind: String, strings: list[String], max_distance_func=None) -> String:
    """
    Aplica refinamento local usando estratégia 2-opt.

    O refinamento 2-opt é inspirado no problema do caixeiro viajante e
    funciona invertendo segmentos da string. Para cada par de posições,
    inverte o segmento entre elas se isso melhorar o fitness.

    ESTRATÉGIA ALGORÍTMICA:
    - Enumera todos os pares de posições que definem segmentos
    - Para cada segmento (mínimo 2 caracteres), tenta a inversão
    - Avalia fitness após cada inversão
    - Mantém inversões que melhoram o fitness
    - Repete até não haver mais melhorias

    Args:
        ind (String): Indivíduo a ser refinado
        strings (list[String]): Lista de strings de referência para cálculo do fitness

    Returns:
        String: Indivíduo refinado através de inversões de segmentos

    Examples:
        >>> strings = ["ABCDEF", "ABFEDC", "ABEDCF"]
        >>> refine_2opt("ABCDEF", strings)
        'ABFEDC'  # Segmentos invertidos para melhor alinhamento

        >>> refine_2opt("ACGT", ["ACGT"])
        'ACGT'  # Já é ótimo, não há mudança

    Note:
        - Preserva a composição de símbolos da string
        - Complexidade O(n²) por iteração, pode ter múltiplas iterações
        - Útil para reorganizar estruturas invertidas
        - Pode encontrar soluções que outras buscas locais não conseguem
        - Uses standardized distance calculation when available
    """
    if max_distance_func is None:
        # Fallback for backward compatibility
        from src.domain.distance import HammingDistanceCalculator

        distance_calc = HammingDistanceCalculator(strings)
        max_distance_func = distance_calc.max_distance

    best_ind = ind
    best_fitness = max_distance_func(best_ind)
    current = list(best_ind)
    L = len(current)
    improved = True

    # Repete até não haver melhorias
    while improved:
        improved = False

        # Tenta inversões de segmentos
        for i in range(L):
            for j in range(i + 2, L + 1):  # Segmento de pelo menos 2 caracteres
                # Inverte segmento de i a j-1
                new_current = current[:i] + current[i:j][::-1] + current[j:]
                test_ind = "".join(new_current)
                test_fitness = max_distance_func(test_ind)

                # Se melhorou, atualiza e marca melhoria
                if test_fitness < best_fitness:
                    best_fitness = test_fitness
                    best_ind = test_ind
                    current = new_current
                    improved = True
                    break

            # Se encontrou melhoria, reinicia busca
            if improved:
                break

    return best_ind
