"""
Operações genéticas para BLF-GA: mutação, crossover, refinamento, diversidade, etc.

Este módulo implementa operadores genéticos utilizados em algoritmos genéticos aplicados ao problema BLF-GA (Block Letter Frequency Genetic Algorithm).
Inclui funções para:
- Medir diversidade populacional (distância de Hamming)
- Realizar diferentes tipos de mutação (multi-ponto, inversão, transposição)
- Realizar diferentes tipos de crossover (um ponto, uniforme, por blocos)
- Refinamento local de indivíduos (busca gulosa, trocas, inserção, 2-opt)

Cada função está documentada individualmente. O objetivo é fornecer ferramentas flexíveis para manipulação de populações de strings, onde cada string representa um indivíduo.
"""

import random

import numpy as np

String = str
Population = list[String]


# --- Diversidade ---
def mean_hamming_distance(pop: Population) -> float:
    """
    Calcula a média das distâncias de Hamming entre todos os pares da população.
    Útil para medir a diversidade genética da população.

    Args:
        pop: Lista de indivíduos (strings)
    Returns:
        float: Média das distâncias de Hamming
    """
    if len(pop) < 2:
        return 0.0
    arr = np.array([list(ind) for ind in pop])
    dists = np.sum(arr[:, None, :] != arr[None, :, :], axis=2)
    iu = np.triu_indices(len(pop), 1)
    return np.mean(dists[iu])


# --- Mutação ---
def mutate_multi(ind: str, alphabet: str, rng: random.Random, n: int = 2) -> str:
    """
    Realiza mutação em múltiplos pontos do indivíduo, alterando até n posições para símbolos diferentes do alfabeto.

    Args:
        ind: String original (indivíduo)
        alphabet: Alfabeto permitido
        rng: Instância de random.Random
        n: Número de mutações
    Returns:
        String mutada
    """
    chars = list(ind)
    L = len(chars)
    for _ in range(n):
        pos = rng.randint(0, L - 1)
        old = chars[pos]
        choices = [c for c in alphabet if c != old]
        if choices:
            chars[pos] = rng.choice(choices)
    return "".join(chars)


def mutate_inversion(ind: str, rng: random.Random) -> str:
    """
    Realiza mutação por inversão: seleciona dois pontos e inverte o segmento entre eles.
    """
    chars = list(ind)
    L = len(chars)
    a, b = sorted(rng.sample(range(L), 2))
    chars[a : b + 1] = chars[a : b + 1][::-1]
    return "".join(chars)


def mutate_transposition(ind: str, rng: random.Random) -> str:
    """
    Realiza mutação por transposição: seleciona um segmento e o move para outra posição.
    """
    chars = list(ind)
    L = len(chars)
    a, b = sorted(rng.sample(range(L), 2))
    seg = chars[a : b + 1]
    del chars[a : b + 1]
    pos = rng.randint(0, len(chars))
    chars[pos:pos] = seg
    return "".join(chars)


# --- Crossover ---
def crossover_one_point(
    p1: String, p2: String, rng: random.Random
) -> tuple[String, String]:
    """
    Crossover de um ponto: corta os pais em um ponto aleatório e troca os sufixos.
    """
    L = len(p1)
    point = rng.randint(1, L - 1)
    c1 = p1[:point] + p2[point:]
    c2 = p2[:point] + p1[point:]
    return c1, c2


def crossover_uniform(
    p1: String, p2: String, rng: random.Random
) -> tuple[String, String]:
    """
    Crossover uniforme: para cada posição, escolhe aleatoriamente o gene de um dos pais.
    """
    L = len(p1)
    c1 = []
    c2 = []
    for i in range(L):
        if rng.random() < 0.5:
            c1.append(p1[i])
            c2.append(p2[i])
        else:
            c1.append(p2[i])
            c2.append(p1[i])
    return "".join(c1), "".join(c2)


# Blending dos blocos pode ser implementado conforme necessário


def crossover_blend_blocks(
    p1: String, p2: String, blocks: list[tuple[int, int]], rng: random.Random
) -> tuple[String, String]:
    """
    Crossover por blocos: para cada bloco definido, troca segmentos entre os pais com 50% de chance.
    Args:
        p1, p2: Strings dos pais
        blocks: Lista de tuplas (início, fim) dos blocos
        rng: Random
    Returns:
        Dois filhos resultantes
    """
    c1 = list(p1)
    c2 = list(p2)
    for l, r in blocks:
        if rng.random() < 0.5:
            c1[l:r] = p2[l:r]
            c2[l:r] = p1[l:r]
    return "".join(c1), "".join(c2)


# --- Refinamento Local ---
def refine_greedy(ind: String, strings: list[String]) -> String:
    """
    Refinamento guloso posição por posição.

    Para cada posição, tenta todos os símbolos do alfabeto e
    escolhe o que minimiza a distância máxima para as strings de referência.
    Repete até não haver mais melhorias.
    """
    from collections import Counter

    from src.utils.distance import max_distance

    # Extrai alfabeto das strings
    alphabet = set()
    for s in strings:
        alphabet.update(s)
    alphabet = sorted(alphabet)
    best_ind = ind
    best_fitness = max_distance(best_ind, strings)
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
                    test_fitness = max_distance(test_ind, strings)
                    if test_fitness < best_fitness:
                        best_fitness = test_fitness
                        best_char = char
                        improved = True
            # Aplica a melhor mudança
            current[pos] = best_char
        best_ind = "".join(current)
    return best_ind


def refine_swap(ind: String, strings: list[String]) -> String:
    """
    Refinamento por troca de posições.

    Tenta trocar pares de posições para melhorar o fitness.

    Args:
        ind: Indivíduo a ser refinado
        strings: Strings de referência

    Returns:
        String: Indivíduo refinado
    """
    from src.utils.distance import max_distance

    best_ind = ind
    best_fitness = max_distance(best_ind, strings)

    current = list(best_ind)
    L = len(current)

    # Tenta todas as trocas de pares
    for i in range(L):
        for j in range(i + 1, L):
            # Troca posições i e j
            current[i], current[j] = current[j], current[i]
            test_ind = "".join(current)
            test_fitness = max_distance(test_ind, strings)

            if test_fitness < best_fitness:
                best_fitness = test_fitness
                best_ind = test_ind
            else:
                # Desfaz a troca se não melhorou
                current[i], current[j] = current[j], current[i]

    return best_ind


def refine_insertion(ind: String, strings: list[String]) -> String:
    """
    Refinamento por inserção/remoção de segmentos.

    Tenta mover segmentos pequenos para diferentes posições.

    Args:
        ind: Indivíduo a ser refinado
        strings: Strings de referência

    Returns:
        String: Indivíduo refinado
    """
    from src.utils.distance import max_distance

    best_ind = ind
    best_fitness = max_distance(best_ind, strings)

    current = list(best_ind)
    L = len(current)

    # Tenta mover segmentos de tamanho 1-3
    for seg_len in range(1, min(4, L)):
        for start in range(L - seg_len + 1):
            segment = current[start : start + seg_len]

            # Remove segmento
            temp = current[:start] + current[start + seg_len :]

            # Tenta inserir em cada posição
            for insert_pos in range(len(temp) + 1):
                new_current = temp[:insert_pos] + segment + temp[insert_pos:]

                if len(new_current) == L:  # Mantém comprimento
                    test_ind = "".join(new_current)
                    test_fitness = max_distance(test_ind, strings)

                    if test_fitness < best_fitness:
                        best_fitness = test_fitness
                        best_ind = test_ind
                        current = new_current

    return best_ind


def refine_2opt(ind: String, strings: list[String]) -> String:
    """
    Refinamento 2-opt: inverte segmentos da string.

    Para cada par de posições, inverte o segmento entre elas
    se isso melhorar o fitness.

    Args:
        ind: Indivíduo a ser refinado
        strings: Strings de referência

    Returns:
        String: Indivíduo refinado
    """
    from src.utils.distance import max_distance

    best_ind = ind
    best_fitness = max_distance(best_ind, strings)

    current = list(best_ind)
    L = len(current)
    improved = True

    while improved:
        improved = False

        # Tenta inversões de segmentos
        for i in range(L):
            for j in range(i + 2, L + 1):  # Segmento de pelo menos 2 caracteres
                # Inverte segmento de i a j-1
                new_current = current[:i] + current[i:j][::-1] + current[j:]
                test_ind = "".join(new_current)
                test_fitness = max_distance(test_ind, strings)

                if test_fitness < best_fitness:
                    best_fitness = test_fitness
                    best_ind = test_ind
                    current = new_current
                    improved = True
                    break

            if improved:
                break

    return best_ind
