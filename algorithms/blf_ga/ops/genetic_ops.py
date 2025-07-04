"""
Operações genéticas para BLF-GA: mutação, crossover, refinamento, diversidade, etc.
"""

import random

import numpy as np

String = str
Population = list[String]


# --- Diversidade ---
def mean_hamming_distance(pop: Population) -> float:
    """Calcula a média das distâncias de Hamming entre todos os pares da população."""
    if len(pop) < 2:
        return 0.0
    arr = np.array([list(ind) for ind in pop])
    dists = np.sum(arr[:, None, :] != arr[None, :, :], axis=2)
    iu = np.triu_indices(len(pop), 1)
    return np.mean(dists[iu])


# --- Mutação ---
def mutate_multi(ind: str, alphabet: str, rng: random.Random, n: int = 2) -> str:
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
    chars = list(ind)
    L = len(chars)
    a, b = sorted(rng.sample(range(L), 2))
    chars[a : b + 1] = chars[a : b + 1][::-1]
    return "".join(chars)


def mutate_transposition(ind: str, rng: random.Random) -> str:
    chars = list(ind)
    L = len(chars)
    a, b = sorted(rng.sample(range(L), 2))
    seg = chars[a : b + 1]
    del chars[a : b + 1]
    pos = rng.randint(0, len(chars))
    chars[pos:pos] = seg
    return "".join(chars)


# --- Crossover ---
def crossover_one_point(p1: String, p2: String, rng: random.Random) -> tuple[String, String]:
    L = len(p1)
    point = rng.randint(1, L - 1)
    c1 = p1[:point] + p2[point:]
    c2 = p2[:point] + p1[point:]
    return c1, c2


def crossover_uniform(p1: String, p2: String, rng: random.Random) -> tuple[String, String]:
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
    c1 = list(p1)
    c2 = list(p2)
    for l, r in blocks:
        if rng.random() < 0.5:
            c1[l:r] = p2[l:r]
            c2[l:r] = p1[l:r]
    return "".join(c1), "".join(c2)


# --- Refinamento (placeholders, implementar depois) ---
def refine_swap(ind: String, strings: list[String]) -> String:
    # Implementar swap refinement
    return ind


def refine_insertion(ind: String, strings: list[String]) -> String:
    # Implementar insertion refinement
    return ind


def refine_2opt(ind: String, strings: list[String]) -> String:
    # Implementar 2-opt refinement
    return ind
