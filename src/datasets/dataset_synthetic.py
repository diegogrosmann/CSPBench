"""
dataset_synthetic.py
====================

Funções para gerar datasets sintéticos (aleatórios ou base + ruído).
Chamado pelo main.py quando o usuário opta por geração.
"""

import logging
import random
from typing import Any

from src.ui.cli.console_manager import console
from src.utils.config import SYNTHETIC_DEFAULTS, safe_input

logger = logging.getLogger(__name__)


def generate_dataset(silent: bool = False) -> tuple[list[str], dict[str, Any]]:
    """Interatividade completa via prompt."""
    defaults = SYNTHETIC_DEFAULTS
    if silent:
        n = defaults["n"]
        L = defaults["L"]
        alphabet = defaults["alphabet"]
        noise = defaults["noise"]
        fully_random = False
        seed = None
    else:
        n_input = safe_input(f"Quantas strings (n)? [{defaults['n']}]: ")
        n = int(n_input) if n_input else defaults["n"]

        L_input = safe_input(f"Comprimento das strings (L)? [{defaults['L']}]: ")
        L = int(L_input) if L_input else defaults["L"]

        alpha_input = safe_input(f"Alfabeto (ex.: ACGT ou ABCDE...)? [{defaults['alphabet']}]: ").upper()
        alphabet = alpha_input if alpha_input else defaults["alphabet"]

        noise_input = safe_input(f"Taxa de ruído por posição (0–1) [{defaults['noise']}]: ")
        noise = float(noise_input) if noise_input else defaults["noise"]

        fully_random_input = safe_input("Gerar strings totalmente aleatórias? (s/n) [n]: ").lower()
        fully_random = fully_random_input.startswith("s")
        # Pergunta por semente para reprodutibilidade (opcional)
        seed_input = safe_input("Semente aleatória (int) [None]: ")
        seed = int(seed_input) if seed_input else None

    # Gerar semente automaticamente se não fornecida
    if seed is None:
        import time

        seed = int(time.time() * 1000000) % (2**32)  # Gerar semente baseada no timestamp
        logger.info(f"Semente gerada automaticamente: {seed}")

    params = {
        "n": n,
        "L": L,
        "alphabet": alphabet,
        "noise": noise,
        "fully_random": fully_random,
    }
    params["seed"] = seed
    logger.info("Iniciando geração do dataset sintético")
    logger.debug(
        f"Gerando dataset sintético com n={n}, L={L}, |Σ|={len(alphabet)}, noise={noise}, fully_random={fully_random}, seed={seed}"
    )
    rng = random.Random(seed)
    data = []

    # Geração da string centro
    base_string = "".join(rng.choices(alphabet, k=L))
    logger.info("String centro gerada para aplicar ruído")

    # Definição do noise
    if noise is None:
        noise = rng.uniform(0.1, 0.5)
        logger.info(f"Noise não informado, sorteado aleatoriamente: {noise:.3f}")
    else:
        logger.info(f"Noise informado: {noise}")

    # Geração das strings com noise
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

    params = {
        "n": n,
        "L": L,
        "alphabet": alphabet,
        "fully_random": fully_random,
        "noise": noise,
        "seed": seed,
        "base_string": base_string,
    }
    logger.info(f"Parâmetros retornados pelo gerador: {params}")
    return data, params


def _unpack_singleton_lists(params: dict) -> dict:
    """
    Converte valores que são listas de um elemento para o valor simples.
    Exemplo: {'n': [80]} -> {'n': 80}
    """
    return {k: (v[0] if isinstance(v, list) and len(v) == 1 else v) for k, v in params.items()}


def generate_dataset_with_params(params: dict) -> tuple[list[str], dict[str, Any]]:
    """
    Gera dataset sintético com parâmetros específicos fornecidos.

    Args:
        params: Dicionário com parâmetros (n, L, alphabet, noise, fully_random)

    Returns:
        Tupla (sequências, parâmetros_usados)
    """
    # Desempacotar listas singleton
    params = _unpack_singleton_lists(params)
    # Merge com defaults
    merged_params = {**SYNTHETIC_DEFAULTS}
    merged_params.update(params)

    n = merged_params["n"]
    L = merged_params["L"]
    alphabet = merged_params["alphabet"]
    # Aceita tanto fully_random quanto aleatorio_total
    fully_random = merged_params.get("fully_random", merged_params.get("aleatorio_total", False))
    seed = merged_params.get("seed", None)

    # Gerar semente automaticamente se não fornecida
    if seed is None:
        import time

        seed = int(time.time() * 1000000) % (2**32)  # Gerar semente baseada no timestamp
        logger.info(f"Semente gerada automaticamente: {seed}")

    noise = merged_params.get("noise", None)  # Pode ser None agora

    console.print(f"Parâmetros usados: n={n}, L={L}, alphabet='{alphabet}', noise={noise}, fully_random={fully_random}")
    rng = random.Random(seed)
    sequences = []
    used_params = {"n": n, "L": L, "alphabet": alphabet, "fully_random": fully_random}

    # Geração da string centro
    base_string = "".join(rng.choices(alphabet, k=L))

    # Incluir a string base nos parâmetros usados
    used_params["base_string"] = base_string

    if fully_random:
        # Geração completamente aleatória
        for _ in range(n):
            seq = "".join(rng.choice(alphabet) for _ in range(L))
            sequences.append(seq)
        used_params["noise"] = None
    else:
        # Definição do noise
        if noise is not None:
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
            # Noise aleatório para cada string
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
    used_params["seed"] = seed

    # Calcular a distância da string base para o conjunto
    from src.utils.distance import max_distance

    distancia_base = max_distance(base_string, sequences)
    used_params["distancia_string_base"] = distancia_base

    console.print(f"Distância da string base: {used_params['distancia_string_base']}, Semente: {seed}")

    return sequences, used_params


def generate_perturbations(base_string: str, n_strings: int, max_distance: int, alphabet: str, rng) -> list[str]:
    """Gera strings com distância exata igual a max_distance da string base."""
    import logging

    logger = logging.getLogger(__name__)

    # Geração silenciosa de strings

    strings = []
    L = len(base_string)

    for i in range(n_strings):
        # Criar string com exatamente max_distance diferenças
        new_string = list(base_string)
        positions = rng.choice(L, size=max_distance, replace=False)

        logger.debug(f"[GENERATOR] String {i}: posições alteradas: {positions}")

        for pos in positions:
            # Escolher um símbolo diferente do atual
            current_char = base_string[pos]
            available_chars = [c for c in alphabet if c != current_char]
            if available_chars:  # Verificação de segurança
                chosen_char = rng.choice(available_chars)
                new_string[pos] = chosen_char
                logger.debug(f"[GENERATOR] String {i}: pos {pos}: {current_char} -> {chosen_char}")

        result_string = "".join(new_string)
        strings.append(result_string)

        # VALIDAÇÃO: Verificar se a distância está correta
        from src.utils.distance import hamming_distance

        actual_distance = hamming_distance(base_string, result_string)
        # String gerada silenciosamente

        if actual_distance != max_distance:
            logger.error(f"[GENERATOR] ERRO: String {i} tem distância {actual_distance}, esperado {max_distance}")

    # VALIDAÇÃO FINAL: Verificar todas as distâncias
    from src.utils.distance import hamming_distance

    all_distances = [hamming_distance(base_string, s) for s in strings]
    logger.info(f"[GENERATOR] Todas as distâncias geradas: {all_distances}")
    logger.info(f"[GENERATOR] Máxima: {max(all_distances)}, Mínima: {min(all_distances)}")

    return strings
