# dataset_synthetic.py
"""
dataset_synthetic.py
====================

Funções para gerar datasets sintéticos (aleatórios ou base + ruído).
Chamado pelo main.py quando o usuário opta por geração.
"""

import random
from typing import List, Tuple, Dict, Any
import logging
from utils.config import SYNTHETIC_DEFAULTS, safe_input
from src.console_manager import console

logger = logging.getLogger(__name__)

def generate_dataset() -> Tuple[List[str], Dict[str, Any]]:
    """Interatividade completa via prompt."""
    defaults = SYNTHETIC_DEFAULTS
    
    n_input = safe_input(f"Quantas strings (n)? [{defaults['n']}]: ")
    n = int(n_input) if n_input else defaults['n']
    
    L_input = safe_input(f"Comprimento das strings (L)? [{defaults['L']}]: ")
    L = int(L_input) if L_input else defaults['L']
    
    alpha_input = safe_input(f"Alfabeto (ex.: ACGT ou ABCDE...)? [{defaults['alphabet']}]: ").upper()
    alphabet = alpha_input if alpha_input else defaults['alphabet']
    
    noise_input = safe_input(f"Taxa de ruído por posição (0–1) [{defaults['noise']}]: ")
    noise = float(noise_input) if noise_input else defaults['noise']
    
    fully_random_input = safe_input(f"Gerar strings totalmente aleatórias? (s/n) [n]: ").lower()
    fully_random = fully_random_input.startswith('s')
    # Pergunta por semente para reprodutibilidade (opcional)
    seed_input = safe_input(f"Semente aleatória (int) [None]: ")
    seed = int(seed_input) if seed_input else None
    params = {'n': n, 'L': L, 'alphabet': alphabet, 'noise': noise, 'fully_random': fully_random}
    params['seed'] = seed
    logger.debug(f"Gerando dataset sintético com n={n}, L={L}, |Σ|={len(alphabet)}, noise={noise}, fully_random={fully_random}, seed={seed}")
    rng = random.Random(seed)
    data = []
    
    if fully_random:
        # Geração completamente aleatória
        for idx in range(n):
            new_s = "".join(rng.choice(alphabet) for _ in range(L))
            data.append(new_s)
            if idx < 3:
                logger.debug(f"String aleatória {idx}: {new_s}")
    else:
        # Geração baseada em uma string + ruído
        base = "".join(rng.choice(alphabet) for _ in range(L))
        logger.debug(f"String base: {base}")

        for idx in range(n):
            s = list(base)
            for i in range(L):
                if rng.random() < noise:
                    old = s[i]
                    s[i] = rng.choice([c for c in alphabet if c != old])
            new_s = "".join(s)
            data.append(new_s)
            if idx < 3:
                logger.debug(f"String gerada {idx}: {new_s}")
                
    logger.info(f"Dataset sintético gerado: n={n}, L={L}, |Σ|={len(alphabet)}, geração {'aleatória' if fully_random else 'base+ruído'}")
    return data, params

def generate_dataset_with_params(params: dict) -> Tuple[List[str], Dict[str, Any]]:
    """
    Gera dataset sintético com parâmetros específicos fornecidos.
    
    Args:
        params: Dicionário com parâmetros (n, L, alphabet, noise, fully_random)
        
    Returns:
        Tupla (sequências, parâmetros_usados)
    """
    # Merge com defaults
    merged_params = {**SYNTHETIC_DEFAULTS}
    merged_params.update(params)

    n = merged_params['n']
    L = merged_params['L']
    alphabet = merged_params['alphabet']
    # Aceita tanto fully_random quanto aleatorio_total
    fully_random = merged_params.get('fully_random', merged_params.get('aleatorio_total', False))
    seed = merged_params.get('seed', None)
    noise = merged_params.get('noise', None)  # Pode ser None agora

    console.print(f"Gerando dataset sintético: n={n}, L={L}, alphabet='{alphabet}', noise={noise}, fully_random={fully_random}")

    rng = random.Random(seed)
    sequences = []
    used_params = {
        'n': n,
        'L': L,
        'alphabet': alphabet,
        'fully_random': fully_random
    }

    if fully_random:
        # Geração completamente aleatória
        for _ in range(n):
            seq = ''.join(rng.choice(alphabet) for _ in range(L))
            sequences.append(seq)
        used_params['noise'] = None
    else:
        # Geração baseada em string base + ruído
        base_string = ''.join(rng.choice(alphabet) for _ in range(L))
        used_params['base_string'] = base_string

        if noise is not None:
            used_params['noise'] = noise
            for i in range(n):
                seq = list(base_string)
                for pos in range(L):
                    if rng.random() < noise:
                        current_char = seq[pos]
                        other_chars = [c for c in alphabet if c != current_char]
                        if other_chars:
                            seq[pos] = rng.choice(other_chars)
                sequences.append(''.join(seq))
        else:
            # Noise aleatório para cada string
            noises = [rng.uniform(0.01, 0.5) for _ in range(n)]
            used_params['noise'] = noises
            for i in range(n):
                seq = list(base_string)
                this_noise = noises[i]
                for pos in range(L):
                    if rng.random() < this_noise:
                        current_char = seq[pos]
                        other_chars = [c for c in alphabet if c != current_char]
                        if other_chars:
                            seq[pos] = rng.choice(other_chars)
                sequences.append(''.join(seq))

    return sequences, used_params
