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
from utils.config import SYNTHETIC_DEFAULTS

logger = logging.getLogger(__name__)

def generate_dataset() -> Tuple[List[str], Dict[str, Any]]:
    """Interatividade completa via prompt."""
    defaults = SYNTHETIC_DEFAULTS
    
    n_input = input(f"Quantas strings (n)? [{defaults['n']}]: ").strip()
    n = int(n_input) if n_input else defaults['n']
    
    L_input = input(f"Comprimento das strings (L)? [{defaults['L']}]: ").strip()
    L = int(L_input) if L_input else defaults['L']
    
    alpha_input = input(f"Alfabeto (ex.: ACGT ou ABCDE...)? [{defaults['alphabet']}]: ").strip().upper()
    alphabet = alpha_input if alpha_input else defaults['alphabet']
    
    noise_input = input(f"Taxa de ruído por posição (0–1) [{defaults['noise']}]: ").strip()
    noise = float(noise_input) if noise_input else defaults['noise']
    
    params = {'n': n, 'L': L, 'alphabet': alphabet, 'noise': noise}
    logger.debug(f"Gerando dataset sintético com n={n}, L={L}, |Σ|={len(alphabet)}, noise={noise}")

    rng = random.Random(42)
    base = "".join(rng.choice(alphabet) for _ in range(L))
    logger.debug(f"String base: {base}")

    data = []
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
    logger.info(f"Dataset sintético gerado: n={n}, L={L}, |Σ|={len(alphabet)}")
    return data, params
