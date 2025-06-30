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

def generate_dataset_with_params(params: dict) -> Tuple[List[str], Dict[str, Any]]:
    """
    Gera dataset sintético com parâmetros específicos fornecidos.
    
    Args:
        params: Dicionário com parâmetros (n, L, alphabet, noise)
        
    Returns:
        Tupla (sequências, parâmetros_usados)
    """
    # Merge com defaults
    merged_params = {**SYNTHETIC_DEFAULTS}
    merged_params.update(params)
    
    n = merged_params['n']
    L = merged_params['L']
    alphabet = merged_params['alphabet']
    noise = merged_params['noise']
    
    console.print(f"Gerando dataset sintético: n={n}, L={L}, alphabet='{alphabet}', noise={noise}")
    
    # Usar seed fixa para reprodutibilidade em batch
    rng = random.Random(42)
    
    # Gerar string base
    base_string = ''.join(rng.choice(alphabet) for _ in range(L))
    
    # Gerar variações com ruído
    sequences = []
    for i in range(n):
        seq = list(base_string)
        num_mutations = int(L * noise)
        mutation_positions = rng.sample(range(L), min(num_mutations, L))
        
        for pos in mutation_positions:
            # Escolhe letra diferente da atual
            current_char = seq[pos]
            other_chars = [c for c in alphabet if c != current_char]
            if other_chars:
                seq[pos] = rng.choice(other_chars)
        
        sequences.append(''.join(seq))
    
    used_params = {
        'n': n,
        'L': L,
        'alphabet': alphabet,
        'noise': noise,
        'base_string': base_string
    }
    
    return sequences, used_params
