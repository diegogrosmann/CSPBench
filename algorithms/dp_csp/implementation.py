# dp_csp.py
"""
DP-CSP – Programação Dinâmica exata para o Closest String Problem
=================================================================

Implementa busca exata em O(L · |Σ| · (d+1)^n).  Útil para instâncias
pequenas (n ≤ 10, d ≤ 8).  Inclui busca incremental sobre d para achar
o menor raio possível até um limite superior fornecido.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Sequence, Tuple, cast, Callable

from utils.distance import max_hamming
from utils.resource_monitor import get_safe_memory_limit, force_garbage_collection
from .config import DP_CSP_DEFAULTS

logger = logging.getLogger(__name__)

String = str
RemVec = Tuple[int, ...]          # vetor de erros restantes p/ cada string


# ----------------------------------------------------------------------
# DP decisório: existe centro com raio ≤ d?
# ----------------------------------------------------------------------

def _dp_decision(strings: Sequence[String],
                 alphabet: str,
                 d: int) -> Optional[String]:
    """
    Se existir string-centro c com raio ≤ d, devolve uma delas; caso
    contrário retorna None.
    """
    n, L = len(strings), len(strings[0])

    # Pré-cálculo δσ[pos] = vetor 0/1 indicando discrepâncias
    delta: List[Dict[String, RemVec]] = []
    for pos in range(L):
        col = [s[pos] for s in strings]
        delta.append({σ: tuple(int(σ != c) for c in col) for σ in alphabet})

    start: RemVec = (d,) * n
    frontier: set[RemVec] = {start}
    parent: Dict[Tuple[int, RemVec],
                 Tuple[Optional[RemVec], String]] = {
        (0, start): (None, '')       # (prev_rem, char usado para chegar aqui)
    }

    for pos in range(L):
        nxt: set[RemVec] = set()
        for rem in frontier:
            for σ, dv in delta[pos].items():
                new_rem = tuple(r - v for r, v in zip(rem, dv))
                if min(new_rem) < 0:
                    continue
                key = (pos + 1, new_rem)
                if key in parent:
                    continue
                parent[key] = (rem, σ)
                nxt.add(new_rem)
        frontier = nxt
        if not frontier:              # nenhum estado viável
            return None

    # Reconstrução -------------------------------------------------------
    final_rem = next(iter(frontier))
    center_chars: List[String] = []
    pos: int = L
    rem: RemVec = final_rem         # rem nunca é None dentro do loop

    while pos > 0:
        prev_rem, σ = parent[(pos, rem)]
        center_chars.append(σ)
        pos -= 1
        if prev_rem is None:         # chegamos ao estado inicial
            break
        rem = cast(RemVec, prev_rem)  # garante tipo para o próximo acesso

    center_chars.reverse()
    return ''.join(center_chars)


# ----------------------------------------------------------------------
# Interface principal
# ----------------------------------------------------------------------

def exact_dp_closest_string(strings: List[String],
                            alphabet: str,
                            max_d: Optional[int] = None,
                            progress_callback: Optional[Callable[[str], None]] = None,
                            warning_callback: Optional[Callable[[str], None]] = None
                            ) -> Tuple[String, int]:
    """
    Busca o menor raio d* ≤ max_d tal que exista string-centro.
    Retorna (centro, d*).  Levanta RuntimeError se não encontrar.
    Monitora memória e tempo para evitar travamento.
    """
    import time
    baseline_val = max_hamming(strings[0], strings)  # cota superior simples
    if max_d is None:
        max_d = baseline_val

    n = len(strings)

    # --- Monitoramento de recursos ---
    safe_mem_mb = get_safe_memory_limit()
    max_time = DP_CSP_DEFAULTS.get('max_time', 300)
    t0 = time.time()

    def check_limits(d):
        # Checa memória, tempo e viabilidade incrementalmente
        import os, gc
        gc.collect()
        mem_mb = 0.0
        try:
            if os.path.exists('/proc/self/status'):
                with open('/proc/self/status', 'r') as f:
                    for line in f:
                        if line.startswith('VmRSS:'):
                            kb = int(line.split()[1])
                            mem_mb = kb / 1024.0
                            break
        except Exception:
            pass
        elapsed = time.time() - t0
        # Limite prático para (d+1)^n: ~2e9 (16GB RAM, 8 bytes/estado)
        state_count_est = (d + 1) ** n
        if state_count_est > 2_000_000_000:
            msg = (f"DP-CSP interrompido: (d+1)^n = {state_count_est:,} excede limite prático "
                   "(~2e9 estados, ~16GB RAM). Tente reduzir n ou d.")
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)
        if mem_mb > safe_mem_mb * 0.95:
            msg = f"DP-CSP interrompido: uso de memória {mem_mb:.1f}MB excedeu limite seguro ({safe_mem_mb:.1f}MB)"
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)
        if elapsed > max_time:
            msg = f"DP-CSP interrompido: tempo de execução excedeu {max_time}s"
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

    for d in range(max_d + 1):
        # Checar limites antes de cada iteração principal
        check_limits(d)
        if progress_callback:
            progress_callback(f"Testando d={d}")
        else:
            logger.debug(f"DP-CSP tentando d={d}")
        center = _dp_decision(strings, alphabet, d)
        if center is not None:
            logger.info(f"DP-CSP encontrou solução com d={d}")
            return center, d

    raise RuntimeError(
        f"Não foi possível encontrar centro com d ≤ {max_d}. "
        "Tente aumentar o limite."
    )
