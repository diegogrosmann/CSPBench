from __future__ import annotations

"""
Implementação exata do DP-CSP (Programação Dinâmica) para o Closest String Problem.

Funções auxiliares:
    _dp_decision(...): DP decisório para existência de centro com raio <= d.
    exact_dp_closest_string(...): Busca incremental do menor raio possível.
"""

# dp_csp.py
"""
DP-CSP – Programação Dinâmica exata para o Closest String Problem
=================================================================

Implementa busca exata em O(L · |Σ| · (d+1)^n).  Útil para instâncias
pequenas (n ≤ 10, d ≤ 8).  Inclui busca incremental sobre d para achar
o menor raio possível até um limite superior fornecido.
"""

import logging
from collections.abc import Callable, Sequence
from typing import cast

from src.utils.distance import max_distance
from src.utils.resource_monitor import get_safe_memory_limit

from .config import DP_CSP_DEFAULTS

logger = logging.getLogger(__name__)

String = str
RemVec = tuple[int, ...]  # vetor de erros restantes p/ cada string


# ----------------------------------------------------------------------
# DP decisório: existe centro com raio ≤ d?
# ----------------------------------------------------------------------


def _dp_decision(strings: Sequence[String], alphabet: str, d: int) -> String | None:
    """
    Se existir string-centro c com raio ≤ d, devolve uma delas; caso
    contrário retorna None.
    """
    n, L = len(strings), len(strings[0])

    # Logs removidos para reduzir verbosidade

    # Pré-cálculo δσ[pos] = vetor 0/1 indicando discrepâncias
    delta: list[dict[String, RemVec]] = []
    for pos in range(L):
        col = [s[pos] for s in strings]
        delta.append({σ: tuple(int(σ != c) for c in col) for σ in alphabet})

    start: RemVec = (d,) * n
    frontier: set[RemVec] = {start}
    parent: dict[tuple[int, RemVec], tuple[RemVec | None, String]] = {
        (0, start): (None, "")  # (prev_rem, char usado para chegar aqui)
    }

    # Estado inicial silencioso

    for pos in range(L):
        # Processamento silencioso da posição
        nxt: set[RemVec] = set()
        for rem in frontier:
            # Estados processados silenciosamente
            for σ, dv in delta[pos].items():
                new_rem = tuple(r - v for r, v in zip(rem, dv))
                # Teste silencioso de caracteres
                if min(new_rem) < 0:
                    # Estado inviável - processamento silencioso
                    continue
                key = (pos + 1, new_rem)
                if key in parent:
                    # Estado já visitado - processamento silencioso
                    continue
                parent[key] = (rem, σ)
                nxt.add(new_rem)
                # Estado adicionado silenciosamente
        frontier = nxt
        # Nova fronteira processada silenciosamente
        if not frontier:  # nenhum estado viável
            # Solução inviável - retorno silencioso
            return None

    # Reconstrução -------------------------------------------------------
    final_rem = next(iter(frontier))
    # Reconstrução silenciosa da solução

    center_chars: list[String] = []
    pos: int = L
    rem: RemVec = final_rem  # rem nunca é None dentro do loop

    while pos > 0:
        prev_rem, σ = parent[(pos, rem)]
        center_chars.append(σ)
        # Caractere escolhido silenciosamente
        pos -= 1
        if prev_rem is None:  # chegamos ao estado inicial
            break
        rem = cast(RemVec, prev_rem)  # garante tipo para o próximo acesso

    center_chars.reverse()
    result = "".join(center_chars)
    # Solução reconstruída silenciosamente

    # Validação silenciosa
    from src.utils.distance import hamming_distance

    max_dist = max(hamming_distance(result, s) for s in strings)
    if max_dist > d:
        logger.error(f"[DP_DECISION] ERRO: Solução inválida! dist={max_dist} > d={d}")

    return result


# ----------------------------------------------------------------------
# Interface principal
# ----------------------------------------------------------------------


def exact_dp_closest_string(
    strings: list[String],
    alphabet: str,
    max_d: int | None = None,
    progress_callback: Callable[[str], None] | None = None,
    warning_callback: Callable[[str], None] | None = None,
) -> tuple[String, int]:
    """
    Busca o menor raio d* ≤ max_d tal que exista string-centro.
    Retorna (centro, d*).  Levanta RuntimeError se não encontrar.
    Monitora memória e tempo para evitar travamento.
    """
    import time

    baseline_val = max_distance(strings[0], strings)  # cota superior simples
    if max_d is None:
        max_d = baseline_val

    n = len(strings)

    # LOG DETALHADO: Parâmetros de entrada
    logger.info(f"[DP_CSP] Iniciando busca exata com max_d={max_d}, baseline={baseline_val}")
    logger.info(f"[DP_CSP] Dataset: n={n}, L={len(strings[0])}, alfabeto={alphabet}")
    for i, s in enumerate(strings):
        logger.info(f"[DP_CSP] String {i}: {s}")

    # --- Monitoramento de recursos ---
    safe_mem_mb = get_safe_memory_limit()
    max_time = DP_CSP_DEFAULTS.get("max_time", 300)
    t0 = time.time()

    logger.info(f"[DP_CSP] Limites: mem={safe_mem_mb:.1f}MB, tempo={max_time}s")

    def check_limits(d):
        # Checa memória, tempo e viabilidade incrementalmente
        import gc
        import os

        gc.collect()
        mem_mb = 0.0
        try:
            if os.path.exists("/proc/self/status"):
                with open("/proc/self/status") as f:
                    for line in f:
                        if line.startswith("VmRSS:"):
                            kb = int(line.split()[1])
                            mem_mb = kb / 1024.0
                            break
        except Exception:
            pass
        elapsed = time.time() - t0
        # Limite prático para (d+1)^n: ~2e9 (16GB RAM, 8 bytes/estado)
        state_count_est = (d + 1) ** n
        # Verificação silenciosa de recursos

        if state_count_est > 2_000_000_000:
            msg = (
                f"DP-CSP interrompido: (d+1)^n = {state_count_est:,} excede limite prático "
                "(~2e9 estados, ~16GB RAM). Tente reduzir n ou d."
            )
            logger.error(f"[DP_CSP] {msg}")
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)
        if mem_mb > safe_mem_mb * 0.95:
            msg = f"DP-CSP interrompido: uso de memória {mem_mb:.1f}MB excedeu limite seguro ({safe_mem_mb:.1f}MB)"
            logger.error(f"[DP_CSP] {msg}")
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)
        if elapsed > max_time:
            msg = f"DP-CSP interrompido: tempo de execução excedeu {max_time}s"
            logger.error(f"[DP_CSP] {msg}")
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

    for d in range(max_d + 1):
        # Checar limites antes de cada iteração principal
        check_limits(d)
        if progress_callback:
            progress_callback(f"Testando d={d}")
        logger.info(f"[DP_CSP] Testando d={d} (tentativa {d+1}/{max_d+1})")

        center = _dp_decision(strings, alphabet, d)
        if center is not None:
            # VALIDAÇÃO FINAL: Verificar se a solução está correta
            from src.utils.distance import hamming_distance

            max_dist = max(hamming_distance(center, s) for s in strings)
            logger.info(f"[DP_CSP] SUCESSO! Encontrou solução com d={d}")
            logger.info(f"[DP_CSP] Centro encontrado: {center}")
            logger.info(f"[DP_CSP] Validação final: distância máxima = {max_dist}")

            if max_dist != d:
                logger.warning(f"[DP_CSP] INCONSISTÊNCIA: d={d} mas distância real={max_dist}")

            return center, d
        else:
            # Nenhuma solução encontrada - processamento silencioso
            pass

    logger.error(f"[DP_CSP] FALHA: Não foi possível encontrar centro com d ≤ {max_d}")
    raise RuntimeError(f"Não foi possível encontrar centro com d ≤ {max_d}. " "Tente aumentar o limite.")
