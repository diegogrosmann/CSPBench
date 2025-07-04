"""
Implementação do H³-CSP (Hybrid Hierarchical Hamming Search) para CSP.

Classes:
    H3CSP: Implementa o algoritmo H³-CSP.

Funções auxiliares:
    split_in_blocks(L): Divide posições em blocos.
    consensus_block(strings, l, r): Consenso por bloco.
    ... (demais funções auxiliares)
"""

from __future__ import annotations

import itertools
import logging
import math
import random
import time
from collections import Counter
from collections.abc import Callable, Sequence

from csp_blfga.utils.distance import max_distance

from .config import H3_CSP_DEFAULTS

logger = logging.getLogger(__name__)

String = str
Block = tuple[int, int]  # (início, fim) 0-based, exclusivo


# ---------------------------------------------------------------------------
# Auxiliares de bloco
# ---------------------------------------------------------------------------


def split_in_blocks(L: int) -> list[Block]:
    """Divide posições 0..L-1 em B ≈ ⌈√L⌉ blocos contíguos."""
    B = math.ceil(math.sqrt(L))
    base_size = math.ceil(L / B)
    blocks = []
    cur = 0
    while cur < L:
        blocks.append((cur, min(cur + base_size, L)))
        cur += base_size
    return blocks


def consensus_block(strings: Sequence[String], l: int, r: int) -> String:
    """Consenso (maioria) para o intervalo [l:r)."""
    counter_cols = [Counter(s[l:r]) for s in strings]
    rep = "".join(Counter(col).most_common(1)[0][0] for col in zip(*[s[l:r] for s in strings]))
    return rep


# ---------------------------------------------------------------------------
# Técnicas por bloco (versão simplificada)
# ---------------------------------------------------------------------------


def _exhaustive_block(strings: Sequence[String], alphabet: str, l: int, r: int, k: int) -> list[String]:
    """
    Busca exaustiva dentro do intervalo se |Σ|^(r-l) < 10 000.
    Caso contrário, retorna os k melhores blocos existentes + consenso.
    """
    m = r - l
    lim = len(alphabet) ** m
    bests: list[tuple[int, String]] = []

    if lim <= 10_000:
        # Bloco processado silenciosamente
        for cand_tuple in itertools.product(alphabet, repeat=m):
            cand = "".join(cand_tuple)
            dist = max_distance(cand, [s[l:r] for s in strings])
            if len(bests) < k or dist < bests[-1][0]:
                bests.append((dist, cand))
                bests.sort(key=lambda x: x[0])
                bests = bests[:k]
    else:
        # Limita-se a candidatos vindos do próprio dataset
        # Bloco grande processado silenciosamente
        seen: set[String] = set()
        for s in strings:
            blk = s[l:r]
            if blk not in seen:
                seen.add(blk)
                dist = max_distance(blk, [t[l:r] for t in strings])
                bests.append((dist, blk))
        bests.sort(key=lambda x: x[0])
        bests = bests[:k]

    # Garante inclusão do consenso
    cons = consensus_block(strings, l, r)
    bests.append((max_distance(cons, [s[l:r] for s in strings]), cons))
    bests.sort(key=lambda x: x[0])
    return [c for _, c in bests[:k]]


def _beam_search_block(
    strings: Sequence[String], alphabet: str, l: int, r: int, beam_width: int, k: int
) -> list[String]:
    """
    Beam search simples posição-a-posição para gerar candidatos de bloco.
    """
    m = r - l
    beam = [""]  # prefixos
    for pos in range(m):
        scored: list[tuple[int, String]] = []
        for prefix in beam:
            for a in alphabet:
                cand = prefix + a
                # avalia parcialmente: dist máx considerando prefixo
                partial_dists = []
                for s in strings:
                    mismatch = sum(1 for i, c in enumerate(cand) if i < len(cand) and c != s[l + i])
                    partial_dists.append(mismatch)
                scored.append((max(partial_dists), cand))
        scored.sort(key=lambda x: x[0])
        beam = [s for _, s in scored[:beam_width]]

    # Score final em todo bloco
    final_scored = [(max_distance(c, [s[l:r] for s in strings]), c) for c in beam]
    final_scored.sort(key=lambda x: x[0])
    return [c for _, c in final_scored[:k]]


# ---------------------------------------------------------------------------
# Busca local global (hill-climbing)
# ---------------------------------------------------------------------------


def _local_search(candidate: String, strings: Sequence[String]) -> String:
    """
    Hill-climbing simples, tentando mudar cada posição para qualquer letra
    presente nas strings originais.
    """
    cand_list = list(candidate)
    L = len(cand_list)
    alphabet_by_pos = [{s[i] for s in strings} for i in range(L)]

    improved = True
    while improved:
        improved = False
        base_val = max_distance("".join(cand_list), list(strings))
        for i in range(L):
            cur_char = cand_list[i]
            for alt in alphabet_by_pos[i]:
                if alt == cur_char:
                    continue
                cand_list[i] = alt
                new_val = max_distance("".join(cand_list), list(strings))
                if new_val < base_val:
                    base_val = new_val
                    improved = True
                    cur_char = alt
                else:
                    cand_list[i] = cur_char  # desfaz
    return "".join(cand_list)


# ---------------------------------------------------------------------------
# Classe principal
# ---------------------------------------------------------------------------


class H3CSP:
    """Interface compatível com BLFGA.run() → retorna (centro, distância)."""

    def __init__(self, strings: Sequence[String], alphabet: str, **params):
        self.strings = list(strings)
        self.alphabet = alphabet
        self.L = len(strings[0])

        # Merge defaults + params
        self.params = dict(H3_CSP_DEFAULTS)
        self.params.update(params)

        self.rng = random.Random(self.params["seed"])
        self.progress_callback: Callable[[str], None] | None = None

        self.blocks = split_in_blocks(self.L)
        # H3CSP inicializado silenciosamente

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Define um callback para relatar o progresso."""
        self.progress_callback = callback

    # ---------------------------------------------------------------------

    def _smart_core(self) -> list[list[String]]:
        """
        Para cada bloco devolve até k candidatos (strings do próprio bloco).
        """
        k = self.params["k_candidates"]
        small_limit = self.params["block_small"]
        medium_limit = self.params["block_medium"]
        beam_width = self.params["beam_width"]

        block_cands: list[list[String]] = []
        for l, r in self.blocks:
            # calcula d_b
            cons = consensus_block(self.strings, l, r)
            d_b = max_distance(cons, [s[l:r] for s in self.strings])

            if d_b <= small_limit:
                cands = _exhaustive_block(self.strings, self.alphabet, l, r, k)
            elif d_b <= medium_limit:
                # para já, beam search médio
                cands = _beam_search_block(self.strings, self.alphabet, l, r, beam_width // 2, k)
            else:
                cands = _beam_search_block(self.strings, self.alphabet, l, r, beam_width, k)

            block_cands.append(cands)
            # Processamento silencioso de candidatos
        return block_cands

    # ---------------------------------------------------------------------

    def _fuse_blocks(self, chosen: list[String]) -> String:
        """Concatena blocos na ordem da lista self.blocks."""
        return "".join(chosen)

    # ---------------------------------------------------------------------

    def run(self) -> tuple[String, int]:
        start_time = time.time()

        try:
            if self.progress_callback:
                self.progress_callback("Analisando blocos...")
            block_cands = self._smart_core()

            if self.progress_callback:
                self.progress_callback("Fusão de blocos...")
            best_by_block = [cands[0] for cands in block_cands]
            center = self._fuse_blocks(best_by_block)
            best_val = max_distance(center, list(self.strings))

            logger.info(f"[Fusão inicial] dist={best_val}")

            if self.progress_callback:
                self.progress_callback("Refinamento global...")
            for it in range(self.params["local_iters"]):
                # Verificar se callback foi chamado (indica possível cancelamento)
                if self.progress_callback:
                    self.progress_callback(f"Refinamento: iteração {it+1}")

                center = _local_search(center, self.strings)
                new_val = max_distance(center, list(self.strings))
                if new_val < best_val:
                    logger.info(f"[Local] it={it+1}  {best_val}->{new_val}")
                    best_val = new_val
                    if self.progress_callback:
                        self.progress_callback(f"Melhoria encontrada: distância={best_val}")
                else:
                    break

                if time.time() - start_time >= self.params["max_time"]:
                    logger.warning("Tempo máximo atingido no H3CSP")
                    if self.progress_callback:
                        self.progress_callback("Timeout atingido")
                    break

            return center, best_val

        except Exception as e:
            logger.error(f"Erro no H3CSP: {e}")
            if self.progress_callback:
                self.progress_callback(f"Erro: {str(e)}")
            raise e  # Re-raise para ser capturado no wrapper
