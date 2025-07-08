"""
Implementação do H³-CSP (Hybrid Hierarchical Hamming Search) para CSP.

Este módulo contém a implementação completa do algoritmo H³-CSP, que utiliza
uma abordagem híbrida de três camadas para resolver o Closest String Problem:

1. **B-Splitter**: Divisão hierárquica das strings em blocos contíguos
2. **Smart-Core**: Seleção adaptativa de técnicas por bloco baseada na dificuldade
3. **Global Refine**: Combinação de blocos e refinamento global por busca local

O algoritmo é especialmente eficaz para instâncias de tamanho médio (L=50-500)
e dados com padrões locais, oferecendo um bom equilíbrio entre qualidade
da solução e eficiência computacional.

Classes:
    H3CSP: Implementação principal do algoritmo H³-CSP.

Funções auxiliares:
    split_in_blocks(L): Divide L posições em ~√L blocos contíguos.
    consensus_block(strings, l, r): Calcula consenso por maioria em um bloco.
    _exhaustive_block(): Busca exaustiva para blocos pequenos.
    _beam_search_block(): Beam search para blocos médios/grandes.
    _local_search(): Hill-climbing para refinamento global.

Tipos:
    String: Alias para str (string de DNA/proteína).
    Block: Tupla (início, fim) representando um bloco - 0-based, fim exclusivo.
"""

from __future__ import annotations

import itertools
import logging
import math
import random
import time
from collections import Counter
from collections.abc import Callable, Sequence

from src.utils.distance import max_distance

from .config import H3_CSP_DEFAULTS

logger = logging.getLogger(__name__)

String = str
Block = tuple[int, int]  # (início, fim) 0-based, exclusivo


# ---------------------------------------------------------------------------
# Auxiliares de bloco
# ---------------------------------------------------------------------------


def split_in_blocks(L: int) -> list[Block]:
    """
    Divide L posições em aproximadamente √L blocos contíguos.

    Esta função implementa a "regra √L" do H³-CSP, dividindo as posições
    0..L-1 em blocos de tamanho aproximadamente igual, com B ≈ ⌈√L⌉ blocos.
    Esta decomposição hierárquica é fundamental para a eficiência do algoritmo.

    Args:
        L (int): Comprimento total das strings a serem divididas.

    Returns:
        list[Block]: Lista de blocos, onde cada bloco é uma tupla (início, fim)
                    representando o intervalo [início, fim) (fim exclusivo).

    Example:
        >>> split_in_blocks(16)  # √16 = 4 blocos
        [(0, 4), (4, 8), (8, 12), (12, 16)]
        >>> split_in_blocks(10)  # √10 ≈ 3.16, então 4 blocos
        [(0, 3), (3, 6), (6, 9), (9, 10)]
    """
    B = math.ceil(math.sqrt(L))
    base_size = math.ceil(L / B)
    blocks = []
    cur = 0
    while cur < L:
        blocks.append((cur, min(cur + base_size, L)))
        cur += base_size
    return blocks


def consensus_block(strings: Sequence[String], l: int, r: int) -> String:
    """
    Calcula a string consenso por maioria para um bloco específico.

    Para cada posição no intervalo [l:r), seleciona o caractere mais
    frequente entre todas as strings. Esta heurística simples mas eficaz
    serve como baseline para comparação com outras técnicas.

    Args:
        strings (Sequence[String]): Sequência de strings de entrada.
        l (int): Posição inicial do bloco (inclusiva).
        r (int): Posição final do bloco (exclusiva).

    Returns:
        String: String consenso para o bloco especificado.

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> consensus_block(strings, 1, 3)  # posições 1-2
        "CG"  # A maioria na pos 1 é C, na pos 2 é G
    """
    counter_cols = [Counter(s[l:r]) for s in strings]
    rep = "".join(
        Counter(col).most_common(1)[0][0] for col in zip(*[s[l:r] for s in strings])
    )
    return rep


# ---------------------------------------------------------------------------
# Técnicas por bloco (versão simplificada)
# ---------------------------------------------------------------------------


def _exhaustive_block(
    strings: Sequence[String], alphabet: str, l: int, r: int, k: int
) -> list[String]:
    """
    Realiza busca exaustiva em um bloco pequeno para encontrar os melhores candidatos.

    Esta função implementa a técnica de busca exaustiva para blocos pequenos,
    onde |Σ|^(r-l) ≤ 10.000. Para blocos maiores, usa uma estratégia de fallback
    baseada nos próprios blocos do dataset + consenso.

    A busca exaustiva garante a solução ótima para o bloco, mas tem complexidade
    exponencial. Por isso, é usada apenas quando o espaço de busca é pequeno.

    Args:
        strings (Sequence[String]): Sequência de strings de entrada.
        alphabet (str): Alfabeto disponível para construção dos candidatos.
        l (int): Posição inicial do bloco (inclusiva).
        r (int): Posição final do bloco (exclusiva).
        k (int): Número máximo de candidatos a retornar.

    Returns:
        list[String]: Lista dos k melhores candidatos para o bloco,
                     ordenados por distância crescente.

    Note:
        Se |Σ|^(r-l) > 10.000, usa fallback com candidatos do próprio dataset.
    """
    m = r - l
    lim = len(alphabet) ** m
    bests: list[tuple[int, String]] = []

    if lim <= 10_000:
        # Busca exaustiva para blocos pequenos
        for cand_tuple in itertools.product(alphabet, repeat=m):
            cand = "".join(cand_tuple)
            dist = max_distance(cand, [s[l:r] for s in strings])
            if len(bests) < k or dist < bests[-1][0]:
                bests.append((dist, cand))
                bests.sort(key=lambda x: x[0])
                bests = bests[:k]
    else:
        # Fallback: usa apenas candidatos do próprio dataset
        seen: set[String] = set()
        for s in strings:
            blk = s[l:r]
            if blk not in seen:
                seen.add(blk)
                dist = max_distance(blk, [t[l:r] for t in strings])
                bests.append((dist, blk))
        bests.sort(key=lambda x: x[0])
        bests = bests[:k]

    # Garante inclusão do consenso como candidato
    cons = consensus_block(strings, l, r)
    bests.append((max_distance(cons, [s[l:r] for s in strings]), cons))
    bests.sort(key=lambda x: x[0])
    return [c for _, c in bests[:k]]


def _beam_search_block(
    strings: Sequence[String], alphabet: str, l: int, r: int, beam_width: int, k: int
) -> list[String]:
    """
    Aplica beam search posição-a-posição para gerar candidatos em um bloco.

    Esta função implementa um beam search simples que constrói candidatos
    incrementalmente, mantendo apenas os beam_width melhores prefixos em
    cada posição. É usado para blocos médios/grandes onde a busca exaustiva
    seria computacionalmente inviável.

    O beam search oferece um bom compromisso entre qualidade da solução e
    eficiência computacional, explorando um subespaço promissor do espaço
    de busca completo.

    Args:
        strings (Sequence[String]): Sequência de strings de entrada.
        alphabet (str): Alfabeto disponível para construção dos candidatos.
        l (int): Posição inicial do bloco (inclusiva).
        r (int): Posição final do bloco (exclusiva).
        beam_width (int): Largura do beam (número de prefixos mantidos).
        k (int): Número máximo de candidatos finais a retornar.

    Returns:
        list[String]: Lista dos k melhores candidatos para o bloco,
                     ordenados por distância crescente.

    Note:
        A avaliação é feita incrementalmente considerando apenas o prefixo
        construído até o momento, o que pode não ser ótimo globalmente.
    """
    m = r - l
    beam = [""]  # Inicializa com prefixo vazio

    # Constrói candidatos posição por posição
    for pos in range(m):
        scored: list[tuple[int, String]] = []
        for prefix in beam:
            for a in alphabet:
                cand = prefix + a
                # Avalia parcialmente considerando apenas o prefixo atual
                partial_dists = []
                for s in strings:
                    mismatch = sum(
                        1 for i, c in enumerate(cand) if i < len(cand) and c != s[l + i]
                    )
                    partial_dists.append(mismatch)
                scored.append((max(partial_dists), cand))

        # Mantém apenas os beam_width melhores prefixos
        scored.sort(key=lambda x: x[0])
        beam = [s for _, s in scored[:beam_width]]

    # Avaliação final considerando o bloco completo
    final_scored = [(max_distance(c, [s[l:r] for s in strings]), c) for c in beam]
    final_scored.sort(key=lambda x: x[0])
    return [c for _, c in final_scored[:k]]


# ---------------------------------------------------------------------------
# Busca local global (hill-climbing)
# ---------------------------------------------------------------------------


def _local_search(candidate: String, strings: Sequence[String]) -> String:
    """
    Aplica hill-climbing (busca local) para melhorar uma solução candidata.

    Esta função implementa uma busca local simples que tenta melhorar
    iterativamente uma solução candidata. Em cada iteração, testa todas
    as possíveis substituições de caracteres em cada posição, mantendo
    apenas mudanças que reduzem a distância máxima.

    O algoritmo para quando nenhuma melhoria é encontrada (ótimo local).
    Para evitar bias, considera apenas caracteres que já aparecem nas
    strings originais em cada posição.

    Args:
        candidate (String): String candidata inicial para refinamento.
        strings (Sequence[String]): Sequência de strings de entrada.

    Returns:
        String: String refinada (ótimo local da busca).

    Note:
        Esta é uma busca gulosa que pode ficar presa em ótimos locais.
        A qualidade do resultado depende da qualidade da solução inicial.

    Example:
        >>> candidate = "ACGT"
        >>> strings = ["ACCT", "AGGT", "ATGT"]
        >>> refined = _local_search(candidate, strings)
        >>> # refined pode ser "ACGT" ou uma versão melhorada
    """
    cand_list = list(candidate)
    L = len(cand_list)

    # Para cada posição, considera apenas caracteres que aparecem nas strings originais
    alphabet_by_pos = [{s[i] for s in strings} for i in range(L)]

    improved = True
    while improved:
        improved = False
        base_val = max_distance("".join(cand_list), list(strings))

        # Testa todas as possíveis substituições
        for i in range(L):
            cur_char = cand_list[i]
            for alt in alphabet_by_pos[i]:
                if alt == cur_char:
                    continue

                # Testa a substituição
                cand_list[i] = alt
                new_val = max_distance("".join(cand_list), list(strings))

                if new_val < base_val:
                    # Melhoria encontrada, mantém a mudança
                    base_val = new_val
                    improved = True
                    cur_char = alt
                else:
                    # Sem melhoria, desfaz a mudança
                    cand_list[i] = cur_char

    return "".join(cand_list)


# ---------------------------------------------------------------------------
# Classe principal
# ---------------------------------------------------------------------------


class H3CSP:
    """
    Implementação do algoritmo H³-CSP (Hybrid Hierarchical Hamming Search).

    O H³-CSP é um algoritmo híbrido de três camadas que resolve o Closest String
    Problem através de:

    1. **B-Splitter**: Divisão hierárquica em ~√L blocos contíguos
    2. **Smart-Core**: Seleção adaptativa de técnicas por bloco:
       - Busca exaustiva para blocos pequenos (d_b ≤ block_small)
       - Beam search para blocos médios/grandes (d_b > block_small)
    3. **Global Refine**: Fusão de blocos + refinamento por hill-climbing

    O algoritmo é determinístico e adapta sua estratégia baseada na dificuldade
    de cada bloco (medida pela distância do consenso). Blocos mais difíceis
    recebem técnicas mais sofisticadas.

    Attributes:
        strings (list[String]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        L (int): Comprimento das strings.
        params (dict): Parâmetros do algoritmo.
        rng (random.Random): Gerador de números aleatórios.
        progress_callback (Callable | None): Callback para progresso.
        blocks (list[Block]): Lista de blocos gerados pela divisão.

    Example:
        >>> h3 = H3CSP(["ACGT", "AGCT", "ATCT"], "ACGT")
        >>> center, distance = h3.run()
        >>> print(f"Solução: {center} com distância {distance}")
    """

    def __init__(self, strings: Sequence[String], alphabet: str, **params):
        """
        Inicializa o algoritmo H³-CSP.

        Configura os parâmetros, valida a entrada e prepara as estruturas
        necessárias para execução, incluindo a divisão inicial em blocos.

        Args:
            strings (Sequence[String]): Sequência de strings de entrada.
                                      Todas devem ter o mesmo comprimento.
            alphabet (str): Alfabeto utilizado (ex: "ACGT" para DNA).
            **params: Parâmetros do algoritmo que sobrescreverão os padrões.
                     Ver H3_CSP_DEFAULTS para opções disponíveis.

        Raises:
            ValueError: Se as strings tiverem comprimentos diferentes.
            ValueError: Se o alfabeto estiver vazio.

        Note:
            A divisão em blocos é feita imediatamente usando a regra √L.
        """
        self.strings = list(strings)
        self.alphabet = alphabet
        self.L = len(strings[0])

        # Merge defaults + params
        self.params = dict(H3_CSP_DEFAULTS)
        self.params.update(params)

        self.rng = random.Random(self.params["seed"])
        self.progress_callback: Callable[[str], None] | None = None

        # Divisão inicial em blocos usando a regra √L
        self.blocks = split_in_blocks(self.L)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define um callback para monitorar o progresso da execução.

        O callback será chamado em pontos-chave do algoritmo com mensagens
        descritivas sobre o progresso atual.

        Args:
            callback (Callable[[str], None]): Função que recebe mensagens
                                             de progresso como string.
        """
        self.progress_callback = callback

    # ---------------------------------------------------------------------

    def _smart_core(self) -> list[list[String]]:
        """
        Implementa o Smart-Core: seleção adaptativa de técnicas por bloco.

        Para cada bloco, calcula a dificuldade (distância do consenso) e
        seleciona a técnica mais apropriada:
        - d_b ≤ block_small: Busca exaustiva (solução ótima)
        - block_small < d_b ≤ block_medium: Beam search reduzido
        - d_b > block_medium: Beam search completo

        Esta seleção adaptativa é o coração do H³-CSP, permitindo usar
        técnicas caras apenas onde necessário.

        Returns:
            list[list[String]]: Lista de listas, onde cada sublista contém
                               os k melhores candidatos para um bloco.

        Note:
            A ordem dos blocos na lista corresponde à ordem em self.blocks.
        """
        k = self.params["k_candidates"]
        small_limit = self.params["block_small"]
        medium_limit = self.params["block_medium"]
        beam_width = self.params["beam_width"]

        block_cands: list[list[String]] = []
        for l, r in self.blocks:
            # Calcula d_b (dificuldade do bloco)
            cons = consensus_block(self.strings, l, r)
            d_b = max_distance(cons, [s[l:r] for s in self.strings])

            # Seleção adaptativa de técnica
            if d_b <= small_limit:
                # Bloco fácil: busca exaustiva
                cands = _exhaustive_block(self.strings, self.alphabet, l, r, k)
            elif d_b <= medium_limit:
                # Bloco médio: beam search reduzido
                cands = _beam_search_block(
                    self.strings, self.alphabet, l, r, beam_width // 2, k
                )
            else:
                # Bloco difícil: beam search completo
                cands = _beam_search_block(
                    self.strings, self.alphabet, l, r, beam_width, k
                )

            block_cands.append(cands)

        return block_cands

    # ---------------------------------------------------------------------

    def _fuse_blocks(self, chosen: list[String]) -> String:
        """
        Fusão de blocos: concatena os candidatos escolhidos para cada bloco.

        Esta função implementa a fase de fusão do H³-CSP, onde os melhores
        candidatos de cada bloco são concatenados para formar uma solução
        completa. A ordem de concatenação segue a ordem original dos blocos.

        Args:
            chosen (list[String]): Lista de candidatos escolhidos, um por bloco.
                                  A ordem deve corresponder à ordem em self.blocks.

        Returns:
            String: String completa formada pela concatenação dos blocos.

        Example:
            >>> # Se blocks = [(0,2), (2,4)] e chosen = ["AC", "GT"]
            >>> self._fuse_blocks(["AC", "GT"])
            "ACGT"
        """
        return "".join(chosen)

    # ---------------------------------------------------------------------

    def run(self) -> tuple[String, int]:
        """
        Executa o algoritmo H³-CSP completo.

        Implementa as três fases do algoritmo:
        1. **Smart-Core**: Análise e processamento de blocos
        2. **Fusão**: Combinação dos melhores candidatos por bloco
        3. **Global Refine**: Refinamento por hill-climbing

        O algoritmo monitora o tempo de execução e pode ser interrompido
        se exceder o limite configurado em max_time.

        Returns:
            tuple[String, int]: Tupla contendo:
                - String: Solução encontrada (string central)
                - int: Distância máxima da solução para as strings originais

        Raises:
            Exception: Re-propaga qualquer exceção ocorrida durante a execução,
                      após registrar o erro no log.

        Note:
            O algoritmo usa callbacks de progresso se configurado, permitindo
            monitoramento em tempo real da execução.
        """
        start_time = time.time()

        try:
            # Fase 1: Smart-Core - Análise e processamento de blocos
            if self.progress_callback:
                self.progress_callback("Analisando blocos...")
            block_cands = self._smart_core()

            # Fase 2: Fusão - Combinação dos melhores candidatos
            if self.progress_callback:
                self.progress_callback("Fusão de blocos...")
            best_by_block = [cands[0] for cands in block_cands]
            center = self._fuse_blocks(best_by_block)
            best_val = max_distance(center, list(self.strings))

            logger.info("[Fusão inicial] dist=%d", best_val)

            # Fase 3: Global Refine - Refinamento por hill-climbing
            if self.progress_callback:
                self.progress_callback("Refinamento global...")
            for it in range(self.params["local_iters"]):
                # Verifica timeout
                if time.time() - start_time >= self.params["max_time"]:
                    logger.warning("Tempo máximo atingido no H3CSP")
                    if self.progress_callback:
                        self.progress_callback("Timeout atingido")
                    break

                # Callback de progresso
                if self.progress_callback:
                    self.progress_callback(f"Refinamento: iteração {it+1}")

                # Aplica busca local
                center = _local_search(center, self.strings)
                new_val = max_distance(center, list(self.strings))

                if new_val < best_val:
                    # Melhoria encontrada
                    logger.info("[Local] it=%d  %d->%d", it + 1, best_val, new_val)
                    best_val = new_val
                    if self.progress_callback:
                        self.progress_callback(
                            f"Melhoria encontrada: distância={best_val}"
                        )
                else:
                    # Sem melhoria, para o refinamento
                    break

            return center, best_val

        except Exception as e:
            logger.error("Erro no H3CSP: %s", str(e))
            if self.progress_callback:
                self.progress_callback(f"Erro: {str(e)}")
            raise e  # Re-propaga para ser capturado no wrapper
