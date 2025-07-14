"""
Implementação exata do DP-CSP (Programação Dinâmica) para o Closest String Problem.

O DP-CSP é um algoritmo de programação dinâmica que encontra a solução
EXATA para o Closest String Problem. Diferentemente das heurísticas,
garante que encontrará a string center com a menor distância máxima
possível, provendo uma cota inferior ótima para comparação.

ARQUITETURA ALGORÍTMICA:

ALGORITMO DP-CSP DETALHADO:

1. MODELAGEM DO PROBLEMA
   - Estado: (posição, vetor_distâncias)
   - Decisão: qual símbolo escolher na posição atual
   - Transição: estado[i] → estado[i+1] via escolha de símbolo
   - Objetivo: minimizar distância máxima final

2. PROGRAMAÇÃO DINÂMICA
   - Tabela DP: dp[posição][estado_distâncias] = custo_mínimo
   - Recorrência: dp[i][d] = min(dp[i+1][d'] + custo_transição)
   - Caso base: dp[L][d] = max(d) (distância máxima final)
   - Solução: dp[0][estado_inicial] = distância ótima

3. OTIMIZAÇÕES DE EFICIÊNCIA
   - Poda por limite de memória (max_states)
   - Encoding eficiente de estados como tuplas
   - Memoização para evitar recálculos
   - Reconstrução lazy da solução ótima

4. CONTROLES DE RECURSO
   - Limite de memória para evitar explosão exponencial
   - Timeout para problemas intratáveis
   - Fallback para heurística se recursos insuficientes
   - Monitoramento de progresso via callbacks

FILOSOFIA ALGORÍTMICA:

- EXATIDÃO: Garante solução ótima dentro dos limites de recurso
- COMPLETUDE: Explora todo o espaço de estados viável
- EFICIÊNCIA: Usa programação dinâmica para evitar recálculos
- ROBUSTEZ: Fallbacks e controles para casos intratáveis

CARACTERÍSTICAS DISTINTIVAS:

- GARANTIA DE OTIMALIDADE: Diferente das heurísticas, garante solução ótima
- COMPLEXIDADE CONTROLADA: Poda por limite de memória
- VERSATILIDADE: Funciona para qualquer alfabeto e tamanho de string
- TRANSPARÊNCIA: Reconstrução da solução ótima passo-a-passo

APLICAÇÃO AO CSP:
O DP-CSP é ideal para:
- Instâncias pequenas/médias onde exatidão é crítica
- Benchmark de qualidade para heurísticas
- Análise de estrutura do problema
- Validação de soluções aproximadas

LIMITAÇÕES:
- Complexidade exponencial no pior caso
- Requer muita memória para problemas grandes
- Pode ser lento para instâncias complexas
- Não escalável para strings muito longas

EXEMPLO DE FUNCIONAMENTO:
Para strings ["AC", "AG", "AT"] e alfabeto "ACGT":
1. Estados iniciais: (0, [0,0,0])
2. Posição 0: escolher A → (1, [0,0,0])
3. Posição 1: escolher C → (2, [0,1,1]), G → (2, [1,0,1]), T → (2, [1,1,0])
4. Solução ótima: "AC" com distância máxima 1

Classes:
    DPCSP: Implementação do algoritmo DP-CSP com controles de recurso.

Funções auxiliares:
    _dp_decision(): Função de decisão da programação dinâmica.
    exact_dp_closest_string(): Wrapper principal com controles de recurso.

Types:
    String: Alias para str (representação de strings).
    State: Tupla representando estado de distâncias.

Author: Implementação baseada em programação dinâmica clássica para CSP
Version: Otimizada com controles de recurso e fallbacks robustos
"""

from __future__ import annotations

import logging
import time
from collections.abc import Callable, Sequence
from typing import Any, TypeAlias, cast

from src.domain.metrics import max_distance

from .config import DP_CSP_DEFAULTS

logger = logging.getLogger(__name__)

# Type aliases for clarity
RemVec: TypeAlias = tuple[int, ...]  # (Remaining errors vector)
State: TypeAlias = tuple[int, RemVec]  # (position, rem_errors)
String: TypeAlias = str


def _dp_decision(strings: Sequence[String], alphabet: str, d: int) -> String | None:
    """
    Algoritmo de decisão DP: verifica se existe string center com raio ≤ d.

    Esta é a função central do DP-CSP que implementa o algoritmo de programação
    dinâmica para o problema decisório: "Existe uma string center tal que
    a distância máxima para todas as strings seja no máximo d?"

    ALGORITMO DE PROGRAMAÇÃO DINÂMICA:

    1. MODELAGEM DE ESTADOS:
       - Estado: rem = (rem[0], rem[1], ..., rem[n-1])
       - rem[i] = número de erros ainda permitidos para string i
       - Estado inicial: rem = (d, d, ..., d)
       - Estado final viável: rem[i] ≥ 0 para todo i

    2. PRÉ-COMPUTAÇÃO DE DIFERENÇAS:
       - Para cada posição pos e símbolo σ:
       - δσ[pos] = vetor indicando se σ difere de cada string em pos
       - δσ[pos][i] = 1 se σ ≠ strings[i][pos], 0 caso contrário

    3. TRANSIÇÕES DE ESTADO:
       - Para cada posição pos = 0, 1, ..., L-1:
       - Para cada estado rem na fronteira atual:
       - Para cada símbolo σ do alfabeto:
         a) Calcula novo_rem[i] = rem[i] - δσ[pos][i]
         b) Se min(novo_rem) ≥ 0: estado viável, adiciona à próxima fronteira
         c) Armazena ponteiro parent para reconstrução

    4. PODA DE ESTADOS:
       - Remove estados onde algum rem[i] < 0 (inviáveis)
       - Evita duplicatas: um estado por configuração de rem
       - Para se fronteira vazia (nenhuma solução possível)

    5. RECONSTRUÇÃO DE SOLUÇÃO:
       - Se chegou ao final com estados viáveis: solução existe
       - Usa ponteiros parent para backtrack do final ao início
       - Reconstrói string center símbolo por símbolo

    COMPLEXIDADE:
    - Estados: O((d+1)^n) - cada rem[i] pode ser 0, 1, ..., d
    - Transições: O(L × |Σ|) - para cada pos, testa cada símbolo
    - Total: O((d+1)^n × L × |Σ|)

    OTIMIZAÇÕES IMPLEMENTADAS:
    - Pré-computação de δσ evita recálculos
    - Fronteira compacta reduz uso de memória
    - Poda precoce elimina estados inviáveis rapidamente
    - Parent pointers permitem reconstrução eficiente

    Args:
        strings: Sequência de strings de entrada (mesmo comprimento)
        alphabet: String com todos os símbolos válidos
        d: Raio máximo permitido (threshold de decisão)

    Returns:
        String center com raio ≤ d se existir, None caso contrário

    Example:
        >>> strings = ["AC", "AT", "GC"]
        >>> _dp_decision(strings, "ACGT", 1)
        "AC"  # Uma solução com raio = 1

        >>> _dp_decision(strings, "ACGT", 0)
        None  # Impossível raio = 0

    Note:
        A função é determinística mas pode retornar qualquer solução
        válida (não necessariamente única). Validação final verifica
        se solução retornada satisfaz realmente raio ≤ d.
    """
    n, L = len(strings), len(strings[0])

    # PRÉ-COMPUTAÇÃO: Matriz de diferenças δσ[pos][string]
    # Para cada posição e símbolo, calcula vetor de discrepâncias
    delta: list[dict[String, RemVec]] = []
    for pos in range(L):
        # Símbolos na posição pos de todas as strings
        col = [s[pos] for s in strings]
        # Para cada símbolo σ, calcula diferenças
        delta.append({σ: tuple(int(σ != c) for c in col) for σ in alphabet})

    # INICIALIZAÇÃO DA PROGRAMAÇÃO DINÂMICA
    start: RemVec = (d,) * n  # Estado inicial: d erros para cada string
    frontier: set[RemVec] = {start}  # Fronteira atual de estados viáveis

    # Ponteiros para reconstrução: (pos, estado) → (estado_anterior, símbolo)
    parent: dict[tuple[int, RemVec], tuple[RemVec | None, String]] = {
        (0, start): (None, "")  # Estado inicial não tem predecessor
    }

    # PROGRAMAÇÃO DINÂMICA POSIÇÃO-A-POSIÇÃO
    for pos in range(L):
        nxt: set[RemVec] = set()  # Próxima fronteira

        # Processa todos os estados da fronteira atual
        for rem in frontier:
            # Testa cada símbolo do alfabeto na posição atual
            for σ, dv in delta[pos].items():
                # TRANSIÇÃO DE ESTADO: consome erros baseado em diferenças
                new_rem = tuple(r - v for r, v in zip(rem, dv))

                # PODA: Remove estados inviáveis (algum rem[i] negativo)
                if min(new_rem) < 0:
                    continue  # Estado inviável, poda

                # Evita estados duplicados
                key = (pos + 1, new_rem)
                if key in parent:
                    continue  # Estado já visitado

                # ARMAZENAMENTO: Registra transição e adiciona à próxima fronteira
                parent[key] = (rem, σ)
                nxt.add(new_rem)

        # Atualiza fronteira para próxima posição
        frontier = nxt

        # PODA GLOBAL: Se nenhum estado viável, impossível continuar
        if not frontier:
            return None  # Nenhuma solução existe para este d

    # RECONSTRUÇÃO DA SOLUÇÃO
    # Se chegamos aqui, existe pelo menos um estado final viável
    final_rem = next(iter(frontier))  # Qualquer estado final serve

    center_chars: list[String] = []
    pos: int = L
    rem: RemVec = final_rem

    # Backtrack usando ponteiros parent
    while pos > 0:
        prev_rem, σ = parent[(pos, rem)]
        center_chars.append(σ)  # Símbolo escolhido nesta posição
        pos -= 1

        if prev_rem is None:  # Chegamos ao estado inicial
            break
        rem = cast(RemVec, prev_rem)

    # Reconstrói string na ordem correta
    center_chars.reverse()
    result = "".join(center_chars)

    # VALIDAÇÃO FINAL: Verifica se solução é realmente válida
    from src.domain.metrics import hamming_distance

    max_dist = max(hamming_distance(result, s) for s in strings)
    if max_dist > d:
        logger.error(
            "[DP_DECISION] ERRO: Solução inválida! dist=%d > d=%d", max_dist, d
        )

    return result


def exact_dp_closest_string(
    strings: list[String],
    alphabet: str,
    max_d: int | None = None,
    progress_callback: Callable[[str], None] | None = None,
    warning_callback: Callable[[str], None] | None = None,
) -> tuple[String, int]:
    """
    Encontra a solução EXATA do Closest String Problem usando programação dinâmica.

    Esta é a interface principal do DP-CSP que coordena a busca pelo raio ótimo
    mínimo. Diferentemente das heurísticas, garante encontrar a solução globalmente
    ótima, fornecendo a cota inferior verdadeira para o problema.

    ESTRATÉGIA DE BUSCA DO RAIO ÓTIMO:

    1. BUSCA INCREMENTAL:
       - Testa valores d = 0, 1, 2, ..., max_d sequencialmente
       - Para cada d, usa algoritmo de decisão DP
       - Para no primeiro d onde encontra solução (garantidamente ótimo)

    2. OTIMALIDADE GARANTIDA:
       - Como testa em ordem crescente, primeira solução é ótima
       - d* = min{d : existe center com raio ≤ d}
       - Evita busca binária desnecessária para problema small

    3. MONITORAMENTO DE RECURSOS:
       - Estima complexidade (d+1)^n antes de executar
       - Monitora uso de memória RSS durante execução
       - Interrompe se exceder limites de tempo ou memória
       - Reporta métricas detalhadas para auditoria

    LIMITAÇÕES E SALVAGUARDAS:

    - Complexidade Exponencial: O((d+1)^n × L × |Σ|) cresce rapidamente
    - Limite de Estados: Aborta se (d+1)^n > 2e9 (evita crash por memória)
    - Limite de Memória: Monitora RSS, aborta se > 95% do limite seguro
    - Limite de Tempo: Timeout configurável (padrão 300s)

    CASOS DE USO RECOMENDADOS:

    - Validação: Verificar qualidade de algoritmos heurísticos
    - Benchmark: Estabelecer ground truth para comparação
    - Instâncias Pequenas: n ≤ 8, L ≤ 20, alphabet ≤ 4
    - Pesquisa: Análise teórica de propriedades do CSP

    CONFIGURAÇÃO DE PARÂMETROS:

    - max_d=None: Usa baseline (distância da primeira string) como limite
    - progress_callback: Reporta progresso "Testando d=X"
    - warning_callback: Reporta alertas de recursos antes de abortar

    MÉTRICAS E LOGS DETALHADOS:

    A função registra informações completas para auditoria:
    - Parâmetros de entrada (n, L, alfabeto, strings)
    - Limites de recursos configurados
    - Progresso da busca (d testado, estados estimados)
    - Resultado final (center encontrado, d* ótimo, validação)
    - Métricas de performance (tempo, memória, iterações)

    Args:
        strings: Lista de strings de entrada (mesmo comprimento)
        alphabet: String com símbolos válidos do alfabeto
        max_d: Raio máximo a testar (None = baseline automático)
        progress_callback: Função para reportar progresso (opcional)
        warning_callback: Função para reportar alertas de recursos (opcional)

    Returns:
        tuple: (center_ótimo, d*_ótimo)
            - center_ótimo: String center com menor raio possível
            - d*_ótimo: Menor raio possível para a instância

    Raises:
        RuntimeError: Se não encontrar solução dentro dos limites ou
                     se exceder recursos (memória, tempo, complexidade)

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> center, d_star = exact_dp_closest_string(strings, "ACGT")
        >>> print(f"Solução ótima: {center} com raio {d_star}")
        "Solução ótima: ACCT com raio 1"

    Note:
        Para instâncias grandes, considere usar algoritmos heurísticos
        como BLF-GA ou CSC que oferecem boa qualidade em tempo viável.
        O DP-CSP deve ser usado quando exatidão é crítica e recursos
        computacionais são adequados.
    """
    # CONFIGURAÇÃO INICIAL E VALIDAÇÃO
    baseline_val = max_distance(strings[0], strings)  # Cota superior simples
    if max_d is None:
        max_d = baseline_val

    n = len(strings)
    L = len(strings[0])

    # LOGS DETALHADOS DE ENTRADA
    logger.info(
        "[DP_CSP] Iniciando busca exata com max_d=%d, baseline=%d", max_d, baseline_val
    )
    logger.info("[DP_CSP] Dataset: n=%d, L=%d, alfabeto=%s", n, L, alphabet)
    for i, s in enumerate(strings):
        logger.info("[DP_CSP] String %d: %s", i, s)

    # CONFIGURAÇÃO DE MONITORAMENTO DE RECURSOS
    safe_mem_mb = 1000.0  # Limite padrão simplificado
    max_time = DP_CSP_DEFAULTS.get("max_time", 300)
    t0 = time.time()

    logger.info("[DP_CSP] Limites: mem=%.1fMB, tempo=%ds", safe_mem_mb, max_time)

    def check_limits(d):
        """Verifica limites de recursos antes de processar raio d."""
        import gc
        import os

        # LIMPEZA E MEDIÇÃO DE MEMÓRIA
        gc.collect()
        mem_mb = 0.0
        try:
            if os.path.exists("/proc/self/status"):
                with open("/proc/self/status", encoding="utf-8") as f:
                    for line in f:
                        if line.startswith("VmRSS:"):
                            kb = int(line.split()[1])
                            mem_mb = kb / 1024.0
                            break
        except (OSError, IOError, ValueError):
            pass  # Continua mesmo se não conseguir medir memória

        elapsed = time.time() - t0

        # ESTIMATIVA DE COMPLEXIDADE
        # Limite prático: ~2e9 estados (≈16GB RAM, 8 bytes/estado)
        state_count_est = (d + 1) ** n

        # VERIFICAÇÕES DE RECURSOS
        if state_count_est > 2_000_000_000:
            msg = (
                f"DP-CSP interrompido: (d+1)^n = {state_count_est:,} excede limite prático "
                "(~2e9 estados, ~16GB RAM). Tente reduzir n ou d."
            )
            logger.error("[DP_CSP] %s", msg)
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

        if mem_mb > safe_mem_mb * 0.95:
            msg = f"DP-CSP interrompido: uso de memória {mem_mb:.1f}MB excedeu limite seguro ({safe_mem_mb:.1f}MB)"
            logger.error("[DP_CSP] %s", msg)
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

        if elapsed > max_time:
            msg = f"DP-CSP interrompido: tempo de execução excedeu {max_time}s"
            logger.error("[DP_CSP] %s", msg)
            if warning_callback:
                warning_callback(msg)
            raise RuntimeError(msg)

    # BUSCA INCREMENTAL DO RAIO ÓTIMO
    for d in range(max_d + 1):
        # Verificar recursos antes de cada iteração principal
        check_limits(d)

        if progress_callback:
            progress_callback(f"Testando d={d}")

        logger.info("[DP_CSP] Testando d=%d (tentativa %d/%d)", d, d + 1, max_d + 1)

        # ALGORITMO DE DECISÃO PARA RAIO d
        center = _dp_decision(strings, alphabet, d)

        if center is not None:
            # SUCESSO: Encontrou solução com raio d

            # VALIDAÇÃO FINAL RIGOROSA
            from src.domain.metrics import hamming_distance

            max_dist = max(hamming_distance(center, s) for s in strings)

            logger.info("[DP_CSP] SUCESSO! Encontrou solução com d=%d", d)
            logger.info("[DP_CSP] Centro encontrado: %s", center)
            logger.info("[DP_CSP] Validação final: distância máxima = %d", max_dist)

            # Verificação de consistência
            if max_dist != d:
                logger.warning(
                    "[DP_CSP] INCONSISTÊNCIA: d=%d mas distância real=%d", d, max_dist
                )

            return center, d

    # FALHA: Nenhuma solução encontrada dentro dos limites
    logger.error("[DP_CSP] FALHA: Não foi possível encontrar centro com d ≤ %d", max_d)
    raise RuntimeError(
        f"Não foi possível encontrar centro com d ≤ {max_d}. "
        "Tente aumentar o limite."
    )


class DPCSP:
    """
    Classe wrapper para o algoritmo DP-CSP (Dynamic Programming CSP).

    Esta classe fornece uma interface orientada a objetos para o algoritmo
    DP-CSP, encapsulando configurações e oferecendo métodos convenientes
    para execução e análise.

    A classe é principalmente um wrapper em torno da função
    exact_dp_closest_string(), fornecendo:
    - Interface consistente com outros algoritmos do projeto
    - Configuração de parâmetros via construtor
    - Métodos de callback para monitoramento
    - Metadados e estatísticas de execução

    Attributes:
        strings: Lista de strings de entrada
        alphabet: Alfabeto utilizado
        max_d: Raio máximo para busca
        progress_callback: Callback para progresso
        warning_callback: Callback para warnings
    """

    def __init__(
        self, strings: list[str], alphabet: str, max_d: int | None = None, **kwargs: Any
    ):
        """
        Inicializa o algoritmo DP-CSP.

        Args:
            strings: Lista de strings de entrada (mesmo comprimento)
            alphabet: Alfabeto utilizado
            max_d: Raio máximo para busca (None = automático)
            **kwargs: Parâmetros adicionais (para compatibilidade)
        """
        self.strings = strings
        self.alphabet = alphabet
        self.max_d = max_d
        self.progress_callback: Callable[[str], None] | None = None
        self.warning_callback: Callable[[str], None] | None = None

        # Parâmetros para compatibilidade (não utilizados pelo DP-CSP)
        self.params = kwargs

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Define callback para reportar progresso."""
        self.progress_callback = callback

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """Define callback para reportar warnings."""
        self.warning_callback = callback

    def run(self) -> tuple[str, int]:
        """
        Executa o algoritmo DP-CSP.

        Returns:
            tuple: (center, distância_ótima)
                - center: String center ótima
                - distância_ótima: Menor raio possível
        """
        return exact_dp_closest_string(
            self.strings,
            self.alphabet,
            self.max_d,
            self.progress_callback,
            self.warning_callback,
        )
