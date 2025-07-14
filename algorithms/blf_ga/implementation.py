"""
Implementação da heurística BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) para CSP.

O BLF-GA é uma metaheurística híbrida avançada que combina três estratégias principais:
1. **Blockwise Learning**: Divide as strings em blocos e aprende padrões locais
2. **Genetic Algorithm**: Usa operações genéticas para evolução global
3. **Fusion**: Combina conhecimento local e global para melhor convergência

ARQUITETURA ALGORÍTMICA:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                           ALGORITMO BLF-GA DETALHADO                            │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. INICIALIZAÇÃO                                                               │
│   ├── Criar população inicial híbrida:                                         │
│   │   ├── Consenso das strings originais (qualidade inicial)                  │
│   │   ├── Variações do consenso por blocos (diversidade inteligente)          │
│   │   └── Strings aleatórias (exploração ampla)                               │
│   └── Dividir strings em blocos adaptativos baseado em initial_blocks         │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. LOOP PRINCIPAL (até critério de parada)                                     │
│   ├── A) APRENDIZADO POR BLOCOS (Blockwise Learning)                          │
│   │   ├── Para cada bloco, encontrar padrão consenso da população atual       │
│   │   ├── Criar repositório de conhecimento local                             │
│   │   └── Usar padrões para guiar operações genéticas                         │
│   ├── B) EVOLUÇÃO GENÉTICA                                                     │
│   │   ├── Seleção por torneio (pressão seletiva balanceada)                   │
│   │   ├── Crossover adaptativo:                                               │
│   │   │   ├── one-point: preserva segmentos grandes                           │
│   │   │   ├── uniform: máxima recombinação                                     │
│   │   │   └── blend-blocks: preserva padrões locais                           │
│   │   ├── Mutação adaptativa:                                                 │
│   │   │   ├── multi: alterações pontuais                                      │
│   │   │   ├── inversion: reorganização estrutural                             │
│   │   │   └── transposition: relocação de segmentos                           │
│   │   └── Elitismo adaptativo (pode ser desabilitado para diversidade)       │
│   ├── C) MECANISMOS ADAPTATIVOS                                               │
│   │   ├── Redivisão de blocos baseada em entropia posicional                  │
│   │   ├── Imigrantes aleatórios para manter diversidade                       │
│   │   ├── Mutação adaptativa baseada em convergência e diversidade            │
│   │   ├── Refinamento local intensivo dos melhores indivíduos                 │
│   │   ├── Niching para preservar múltiplas soluções                           │
│   │   └── Restart populacional em caso de estagnação                          │
│   └── D) CRITÉRIOS DE PARADA                                                  │
│       ├── Solução ótima encontrada (distância = 0)                            │
│       ├── Limite de tempo/gerações alcançado                                  │
│       └── Early stopping por estagnação (no_improve_patience)                 │
└─────────────────────────────────────────────────────────────────────────────────┘

CARACTERÍSTICAS DISTINTIVAS:

• **HIBRIDIZAÇÃO INTELIGENTE**: Combina busca global (AG) com aprendizado local (blocos)
• **ADAPTATIVIDADE DINÂMICA**: Parâmetros se ajustam automaticamente durante evolução
• **DIVERSIDADE CONTROLADA**: Múltiplos mecanismos para evitar convergência prematura
• **PARALELIZAÇÃO EFICIENTE**: Aproveita múltiplos cores para avaliação e refinamento
• **CONFIGURABILIDADE AVANÇADA**: Parâmetros dinâmicos e estratégias intercambiáveis

APLICAÇÃO AO CSP:
O BLF-GA é especialmente adequado para o Closest String Problem pois:
- Blocos capturam padrões locais comuns entre strings
- Aprendizado adapta-se à estrutura específica do problema
- Hibridização balanceia exploração global com intensificação local
- Mecanismos adaptativos evitam armadilhas de convergência prematura

EXEMPLO DE FLUXO:
Strings: ["ACGT", "AGCT", "ATCT"]
1. Consenso inicial: "ACCT"
2. Blocos: [(0,2), (2,4)] baseado em entropia
3. Aprendizado: bloco1="AC", bloco2="CT"
4. Evolução: crossover_blend_blocks usa padrões aprendidos
5. Adaptação: redivide blocos se entropia mudar
6. Convergência: refina soluções até ótimo ou critério de parada

Classes:
    BLFGA: Implementa o algoritmo BLF-GA com todos os mecanismos adaptativos.

Funções auxiliares:
    hamming_dist(a, b): Wrapper para distância de Hamming.
    max_distance(center, strings): Wrapper para distância máxima.

Author: Implementação baseada em pesquisa de metaheurísticas híbridas
Version: Otimizada para CSP com mecanismos adaptativos avançados
"""

import logging
import os
import random
import time
from collections import Counter
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

import numpy as np

from src.domain.metrics import hamming_distance, max_distance

from .config import BLF_GA_DEFAULTS
from .ops import genetic_ops

logger = logging.getLogger(__name__)

String = str
Population = list[String]


def hamming_dist(a: String, b: String) -> int:
    """
    Calcula a distância de Hamming entre duas strings.

    Wrapper para manter compatibilidade com a interface esperada pelo BLF-GA,
    delegando o cálculo para a implementação otimizada em src.utils.distance.

    A distância de Hamming é o número de posições onde os caracteres diferem
    entre duas strings de mesmo comprimento. É fundamental para o CSP pois
    define a função de fitness do problema.

    Args:
        a: Primeira string para comparação
        b: Segunda string para comparação (deve ter mesmo comprimento que a)

    Returns:
        int: Número de posições onde as strings diferem

    Example:
        >>> hamming_dist("ACGT", "AGTC")
        2  # Diferem nas posições 1 e 2
    """
    return hamming_distance(a, b)


class BLFGA:
    """
    Implementação do algoritmo BLF-GA (Blockwise Learning Fusion + Genetic Algorithm).

    O BLF-GA é uma metaheurística híbrida que resolve o Closest String Problem através
    da combinação de três componentes principais:

    1. **BLOCKWISE LEARNING (Aprendizado por Blocos)**:
       - Divide as strings de entrada em blocos menores
       - Aprende padrões locais em cada bloco através de consenso
       - Permite exploração mais eficiente do espaço de busca

    2. **GENETIC ALGORITHM (Algoritmo Genético)**:
       - Mantém uma população de soluções candidatas
       - Usa operações genéticas (seleção, crossover, mutação) para evolução
       - Implementa elitismo para preservar melhores soluções

    3. **FUSION (Fusão)**:
       - Combina conhecimento local (blocos) com busca global (GA)
       - Usa refinamento local intensivo nos melhores indivíduos
       - Aplica mecanismos adaptativos para balancear exploração/exploração

    CARACTERÍSTICAS PRINCIPAIS:
    - Blocking adaptativo baseado em entropia das posições
    - Múltiplos operadores de crossover e mutação
    - Controle de diversidade com imigrantes aleatórios
    - Mutação adaptativa baseada em convergência
    - Refinamento local dos elites
    - Múltiplos critérios de parada

    Args:
        strings: Lista de strings de entrada para encontrar o centro
        alphabet: Alfabeto das strings (ex: 'ATCG' para DNA)
        pop_size: Tamanho da população (int fixo ou float para multiplicador de n)
        min_pop_size: Tamanho mínimo da população quando usar proporção
        initial_blocks: Número inicial de blocos (int fixo ou float para proporção)
        min_block_len: Tamanho mínimo de cada bloco
        cross_prob: Probabilidade de crossover [0.6-0.95]
        mut_prob: Probabilidade de mutação [0.01-0.2]
        elite_rate: Taxa de elitismo [0.01-0.1]
        rediv_freq: Frequência de redivisão de blocos (gerações)
        max_gens: Número máximo de gerações
        max_time: Tempo máximo de execução (segundos)
        seed: Semente para reprodutibilidade

    Demais parâmetros controlam aspectos específicos como:
    - Diversidade e imigrantes
    - Mutação adaptativa
    - Tipos de operadores genéticos
    - Refinamento local
    - Critérios de parada e reinício
    """

    def __init__(
        self,
        strings: list[String],
        alphabet: str,
        pop_size: int | float | None = None,
        min_pop_size: int | None = None,  # pylint: disable=unused-argument
        initial_blocks: int | float | None = None,
        min_block_len: int | None = None,
        cross_prob: float | None = None,
        mut_prob: float | None = None,
        elite_rate: float | None = None,
        rediv_freq: int | None = None,
        max_gens: int | None = None,
        max_time: float | None = None,
        seed: int | None = None,
        # Novos parâmetros
        immigrant_freq: int | None = None,  # pylint: disable=unused-argument
        immigrant_ratio: float | None = None,  # pylint: disable=unused-argument
        diversity_threshold: float | None = None,  # pylint: disable=unused-argument
        mutation_adapt_N: int | None = None,  # pylint: disable=unused-argument
        mutation_adapt_factor: float | None = None,  # pylint: disable=unused-argument
        mutation_adapt_duration: int | None = None,  # pylint: disable=unused-argument
        mutation_type: str | None = None,  # pylint: disable=unused-argument
        mutation_multi_n: int | None = None,  # pylint: disable=unused-argument
        tournament_k: int | None = None,  # pylint: disable=unused-argument
        crossover_type: str | None = None,  # pylint: disable=unused-argument
        niching: bool | None = None,  # pylint: disable=unused-argument
        niching_radius: int | None = None,  # pylint: disable=unused-argument
        refinement_type: str | None = None,  # pylint: disable=unused-argument
        refine_elites: str | None = None,  # pylint: disable=unused-argument
        refine_iter_limit: int | None = None,  # pylint: disable=unused-argument
        restart_patience: int | None = None,  # pylint: disable=unused-argument
        restart_ratio: float | None = None,  # pylint: disable=unused-argument
        disable_elitism_gens: int | None = None,  # pylint: disable=unused-argument
        no_improve_patience: int | None = None,  # pylint: disable=unused-argument
    ):
        self.strings = strings
        self.n = len(strings)
        self.L = len(strings[0])
        self.alphabet = alphabet

        # Configurar workers internos a partir da variável de ambiente
        self.internal_workers = int(os.environ.get("INTERNAL_WORKERS", "1"))

        # Carrega parâmetros do dicionário de defaults, permite sobrescrever via argumentos
        params = {**BLF_GA_DEFAULTS}
        if pop_size is not None:
            params["pop_size"] = pop_size
        if initial_blocks is not None:
            params["initial_blocks"] = initial_blocks
        if min_block_len is not None:
            params["min_block_len"] = min_block_len
        if cross_prob is not None:
            params["cross_prob"] = cross_prob
        if mut_prob is not None:
            params["mut_prob"] = mut_prob
        if elite_rate is not None:
            params["elite_rate"] = elite_rate
        if rediv_freq is not None:
            params["rediv_freq"] = rediv_freq
        if max_gens is not None:
            params["max_gens"] = max_gens
        if max_time is not None:
            params["max_time"] = max_time
        if seed is not None:
            params["seed"] = seed

        # Incluir parâmetro min_pop_size se fornecido
        if min_pop_size is not None:
            params["min_pop_size"] = min_pop_size

        # Novos parâmetros
        self.immigrant_freq = params["immigrant_freq"]
        self.immigrant_ratio = params["immigrant_ratio"]
        self.diversity_threshold = params["diversity_threshold"]
        self.mutation_adapt_N = params["mutation_adapt_N"]
        self.mutation_adapt_factor = params["mutation_adapt_factor"]
        self.mutation_adapt_duration = params["mutation_adapt_duration"]
        self.mutation_type = params["mutation_type"]
        self.mutation_multi_n = params["mutation_multi_n"]
        self.tournament_k = params["tournament_k"]
        self.crossover_type = params["crossover_type"]
        self.min_pop_size = params["min_pop_size"]  # Tamanho mínimo da população
        self.niching = params["niching"]
        self.niching_radius = params["niching_radius"]
        self.refinement_type = params["refinement_type"]
        self.refine_elites = params["refine_elites"]
        self.refine_iter_limit = params["refine_iter_limit"]
        self.restart_patience = params["restart_patience"]
        self.restart_ratio = params["restart_ratio"]
        self.disable_elitism_gens = params["disable_elitism_gens"]

        # Conversão dinâmica de pop_size e initial_blocks se forem proporções
        self.pop_size = self._resolve_dynamic_param(
            params["pop_size"],
            self.n,
            mode="multiplier",
            min_value=params["min_pop_size"],
        )
        self.initial_blocks = self._resolve_dynamic_param(
            params["initial_blocks"], self.L, mode="proportion"
        )
        self.min_block_len = params["min_block_len"]
        self.cross_prob = params["cross_prob"]
        self.mut_prob = params["mut_prob"]
        self.elite_rate = params["elite_rate"]
        self.rediv_freq = params["rediv_freq"]
        self.max_gens = params["max_gens"]
        self.max_time = params["max_time"]
        self.rng = random.Random(params["seed"])
        self.progress_callback: Callable[[str], None] | None = None
        self.history_callback: Callable[[int, dict], None] | None = (
            None  # Callback para eventos dinâmicos
        )

        # Inicializa os blocos após todos os parâmetros necessários
        self.blocks = self._initial_blocking()
        self.history = []  # Histórico de distâncias por geração

        # Early stopping: converte no_improve_patience para int se for proporção
        patience_param = params.get("no_improve_patience", 0)
        if isinstance(patience_param, float) and 0 < patience_param < 1:
            self.no_improve_patience = max(1, int(patience_param * self.max_gens))
        else:
            self.no_improve_patience = int(patience_param)

    @staticmethod
    def _resolve_dynamic_param(param, ref_value, mode="proportion", min_value=None):
        """
        Resolve parâmetros que podem ser especificados como valores absolutos ou relativos.

        Esta função centraliza a lógica de conversão de parâmetros dinâmicos,
        permitindo flexibilidade na configuração do algoritmo. Por exemplo,
        pop_size pode ser especificado como um número fixo (50) ou como
        multiplicador do número de strings (2.0 = duas vezes o número de strings).

        Modos de operação:
        - 'multiplier': Para pop_size - float >= 1 é multiplicador, int é absoluto
        - 'proportion': Para initial_blocks - float 0-1 é proporção, int é absoluto

        Args:
            param: Parâmetro a ser resolvido (int ou float)
            ref_value: Valor de referência para cálculo (n para pop_size, L para blocks)
            mode: Modo de resolução ('multiplier' ou 'proportion')
            min_value: Valor mínimo a ser garantido (opcional)

        Returns:
            int: Valor resolvido do parâmetro

        Examples:
            >>> _resolve_dynamic_param(2.0, 10, 'multiplier')  # 20
            >>> _resolve_dynamic_param(50, 10, 'multiplier')   # 50
            >>> _resolve_dynamic_param(0.2, 100, 'proportion') # 20
            >>> _resolve_dynamic_param(5, 100, 'proportion')   # 5
        """
        if mode == "multiplier":
            # Para pop_size: float >= 1 é multiplicador, int é valor absoluto
            if isinstance(param, float) and param >= 1:
                result = int(param * ref_value)
            else:
                result = int(param)
        else:  # mode == "proportion"
            # Para initial_blocks: float 0-1 é proporção, int é valor absoluto
            if isinstance(param, float) and 0 < param <= 1:
                result = int(param * ref_value)
            else:
                result = int(param)

        # Garante valor mínimo se especificado
        if min_value is not None:
            result = max(min_value, result)
        else:
            result = max(1, result)  # Sempre pelo menos 1

        return result

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define um callback para relatar o progresso da execução do algoritmo.

        O callback é chamado em pontos-chave durante a execução para fornecer
        informações sobre o progresso atual. Isso permite implementar interfaces
        gráficas, logs detalhados ou sistemas de monitoramento.

        O callback recebe mensagens como:
        - "Criando população inicial..."
        - "Geração 50/1000, melhor=3"
        - "Timeout após 120.5s"
        - "Solução ótima encontrada!"

        Args:
            callback: Função que aceita uma string e não retorna nada.
                     Será chamada periodicamente durante a execução.

        Examples:
            >>> def my_progress(msg):
            ...     print(f"[BLF-GA] {msg}")
            >>> algorithm.set_progress_callback(my_progress)
        """
        self.progress_callback = callback

    def update_params(self, **params):
        """
        Atualiza parâmetros do algoritmo dinamicamente durante a execução.

        Esta função permite modificar configurações do BLF-GA em tempo de execução,
        incluindo parâmetros que requerem recálculo (como pop_size e initial_blocks).
        É útil para ajuste fino automático ou otimização de hiperparâmetros.

        Funcionalidades:
        - Atualiza parâmetros simples diretamente
        - Recalcula parâmetros dinâmicos (pop_size baseado em n)
        - Reinicializa estruturas dependentes (blocos)
        - Mantém consistência entre parâmetros relacionados

        Args:
            **params: Parâmetros a serem atualizados. Pode incluir qualquer
                     parâmetro válido do BLF-GA como pop_size, mut_prob,
                     cross_prob, initial_blocks, etc.

        Examples:
            >>> algorithm.update_params(mut_prob=0.05, cross_prob=0.9)
            >>> algorithm.update_params(pop_size=2.0)  # 2x o número de strings
            >>> algorithm.update_params(initial_blocks=0.2)  # 20% do comprimento
        """
        # Atualizar parâmetros básicos
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)

        # Reprocessar parâmetros dinâmicos se necessário
        # pop_size pode ser int fixo ou float (multiplicador do número de strings)
        if "pop_size" in params:
            min_pop_size = params.get("min_pop_size", self.min_pop_size)
            self.pop_size = self._resolve_dynamic_param(
                params["pop_size"], self.n, mode="multiplier", min_value=min_pop_size
            )

        # initial_blocks pode ser int fixo ou float (proporção do comprimento)
        if "initial_blocks" in params:
            self.initial_blocks = self._resolve_dynamic_param(
                params["initial_blocks"], self.L, mode="proportion"
            )
            # Reinicializar blocos com nova configuração
            self.blocks = self._initial_blocking()

        # no_improve_patience pode ser int fixo ou float (proporção de max_gens)
        if "no_improve_patience" in params:
            patience_param = params["no_improve_patience"]
            if isinstance(patience_param, float) and 0 < patience_param < 1:
                self.no_improve_patience = max(1, int(patience_param * self.max_gens))
            else:
                self.no_improve_patience = int(patience_param)

        # Reinicializar gerador de números aleatórios se semente mudou
        if "seed" in params:
            self.rng = random.Random(params["seed"])

    def set_history_callback(
        self, callback: Optional[Callable[[int, dict], None]]
    ) -> None:
        """
        Define o callback para registrar eventos dinâmicos durante a evolução.

        Args:
            callback: Função que recebe (generation, event_data) onde event_data
                     contém informações sobre operações dinâmicas como:
                     - disable_elitism_gens: quando elitismo é desabilitado
                     - adaptive_mutation: mudanças na taxa de mutação
                     - immigrant_injection: quando imigrantes são adicionados
                     - block_redivision: quando blocos são redivididos
                     - restart_mechanism: quando população é reiniciada
                     - early_stopping: quando algoritmo para por convergência
                     - elite_refinement: quando melhores soluções são refinadas
                     - adaptive_blocking: mudanças dinâmicas nos blocos
        """
        self.history_callback = callback

    def run(self) -> tuple[String, int, list]:
        """
        Executa o algoritmo BLF-GA para encontrar a string mais próxima.

        FLUXO DO ALGORITMO:

        1. **INICIALIZAÇÃO**:
           - Cria população inicial baseada em consenso + variações + aleatórias
           - Avalia fitness inicial (distância máxima para todas as strings)
           - Inicializa controles de adaptação e histórico

        2. **LOOP PRINCIPAL** (até critério de parada):
           a) **CONTROLE DE TEMPO**: Verifica limite de tempo
           b) **MECANISMOS ADAPTATIVOS**:
              - Gera imigrantes aleatórios para diversidade
              - Calcula diversidade populacional (distância de Hamming)
              - Ajusta mutação baseada em convergência/diversidade
           c) **APRENDIZADO POR BLOCOS**:
              - Extrai padrões consenso de cada bloco
              - Cria repositório de conhecimento local
           d) **EVOLUÇÃO GENÉTICA**:
              - Gera nova geração via seleção, crossover e mutação
              - Aplica elitismo adaptativo
              - Avalia fitness da nova população
           e) **REFINAMENTO LOCAL**:
              - Aplica busca local intensiva nos melhores indivíduos
              - Melhora qualidade das soluções elite
           f) **CONTROLES DE PARADA**:
              - Verifica early stopping por estagnação
              - Executa restart se necessário
              - Redivide blocos adaptativamente

        3. **FINALIZAÇÃO**:
           - Retorna melhor solução encontrada + fitness + histórico

        Returns:
            tuple[String, int, list]: (melhor_string, distância_máxima, histórico)
                - melhor_string: String center encontrada
                - distância_máxima: Fitness da melhor solução
                - histórico: Lista com fitness por geração
        """
        start = time.time()

        if self.progress_callback:
            self.progress_callback("Criando população inicial...")

        # === FASE 1: INICIALIZAÇÃO ===
        # Cria população inicial usando estratégia híbrida:
        # - Consenso das strings de entrada
        # - Variações do consenso com blocos das strings originais
        # - Strings completamente aleatórias
        pop = self._init_population()
        best = min(pop, key=lambda s: max_distance(s, self.strings))
        best_val = max_distance(best, self.strings)
        self.history = [best_val]
        no_improve = 0  # Contador para early stopping
        mut_prob_backup = self.mut_prob  # Backup para mutação adaptativa
        mut_adapt_timer = 0  # Timer para mutação adaptativa

        # === FASE 2: LOOP PRINCIPAL DE EVOLUÇÃO ===
        for gen in range(1, self.max_gens + 1):
            # --- CONTROLE DE TEMPO ---
            elapsed = time.time() - start
            if elapsed >= self.max_time:
                if self.progress_callback:
                    self.progress_callback(f"Timeout após {elapsed:.1f}s")
                break

            # Progresso via callback (que pode lançar exceção se cancelado)
            if self.progress_callback:
                progress_msg = f"Geração {gen}/{self.max_gens}, melhor={best_val}"
                self.progress_callback(progress_msg)

            # --- MECANISMO 1: IMIGRANTES ALEATÓRIOS ---
            # Injeta diversidade substituindo os piores indivíduos por strings aleatórias
            # Previne convergência prematura e mantém exploração
            if self.immigrant_freq and gen % self.immigrant_freq == 0:
                n_imm = int(self.immigrant_ratio * self.pop_size)
                for _ in range(n_imm):
                    rand_s = "".join(
                        self.rng.choice(self.alphabet) for _ in range(self.L)
                    )
                    pop[-(_ + 1)] = rand_s

                # Log da operação dinâmica
                if self.history_callback:
                    self.history_callback(
                        gen,
                        {
                            "event": "immigrant_injection",
                            "description": f"Injetados {n_imm} imigrantes na geração {gen} (a cada {self.immigrant_freq} gerações)",
                            "immigrant_count": n_imm,
                            "immigrant_ratio": self.immigrant_ratio,
                            "pop_size": self.pop_size,
                            "replaced_positions": list(
                                range(self.pop_size - n_imm, self.pop_size)
                            ),
                        },
                    )

            # --- MECANISMO 2: MUTAÇÃO ADAPTATIVA ---
            # Ajusta taxa de mutação baseada na diversidade populacional e convergência
            diversity = genetic_ops.mean_hamming_distance(pop)

            # Se diversidade baixa, aumenta mutação temporariamente
            if diversity < self.diversity_threshold * self.L:
                old_mut_prob = self.mut_prob
                self.mut_prob = mut_prob_backup * self.mutation_adapt_factor
                mut_adapt_timer = self.mutation_adapt_duration

                # Log da operação dinâmica
                if self.history_callback:
                    self.history_callback(
                        gen,
                        {
                            "event": "adaptive_mutation_diversity",
                            "description": f"Mutação aumentada por baixa diversidade: {old_mut_prob:.4f} → {self.mut_prob:.4f}",
                            "trigger": "low_diversity",
                            "diversity": diversity,
                            "threshold": self.diversity_threshold * self.L,
                            "old_mutation_rate": old_mut_prob,
                            "new_mutation_rate": self.mut_prob,
                            "adaptation_factor": self.mutation_adapt_factor,
                            "timer_duration": self.mutation_adapt_duration,
                        },
                    )

            # Se fitness estagnado por N gerações, aumenta mutação
            if len(self.history) > self.mutation_adapt_N and all(
                self.history[-i] == self.history[-1]
                for i in range(1, self.mutation_adapt_N + 1)
            ):
                old_mut_prob = self.mut_prob
                self.mut_prob = mut_prob_backup * self.mutation_adapt_factor
                mut_adapt_timer = self.mutation_adapt_duration

                # Log da operação dinâmica
                if self.history_callback:
                    self.history_callback(
                        gen,
                        {
                            "event": "adaptive_mutation_stagnation",
                            "description": f"Mutação aumentada por estagnação de {self.mutation_adapt_N} gerações: {old_mut_prob:.4f} → {self.mut_prob:.4f}",
                            "trigger": "fitness_stagnation",
                            "stagnation_gens": self.mutation_adapt_N,
                            "stagnant_fitness": self.history[-1],
                            "old_mutation_rate": old_mut_prob,
                            "new_mutation_rate": self.mut_prob,
                            "adaptation_factor": self.mutation_adapt_factor,
                            "timer_duration": self.mutation_adapt_duration,
                        },
                    )

            # Decrementa timer da mutação adaptativa
            if mut_adapt_timer > 0:
                mut_adapt_timer -= 1
                if mut_adapt_timer == 0:
                    old_mut_prob = self.mut_prob
                    self.mut_prob = mut_prob_backup

                    # Log da operação dinâmica
                    if self.history_callback:
                        self.history_callback(
                            gen,
                            {
                                "event": "adaptive_mutation_reset",
                                "description": f"Mutação voltou ao normal após período adaptativo: {old_mut_prob:.4f} → {self.mut_prob:.4f}",
                                "old_mutation_rate": old_mut_prob,
                                "normal_mutation_rate": self.mut_prob,
                                "adaptation_ended": True,
                            },
                        )

            # --- MECANISMO 3: APRENDIZADO POR BLOCOS ---
            # Extrai conhecimento local de cada bloco através de consenso
            # Cria repositório de padrões para guiar a evolução
            repo = self._learn_blocks(pop)

            # --- MECANISMO 4: EVOLUÇÃO GENÉTICA ---
            # Gera nova geração através de seleção, crossover e mutação
            pop = self._next_generation(pop, repo)

            # Ordenar população por fitness usando paralelismo
            pop = self._sort_population_parallel(pop)

            # --- MECANISMO 5: ELITISMO ADAPTATIVO ---
            # Preserva os melhores indivíduos, mas pode desabilitar periodicamente
            k = max(1, int(self.elite_rate * self.pop_size))
            if self.disable_elitism_gens and (gen % self.disable_elitism_gens == 0):
                elites = []  # Desabilita elitismo para aumentar diversidade
                # Log da operação dinâmica
                if self.history_callback:
                    self.history_callback(
                        gen,
                        {
                            "event": "disable_elitism",
                            "description": f"Elitismo desabilitado na geração {gen} (a cada {self.disable_elitism_gens} gerações)",
                            "elite_rate": self.elite_rate,
                            "pop_size": self.pop_size,
                            "expected_elites": k,
                        },
                    )
            else:
                elites = pop[:k]

            # --- MECANISMO 6: REFINAMENTO LOCAL ---
            # Aplica busca local intensiva nos melhores indivíduos
            if self.refine_elites == "all":
                old_elites_fitness = [max_distance(e, self.strings) for e in elites]
                pop[:k] = self._refine_elites(elites)  # Refina todos os elites
                new_elites_fitness = [
                    max_distance(pop[i], self.strings) for i in range(k)
                ]

                # Log da operação dinâmica
                if self.history_callback:
                    improvements = [
                        old - new
                        for old, new in zip(old_elites_fitness, new_elites_fitness)
                    ]
                    self.history_callback(
                        gen,
                        {
                            "event": "elite_refinement_all",
                            "description": f"Refinamento aplicado a todos os {k} elites",
                            "elite_count": k,
                            "refinement_type": self.refine_elites,
                            "old_fitness": old_elites_fitness,
                            "new_fitness": new_elites_fitness,
                            "improvements": improvements,
                            "total_improvement": sum(improvements),
                        },
                    )
            else:
                if elites:
                    old_best_fitness = max_distance(elites[0], self.strings)
                    pop[0] = self._refine_elites([elites[0]])[
                        0
                    ]  # Refina apenas o melhor
                    new_best_fitness = max_distance(pop[0], self.strings)

                    # Log da operação dinâmica
                    if self.history_callback:
                        self.history_callback(
                            gen,
                            {
                                "event": "elite_refinement_best",
                                "description": "Refinamento aplicado ao melhor elite",
                                "refinement_type": self.refine_elites,
                                "old_fitness": old_best_fitness,
                                "new_fitness": new_best_fitness,
                                "improvement": old_best_fitness - new_best_fitness,
                            },
                        )

            cur_best = pop[0]
            cur_val = max_distance(cur_best, self.strings)
            self.history.append(cur_val)

            # Atualiza melhor solução global
            if cur_val < best_val:
                best, best_val = cur_best, cur_val
                no_improve = 0
            else:
                no_improve += 1

            # --- MECANISMO 7: RESTART ---
            # Reinicia parte da população se estagnada por muito tempo
            if self.restart_patience and no_improve >= self.restart_patience:
                n_restart = int(self.restart_ratio * self.pop_size)
                restarted_positions = []
                for i in range(n_restart):
                    pos = -(i + 1)
                    pop[pos] = "".join(
                        self.rng.choice(self.alphabet) for _ in range(self.L)
                    )
                    restarted_positions.append(self.pop_size + pos)  # Posição absoluta

                # Log da operação dinâmica
                if self.history_callback:
                    self.history_callback(
                        gen,
                        {
                            "event": "restart_mechanism",
                            "description": f"Restart de {n_restart} indivíduos após {no_improve} gerações de estagnação",
                            "restarted_count": n_restart,
                            "restart_ratio": self.restart_ratio,
                            "stagnation_count": no_improve,
                            "restart_patience": self.restart_patience,
                            "restarted_positions": restarted_positions,
                            "pop_size": self.pop_size,
                        },
                    )

                no_improve = 0

            # --- CRITÉRIO DE PARADA 1: EARLY STOPPING ---
            # Encerra se não houver melhoria por X gerações
            if self.no_improve_patience and no_improve >= self.no_improve_patience:
                if self.progress_callback:
                    self.progress_callback(
                        f"Encerrando por early stopping após {no_improve} gerações sem melhoria."
                    )

                # Log da operação dinâmica
                if self.history_callback:
                    self.history_callback(
                        gen,
                        {
                            "event": "early_stopping",
                            "description": f"Algoritmo parado por early stopping após {no_improve} gerações sem melhoria",
                            "no_improve_count": no_improve,
                            "patience_threshold": self.no_improve_patience,
                            "current_generation": gen,
                            "best_fitness": best_val,
                            "reason": "fitness_stagnation",
                        },
                    )
                break

            # Log apenas a cada 50 gerações ou na última
            if gen % 50 == 0 or gen == self.max_gens or best_val == 0:
                logger.debug(
                    "Geração %s: melhor_dist=%s, diversidade=%.2f",
                    gen,
                    best_val,
                    diversity,
                )

            # --- CRITÉRIO DE PARADA 2: SOLUÇÃO ÓTIMA ---
            # Encerra se encontrar solução perfeita (distância 0)
            if best_val == 0:
                if self.progress_callback:
                    self.progress_callback("Solução ótima encontrada!")
                break

            # --- MECANISMO 8: REDIVISÃO ADAPTATIVA ---
            # Redefine blocos baseado na entropia das posições
            # Permite adaptação dinâmica da estrutura de blocos
            if gen % self.rediv_freq == 0:
                old_blocks = self.blocks.copy()
                self.blocks = self._adaptive_blocking(pop)

                # Log da operação dinâmica
                if self.history_callback:
                    self.history_callback(
                        gen,
                        {
                            "event": "block_redivision",
                            "description": f"Blocos redivididos na geração {gen} (a cada {self.rediv_freq} gerações)",
                            "old_blocks": old_blocks,
                            "new_blocks": self.blocks,
                            "old_block_count": len(old_blocks),
                            "new_block_count": len(self.blocks),
                            "redivision_frequency": self.rediv_freq,
                        },
                    )

        # === FASE 3: FINALIZAÇÃO ===
        return best, best_val, self.history

    def _init_population(self) -> Population:
        """
        Cria população inicial usando estratégia híbrida inteligente.

        A inicialização é crucial para o sucesso do algoritmo genético. Esta função
        implementa uma estratégia sofisticada que combina qualidade inicial com
        diversidade populacional para evitar convergência prematura.

        ESTRATÉGIA DE INICIALIZAÇÃO HÍBRIDA:

        1. **CONSENSO GLOBAL** (1 indivíduo):
           - Calcula string consenso baseada na moda de cada posição
           - Garante que pelo menos um indivíduo tenha qualidade inicial alta
           - Funciona como "semente" de boa qualidade para a evolução

        2. **VARIAÇÕES INTELIGENTES** (~1/3 da população):
           - Parte do consenso e modifica blocos específicos
           - Para cada bloco, com 50% de chance substitui por bloco correspondente
             de uma string original aleatória
           - Mantém estrutura global boa mas introduz diversidade local
           - Explora combinações de padrões locais das strings originais

        3. **DIVERSIDADE ALEATÓRIA** (restante da população):
           - Strings completamente aleatórias para máxima diversidade
           - Previne viés inicial e garante exploração ampla do espaço
           - Funciona como fonte de "material genético" novo

        VANTAGENS DESTA ABORDAGEM:
        - Qualidade inicial: Consenso fornece ponto de partida próximo ao ótimo
        - Diversidade estrutural: Variações exploram combinações de padrões reais
        - Exploração ampla: Aleatoriedade evita convergência prematura
        - Eficiência: Reduz gerações necessárias comparado à inicialização puramente aleatória

        Returns:
            Population: Lista de strings representando a população inicial ordenada
                       implicitamente por qualidade (consenso primeiro)

        Example:
            Para strings ["ACGT", "AGCT", "ATCT"] com blocos [(0,2), (2,4)]:
            - Consenso: "ACCT" (A=3/3, C=2/3, C=2/3, T=3/3)
            - Variação 1: "ACGT" (bloco 2 de "ACGT")
            - Variação 2: "AGCT" (bloco 1 de "AGCT")
            - Aleatórias: "TTAA", "CGAT", etc.
        """
        # 1. CONSENSO GLOBAL: Criar string baseada na moda de cada posição
        # Counter conta frequência de cada símbolo, most_common(1) pega o mais frequente
        consensus = "".join(
            Counter(pos).most_common(1)[0][0] for pos in zip(*self.strings)
        )
        pop = [consensus]

        # 2. VARIAÇÕES INTELIGENTES: Modificar consenso por blocos
        for _ in range(self.pop_size // 3):
            s = list(consensus)  # Começar com consenso
            # Para cada bloco definido na divisão atual
            for l, r in self.blocks:
                if self.rng.random() < 0.5:  # 50% chance de modificar bloco
                    # Escolher string original aleatória como fonte
                    src = self.rng.choice(self.strings)
                    # Substituir bloco por bloco correspondente da string fonte
                    s[l:r] = src[l:r]
            pop.append("".join(s))

        # 3. DIVERSIDADE ALEATÓRIA: Preencher restante com strings aleatórias
        while len(pop) < self.pop_size:
            # Gerar string completamente aleatória
            rand_s = "".join(self.rng.choice(self.alphabet) for _ in range(self.L))
            pop.append(rand_s)

        return pop

    def _initial_blocking(self) -> list[tuple[int, int]]:
        """
        Cria divisão inicial das strings em blocos de tamanho aproximadamente uniforme.

        A divisão em blocos é fundamental para o componente "Blockwise Learning"
        do BLF-GA. Blocos permitem que o algoritmo aprenda e otimize padrões
        locais independentemente, melhorando a eficiência da busca.

        ALGORITMO DE DIVISÃO:
        1. Calcula tamanho base do bloco respeitando tamanho mínimo
        2. Divide string em blocos consecutivos de tamanho uniforme
        3. Último bloco pode ser menor se L não for divisível

        DESIGN CHOICES:
        - Blocos contíguos preservam localidade espacial
        - Tamanho uniforme simplifica operações genéticas
        - Tamanho mínimo evita blocos muito pequenos (ineficientes)

        Returns:
            list[tuple[int, int]]: Lista de tuplas (início, fim) definindo os blocos.
                                  Início é inclusivo, fim é exclusivo (padrão Python).

        Example:
            Para L=10, initial_blocks=3, min_block_len=2:
            - size = max(2, 10//3) = max(2, 3) = 3
            - Blocos: [(0,3), (3,6), (6,9), (9,10)]
            - Resultado: 4 blocos de tamanhos [3,3,3,1]
        """
        # Calcular tamanho base respeitando mínimo
        size = max(self.min_block_len, self.L // self.initial_blocks)

        # Criar blocos consecutivos
        blocks = []
        for i in range(0, self.L, size):
            # Fim é min(i + size, L) para não ultrapassar comprimento
            blocks.append((i, min(i + size, self.L)))

        return blocks

    def _learn_blocks(self, pop: Population) -> list[String]:
        """
        Extrai conhecimento local de cada bloco através de aprendizado por consenso.

        Esta é a implementação do componente "LEARNING" do BLF-GA. A função
        analisa a população atual para identificar padrões locais eficazes
        em cada bloco, criando um repositório de conhecimento que pode ser
        usado para guiar a evolução.

        ALGORITMO DE APRENDIZADO:

        1. **SEGMENTAÇÃO**: Para cada bloco, extrai segmentos correspondentes
           de todos os indivíduos da população

        2. **ANÁLISE DE CONSENSO**: Calcula consenso posição-a-posição dentro
           do bloco usando votação por maioria

        3. **CONSTRUÇÃO DO REPOSITÓRIO**: Armazena padrões aprendidos para
           uso em operações genéticas futuras

        VANTAGENS DO APRENDIZADO POR BLOCOS:
        - **Localidade**: Identifica padrões em regiões específicas
        - **Adaptabilidade**: Evolui conforme população melhora
        - **Eficiência**: Guia busca com informação estruturada
        - **Hibridização**: Combina consenso local com evolução global

        EXEMPLO DE FUNCIONAMENTO:
        Para bloco (0,2) em população ["AC...", "AG...", "AT..."]:
        - Segmentos: ["AC", "AG", "AT"]
        - Posição 0: A=3/3 → consenso='A'
        - Posição 1: C=1/3, G=1/3, T=1/3 → consenso='C' (primeiro mais comum)
        - Padrão aprendido: "AC"

        Args:
            pop: População atual para análise de padrões

        Returns:
            list[String]: Repositório de padrões consenso, um por bloco.
                         Índice corresponde ao índice do bloco em self.blocks.

        Note:
            O repositório pode ser usado posteriormente em operações como
            crossover_blend_blocks ou para guiar mutações inteligentes.
        """
        repo = []

        # Processar cada bloco independentemente
        for l, r in self.blocks:
            # 1. SEGMENTAÇÃO: Coletar segmentos do bloco atual
            # Filtrar indivíduos com comprimento suficiente para evitar erros
            segments = [ind[l:r] for ind in pop if len(ind) >= r]

            if not segments:  # Se nenhum segmento válido
                # Criar padrão aleatório como fallback
                block_pattern = "".join(
                    self.rng.choice(self.alphabet) for _ in range(r - l)
                )
            else:
                # 2. ANÁLISE DE CONSENSO: Votar posição por posição
                # zip(*segments) transpõe matriz: de [seg1, seg2] para [pos1_chars, pos2_chars]
                consensus_chars = []
                for position_chars in zip(*segments):
                    # Contar frequência e pegar mais comum
                    char_counter = Counter(position_chars)
                    most_common_char = char_counter.most_common(1)[0][0]
                    consensus_chars.append(most_common_char)

                # 3. CONSTRUIR PADRÃO
                block_pattern = "".join(consensus_chars)

            repo.append(block_pattern)

        return repo

    def _tournament_selection(self, pop: Population, k: int = 2) -> String:
        """
        Seleciona um indivíduo da população através de seleção por torneio.

        A seleção por torneio é um dos métodos mais eficazes em algoritmos genéticos
        pois mantém pressão seletiva moderada e permite controle fino através do
        parâmetro k. É especialmente adequada para o CSP pois preserva diversidade
        enquanto favorece soluções de melhor qualidade.

        ALGORITMO DE SELEÇÃO POR TORNEIO:

        1. **AMOSTRAGEM ALEATÓRIA**:
           - Escolhe k indivíduos aleatoriamente da população (sem reposição)
           - Evita viés de posição e garante chance igual a todos

        2. **COMPETIÇÃO LOCAL**:
           - Avalia fitness (distância máxima) de cada participante
           - Menor distância = melhor fitness = maior chance de reprodução

        3. **SELEÇÃO DO VENCEDOR**:
           - Retorna o indivíduo com menor distância máxima do torneio
           - Garante que melhores soluções tenham maior probabilidade de reprodução

        CARACTERÍSTICAS E VANTAGENS:
        - **Pressão Seletiva Ajustável**: k maior = pressão maior
        - **Preservação de Diversidade**: Não elimina completamente soluções piores
        - **Simplicidade Computacional**: O(k) comparações por seleção
        - **Robustez**: Funciona bem independente da distribuição de fitness

        CONFIGURAÇÃO RECOMENDADA:
        - k=2: Pressão moderada, boa para exploration/exploitation
        - k=3-5: Pressão maior, acelera convergência
        - k>5: Pressão muito alta, risco de convergência prematura

        Args:
            pop: População atual para seleção de pais
            k: Tamanho do torneio (número de competidores). Se None, usa
               self.tournament_k ou 2 como padrão.

        Returns:
            String: Indivíduo vencedor do torneio, adequado para reprodução

        Example:
            Em população ["AAAA", "TTTT", "CCCC"] com fitnesses [5, 3, 8]:
            - Para k=2, escolhe 2 aleatórios, ex: ["TTTT", "CCCC"]
            - Compara fitness: TTTT(3) < CCCC(8)
            - Retorna "TTTT" como vencedor
        """
        # Determina tamanho do torneio usando hierarquia de parâmetros
        if k is None:
            k = self.tournament_k if self.tournament_k is not None else 2

        # 1. AMOSTRAGEM: Escolher k competidores aleatórios
        # sample() garante que não há repetições no torneio
        tournament_pool = self.rng.sample(pop, k)

        # 2. COMPETIÇÃO: Avaliar todos os competidores
        # min() com key= encontra o competidor com menor distância máxima
        # Menor distância = melhor solução para o CSP
        winner = min(tournament_pool, key=lambda s: max_distance(s, self.strings))

        return winner

    def _apply_crossover(self, p1, p2):
        """
        Aplica operação de crossover entre dois pais para gerar descendentes.

        O crossover é o operador principal de recombinação genética no BLF-GA,
        responsável por combinar características dos pais para criar filhos
        potencialmente melhores. Diferentes tipos de crossover exploram
        diferentes estratégias de combinação genética.

        TIPOS DE CROSSOVER IMPLEMENTADOS:

        1. **ONE_POINT** (Clássico):
           - Escolhe ponto de corte aleatório
           - Troca segmentos entre pais a partir do ponto
           - Simples e eficaz para problemas de permutação
           - Preserva blocos grandes de características

        2. **UNIFORM** (Gene-a-gene):
           - Para cada posição, escolhe aleatoriamente gene do pai1 ou pai2
           - Máxima mistura genética possível
           - Melhor exploração, mas pode quebrar bons padrões
           - Recomendado quando não há estrutura clara no problema

        3. **BLEND_BLOCKS** (Específico do BLF-GA):
           - Troca blocos inteiros entre os pais
           - Preserva padrões locais aprendidos
           - Alinhado com a filosofia de aprendizado por blocos
           - Combina qualidade local com diversidade global

        DESIGN ADAPTADO PARA CSP:
        - Preserva alinhamento de posições (importante para distância)
        - Mantém comprimento constante das strings
        - Explora combinações de padrões locais eficazes

        Args:
            p1: Primeiro pai (string representando candidato)
            p2: Segundo pai (string representando candidato)

        Returns:
            tuple: (filho1, filho2) - Dois descendentes gerados pela recombinação

        Examples:
            ONE_POINT com ponto=2:
            p1="AAAA", p2="TTTT" → ("AATT", "TTAA")

            UNIFORM com máscara=[1,0,1,0]:
            p1="AAAA", p2="TTTT" → ("ATAT", "TATA")

            BLEND_BLOCKS com blocos=[(0,2),(2,4)]:
            p1="AATT", p2="CCGG" → ("AAGG", "CCTT")
        """
        # Aplicar crossover baseado no tipo configurado
        if self.crossover_type == "one_point":
            # Crossover de ponto único - preserva grandes segmentos
            c1, c2 = genetic_ops.crossover_one_point(p1, p2, self.rng)
        elif self.crossover_type == "uniform":
            # Crossover uniforme - máxima recombinação gene-a-gene
            c1, c2 = genetic_ops.crossover_uniform(p1, p2, self.rng)
        elif self.crossover_type == "blend_blocks":
            # Crossover por blocos - preserva padrões locais aprendidos
            # Usa self.blocks para determinar divisões de troca
            c1, c2 = genetic_ops.crossover_blend_blocks(p1, p2, self.blocks, self.rng)
        else:
            # Fallback: retorna pais inalterados se tipo inválido
            c1, c2 = p1, p2

        return c1, c2

    def _apply_mutation(self, ind):
        """
        Aplica operação de mutação a um indivíduo para introduzir variabilidade.

        A mutação é essencial para manter diversidade genética e explorar
        novas regiões do espaço de busca. No contexto do CSP, diferentes
        tipos de mutação oferecem estratégias complementares de refinamento
        e exploração local.

        TIPOS DE MUTAÇÃO IMPLEMENTADOS:

        1. **MULTI** (Mutação Multi-Ponto):
           - Altera N posições aleatórias simultaneamente
           - Explora vizinhança local ampla
           - Controlado por mutation_multi_n
           - Eficaz para escape de ótimos locais

        2. **INVERSION** (Inversão de Segmento):
           - Seleciona segmento aleatório e inverte sua ordem
           - Reorganiza estrutura local preservando composição
           - Útil quando ordem das posições importa
           - Mantém distribuição de símbolos

        3. **TRANSPOSITION** (Transposição de Segmentos):
           - Escolhe dois segmentos e troca suas posições
           - Realoca padrões para novas posições
           - Explora permutações estruturais
           - Preserva padrões locais em novas posições

        FILOSOFIA DE MUTAÇÃO NO BLF-GA:
        - **Exploração Controlada**: Taxa ajustável dinamicamente
        - **Preservação Estrutural**: Tipos que mantêm padrões úteis
        - **Escape Local**: Capacidade de sair de ótimos locais
        - **Complementaridade**: Diferentes tipos para diferentes situações

        TAXA DE MUTAÇÃO ADAPTATIVA:
        A taxa é ajustada dinamicamente baseada em:
        - Diversidade populacional atual
        - Histórico de stagnação
        - Progresso da evolução

        Args:
            ind: Indivíduo (string) a ser mutado

        Returns:
            String: Indivíduo mutado com possíveis alterações

        Examples:
            MULTI com n=2: "AAAA" → "ATAA" (posições 1,3 alteradas)
            INVERSION: "ABCD" → "ACBD" (segmento [1:3] invertido)
            TRANSPOSITION: "ABCD" → "CDAB" (segmentos [0:2] e [2:4] trocados)
        """
        # Aplicar mutação baseada no tipo configurado
        if self.mutation_type == "multi":
            # Mutação multi-ponto: altera N posições aleatórias
            # Parâmetro mutation_multi_n controla intensidade
            return genetic_ops.mutate_multi(
                ind, self.alphabet, self.rng, n=self.mutation_multi_n
            )
        elif self.mutation_type == "inversion":
            # Mutação por inversão: inverte segmento aleatório
            # Preserva composição mas reorganiza estrutura
            return genetic_ops.mutate_inversion(ind, self.rng)
        elif self.mutation_type == "transposition":
            # Mutação por transposição: troca dois segmentos
            # Explora reorganizações estruturais mantendo padrões
            return genetic_ops.mutate_transposition(ind, self.rng)
        else:
            # Fallback: retorna indivíduo inalterado se tipo inválido
            return ind

    def _apply_refinement(self, ind):
        """
        Aplica refinamento local (busca local intensiva) a um indivíduo.

        O refinamento é uma componente crucial dos algoritmos híbridos genético-busca
        local. No BLF-GA, serve para intensificar a busca em torno de soluções
        promissoras, melhorando a qualidade final e acelerando a convergência
        para ótimos locais de alta qualidade.

        TIPOS DE REFINAMENTO IMPLEMENTADOS:

        1. **GREEDY** (Busca Gulosa):
           - Testa cada posição sistematicamente
           - Para cada posição, tenta todos os símbolos do alfabeto
           - Aceita primeira melhoria encontrada
           - Determinístico e eficiente para melhorias locais

        2. **SWAP** (Troca de Posições):
           - Testa todas as trocas possíveis entre pares de posições
           - Explora reorganizações que preservam composição
           - Útil quando distribuição de símbolos está correta
           - Complexidade O(L²) onde L é comprimento da string

        3. **INSERTION** (Inserção de Segmentos):
           - Remove segmento de uma posição e insere em outra
           - Testa reposicionamento de padrões locais
           - Explora diferentes organizações estruturais
           - Mantém todos os símbolos originais

        4. **2OPT** (Otimização 2-opt):
           - Remove duas arestas e reconecta de forma diferente
           - Clássica de otimização combinatória
           - Elimina cruzamentos em representações espaciais
           - Convergência rápida para ótimos 2-opt

        ESTRATÉGIA DE APLICAÇÃO:
        - **Seletiva**: Aplicada apenas aos melhores indivíduos (elites)
        - **Limitada**: Controlada por refine_iter_limit para evitar custo excessivo
        - **Complementar**: Combina com evolução genética global

        CONFIGURAÇÃO DE USO:
        - refine_elites="all": Refina todos os indivíduos elite
        - refine_elites="best": Refina apenas o melhor indivíduo
        - Refinamento balanceia intensificação vs. custo computacional

        Args:
            ind: Indivíduo (string) a ser refinado localmente

        Returns:
            String: Indivíduo refinado com melhorias locais aplicadas

        Note:
            O refinamento pode não melhorar o indivíduo se já estiver
            em um ótimo local para o tipo de busca aplicado.
        """
        # Aplicar refinamento baseado no tipo configurado
        if self.refinement_type == "greedy":
            # Busca gulosa posição-a-posição
            # Testa todos os símbolos em cada posição sistematicamente
            return genetic_ops.refine_greedy(ind, self.strings)
        elif self.refinement_type == "swap":
            # Refinamento por troca de posições
            # Explora todas as permutações de pares de posições
            return genetic_ops.refine_swap(ind, self.strings)
        elif self.refinement_type == "insertion":
            # Refinamento por inserção de segmentos
            # Testa reposicionamento de blocos para novas posições
            return genetic_ops.refine_insertion(ind, self.strings)
        elif self.refinement_type == "2opt":
            # Otimização 2-opt clássica
            # Remove e reconecta arestas para eliminar cruzamentos
            return genetic_ops.refine_2opt(ind, self.strings)
        else:
            # Fallback: retorna indivíduo inalterado se tipo inválido
            return ind

    def _next_generation(
        self, pop: Population, repo: list[String]
    ) -> Population:  # pylint: disable=unused-argument
        """
        Gera nova população através de operações genéticas clássicas.

        Esta função implementa o core do algoritmo genético, coordenando
        os operadores de seleção, crossover e mutação para produzir uma
        nova geração de candidatos. É o motor principal da evolução no BLF-GA.

        ALGORITMO DE GERAÇÃO (Modelo Geracional):

        1. **AVALIAÇÃO E ORDENAÇÃO**:
           - Calcula fitness de todos os indivíduos (paralelizado)
           - Ordena população por qualidade (melhor → pior)
           - Usa max_distance como função de fitness

        2. **ELITISMO ADAPTATIVO**:
           - Preserva elite_rate % dos melhores indivíduos automaticamente
           - Garante que avanços não sejam perdidos
           - Pode ser desabilitado periodicamente para aumentar diversidade

        3. **REPRODUÇÃO PROBABILÍSTICA**:
           - Seleciona pais via torneio (pressão seletiva balanceada)
           - Aplica crossover com probabilidade cross_prob
           - Aplica mutação com probabilidade mut_prob
           - Gera filhos até completar população alvo

        4. **CONTROLE DE DIVERSIDADE** (Opcional):
           - Aplica niching se habilitado para evitar convergência prematura
           - Mantém distância mínima entre indivíduos
           - Preserva múltiplas regiões promissoras do espaço de busca

        CARACTERÍSTICAS DO MODELO:
        - **Sobreposição Geracional**: Elites sempre preservados
        - **Seleção Duplicada**: Mesmos pais podem gerar múltiplos filhos
        - **Operadores Independentes**: Crossover e mutação aplicados sequencialmente
        - **Tamanho Fixo**: População mantém tamanho constante

        BALANCEAMENTO EXPLORATION/EXPLOITATION:
        - Seleção por torneio: Favorece melhores mas preserva diversidade
        - Crossover probabilístico: Combina características prometedoras
        - Mutação: Introduz novidade e escape de ótimos locais
        - Elitismo: Preserva progresso já alcançado

        Args:
            pop: População atual a ser evoluída
            repo: Repositório de conhecimento por blocos (não usado nesta versão,
                  mas mantido para compatibilidade com versões futuras)

        Returns:
            Population: Nova população gerada com indivíduos evoluídos

        Note:
            O parâmetro 'repo' está marcado como unused pois esta implementação
            não usa diretamente o conhecimento de blocos na geração. Versões
            futuras podem integrar repo nos operadores genéticos.
        """
        # 1. AVALIAÇÃO E ORDENAÇÃO PARALELA
        # Calcula fitness de todos os indivíduos usando múltiplas threads
        # Ordena do melhor (menor distância) para o pior (maior distância)
        evaluated_pop = self._evaluate_population_parallel(pop)
        pop_sorted = [s for s, _ in evaluated_pop]

        # 2. ELITISMO - PRESERVAÇÃO DOS MELHORES
        # Calcula número de elites baseado na taxa configurada
        elite_n = max(1, int(self.elite_rate * self.pop_size))
        # Inicializa nova população com os melhores indivíduos da geração anterior
        new_pop = pop_sorted[:elite_n]

        # 3. REPRODUÇÃO - GERAÇÃO DE DESCENDENTES
        # Completa população através de reprodução sexual
        while len(new_pop) < self.pop_size:
            # SELEÇÃO: Escolhe dois pais via torneio
            p1 = self._tournament_selection(pop_sorted, k=self.tournament_k)
            p2 = self._tournament_selection(pop_sorted, k=self.tournament_k)

            # CROSSOVER: Aplica recombinação baseado em probabilidade
            if self.rng.random() < self.cross_prob:
                c1, c2 = self._apply_crossover(p1, p2)
            else:
                # Se não houver crossover, filhos são cópias dos pais
                c1, c2 = p1, p2

            # MUTAÇÃO: Aplica mutação a ambos os filhos
            # Probabilidade de mutação é aplicada internamente nos operadores
            c1 = self._apply_mutation(c1)
            c2 = self._apply_mutation(c2)

            # INSERÇÃO: Adiciona filhos à nova população
            new_pop.append(c1)
            if len(new_pop) < self.pop_size:  # Evita ultrapassar tamanho alvo
                new_pop.append(c2)

        # 4. CONTROLE DE DIVERSIDADE (NICHING OPCIONAL)
        # Aplica niching se habilitado para preservar diversidade local
        if self.niching:
            new_pop = self._apply_niching(new_pop)

        # 5. GARANTIA DE TAMANHO
        # Trunca população para tamanho exato (precaução contra bugs)
        return new_pop[: self.pop_size]

    def _refine_elites(self, pop: Population) -> Population:
        """
        Aplica refinamento local intensivo aos indivíduos elite da população.

        Esta função implementa a hibridização genético-busca local específica
        para os melhores indivíduos, combinando a capacidade de exploração
        global do algoritmo genético com a intensificação local da busca
        gulosa. É uma das características distintivas do BLF-GA.

        ESTRATÉGIA DE REFINAMENTO ELITE:

        1. **SELEÇÃO SELETIVA**:
           - Aplica refinamento apenas aos melhores indivíduos
           - Evita desperdício computacional em soluções de baixa qualidade
           - Concentra esforço onde há maior potencial de melhoria

        2. **PARALELIZAÇÃO EFICIENTE**:
           - Usa ThreadPoolExecutor para refinamento simultâneo
           - Aproveita múltiplos cores para acelerar busca local
           - Mantém eficiência mesmo com população grande

        3. **INTEGRAÇÃO BALANCEADA**:
           - Balança custo computacional vs. qualidade de solução
           - Evita dominância da busca local sobre evolução genética
           - Preserva diversidade populacional

        FILOSOFIA HÍBRIDA:
        - **Global + Local**: AG explora, busca local refina
        - **Qualidade + Diversidade**: Melhora elites sem homogeneizar
        - **Eficiência + Eficácia**: Paralelização sem overhead excessivo

        CONFIGURAÇÃO ADAPTATIVA:
        - Frequência controlada por refine_elites ("all", "best", ou None)
        - Intensidade limitada por refine_iter_limit
        - Tipo de refinamento configurável (greedy, swap, insertion, 2opt)

        Args:
            pop: População de indivíduos elite a serem refinados

        Returns:
            Population: População refinada com melhorias locais aplicadas

        Note:
            O refinamento é aplicado de forma independente a cada indivíduo,
            permitindo paralelização eficiente sem conflitos de recursos.
        """

        # Função local para refinamento individual
        # Encapsula lógica de refinamento para uso com ThreadPoolExecutor
        def refine(ind):
            """Aplica refinamento local a um indivíduo específico."""
            return self._apply_refinement(ind)

        # PARALELIZAÇÃO DO REFINAMENTO
        # Usa ThreadPoolExecutor para processar múltiplos indivíduos simultaneamente
        # Cada thread refina um indivíduo independentemente
        with ThreadPoolExecutor() as executor:
            # map() aplica função refine a cada elemento de pop
            # list() coleta resultados mantendo ordem original
            refined: Population = list(executor.map(refine, pop))

        return refined

    def _apply_niching(self, pop: Population) -> Population:
        """
        Aplica niching para preservar diversidade local e evitar convergência prematura.

        O niching é uma técnica avançada de diversidade que mantém múltiplas
        subpopulações em diferentes regiões do espaço de busca. No contexto
        do CSP, ajuda a preservar diferentes "famílias" de soluções que podem
        convergir para diferentes ótimos locais de alta qualidade.

        ALGORITMO DE NICHING POR DISTÂNCIA:

        1. **ORDENAÇÃO POR FITNESS**:
           - Organiza indivíduos do melhor para o pior
           - Garante que melhores soluções tenham prioridade na seleção
           - Base para seleção gulosa com restrição de distância

        2. **SELEÇÃO COM RESTRIÇÃO DE DISTÂNCIA**:
           - Itera pelos indivíduos em ordem de qualidade
           - Aceita indivíduo apenas se estiver suficientemente distante
           - Distância medida por hamming_dist entre strings
           - Threshold controlado por niching_radius

        3. **PREENCHIMENTO ADAPTATIVO**:
           - Relaxa restrições se população fica pequena demais
           - Adiciona indivíduos aleatórios se necessário
           - Garante manutenção do tamanho populacional

        VANTAGENS DO NICHING:
        - **Múltiplos Ótimos**: Preserva diferentes soluções de alta qualidade
        - **Diversidade Estrutural**: Evita homogeneização da população
        - **Exploração Paralela**: Mantém busca em múltiplas regiões
        - **Robustez**: Reduz risco de convergência prematura

        CONFIGURAÇÃO DE PARÂMETROS:
        - niching_radius pequeno: Maior diversidade, população mais esparsa
        - niching_radius grande: Menor diversidade, convergência mais rápida
        - Balanço exploration/exploitation através do raio

        APLICAÇÃO NO CSP:
        - Preserva strings com diferentes padrões estruturais
        - Mantém diversidade em regiões promissoras do espaço
        - Evita que população colapse em único ótimo local

        Args:
            pop: População a ser processada com niching

        Returns:
            Population: População com niching aplicado, mantendo diversidade
                       local e tamanho original

        Note:
            O niching pode reduzir taxa de convergência em troca de maior
            diversidade e melhor cobertura do espaço de busca.
        """
        # Validação de precondições
        if not self.niching or len(pop) <= 1:
            return pop

        # 1. ORDENAÇÃO POR FITNESS
        # Cria lista de tuplas (indivíduo, fitness) para seleção eficiente
        pop_with_fitness = [(ind, max_distance(ind, self.strings)) for ind in pop]
        # Ordena do melhor (menor distância) para o pior (maior distância)
        pop_with_fitness.sort(key=lambda x: x[1])

        # 2. SELEÇÃO COM RESTRIÇÃO DE DISTÂNCIA
        niched_pop = []

        for ind, fitness in pop_with_fitness:
            # Verifica distância para todos os indivíduos já aceitos
            too_close = False
            for accepted_ind in niched_pop:
                # Calcula distância de Hamming entre strings
                if hamming_dist(ind, accepted_ind) < self.niching_radius:
                    too_close = True
                    break

            # CRITÉRIO DE ACEITAÇÃO:
            # - Aceita se distante o suficiente de todos os aceitos
            # - OU se população ainda está pequena (relaxa restrições)
            if not too_close or len(niched_pop) < self.pop_size // 2:
                niched_pop.append(ind)

            # CONTROLE DE TAMANHO: Para quando tiver indivíduos suficientes
            if len(niched_pop) >= self.pop_size:
                break

        # 3. PREENCHIMENTO COM DIVERSIDADE ALEATÓRIA
        # Se niching removeu muitos indivíduos, adiciona diversidade aleatória
        while len(niched_pop) < self.pop_size:
            # Gera string completamente aleatória
            rand_s = "".join(self.rng.choice(self.alphabet) for _ in range(self.L))
            niched_pop.append(rand_s)

        return niched_pop

    def _adaptive_blocking(self, pop: Population) -> list[tuple[int, int]]:
        """
        Redivide strings em blocos adaptativos baseado na entropia de cada posição.

        Esta é uma das características mais avançadas do BLF-GA: a capacidade
        de adaptar dinamicamente a estrutura de blocos durante a evolução.
        A redivisão é baseada na análise de entropia, permitindo que o algoritmo
        identifique regiões com diferentes graus de consenso e variabilidade.

        ALGORITMO DE REDIVISÃO ADAPTATIVA:

        1. **ANÁLISE DE ENTROPIA POSICIONAL**:
           - Para cada posição i, calcula entropia H(i) = -Σ p(c) * log₂(p(c))
           - p(c) é probabilidade do símbolo c na posição i
           - Entropia alta = muita variação, Entropia baixa = consenso

        2. **THRESHOLD ADAPTATIVO**:
           - Define limiar como 70% da entropia máxima observada
           - Permite adaptação automática ao grau de diversidade atual
           - Evita hardcoding de valores fixos que podem ser inadequados

        3. **ESTRATÉGIA DE DIVISÃO INTELIGENTE**:
           - Posições com alta entropia (>threshold) → blocos menores
           - Posições com baixa entropia (≤threshold) → blocos maiores
           - Balanceia granularidade de aprendizado com eficiência

        INTUIÇÃO ALGORÍTMICA:
        - **Alta Entropia**: Região ainda em evolução, precisa blocos pequenos
          para capturar detalhes e permitir recombinação fina
        - **Baixa Entropia**: Região convergindo, pode usar blocos maiores
          para preservar padrões consenso já estabelecidos

        BENEFÍCIOS DA ADAPTAÇÃO:
        - **Evolução Estrutural**: Blocos evoluem com a população
        - **Eficiência Dinâmica**: Foco computacional em regiões que precisam
        - **Convergência Inteligente**: Acelera convergência preservando qualidade
        - **Robustez**: Adapta-se automaticamente a diferentes tipos de problema

        EXEMPLO DE FUNCIONAMENTO:
        Para posições com entropias [0.1, 0.8, 0.9, 0.2, 0.1] e threshold=0.6:
        - Pos 0: 0.1 ≤ 0.6 → bloco grande (conservada)
        - Pos 1: 0.8 > 0.6 → bloco pequeno (em evolução)
        - Pos 2: 0.9 > 0.6 → bloco pequeno (em evolução)
        - Pos 3: 0.2 ≤ 0.6 → bloco grande (conservada)
        - Pos 4: 0.1 ≤ 0.6 → bloco grande (conservada)

        Args:
            pop: População atual para análise de entropia

        Returns:
            list[tuple[int, int]]: Nova divisão em blocos adaptada ao estado
                                  evolutivo atual da população

        Note:
            Esta função é chamada periodicamente (rediv_freq) para manter
            a estrutura de blocos sincronizada com o progresso evolutivo.
        """
        # 1. ANÁLISE DE ENTROPIA POSICIONAL
        # Inicializa array para armazenar entropia de cada posição
        ent = np.zeros(self.L)

        for pos in range(self.L):
            # Coleta símbolos válidos na posição atual
            # Filtra strings que são longas suficientes para ter essa posição
            valid_chars = [ind[pos] for ind in pop if pos < len(ind)]

            if valid_chars:  # Se temos pelo menos um caractere válido
                # Conta frequência de cada símbolo
                cnt = Counter(valid_chars)
                # Calcula probabilidades normalizadas
                probs = np.array(list(cnt.values())) / len(valid_chars)
                # Calcula entropia: H = -Σ p * log₂(p)
                # Usa logaritmo base 2 para entropia em bits
                ent[pos] = -np.sum(probs * np.log2(probs))

        # 2. THRESHOLD ADAPTATIVO
        # Define limiar como 70% da entropia máxima observada
        # Se entropia máxima é 0 (população homogênea), threshold = 0
        thr = 0.7 * ent.max() if ent.max() > 0 else 0.0

        # 3. CRIAÇÃO DE BLOCOS ADAPTATIVOS
        blocks = []
        cur = 0  # Posição atual sendo processada

        while cur < self.L:
            # ESTRATÉGIA DE TAMANHO BASEADA EM ENTROPIA:
            # - Alta entropia (>threshold): bloco pequeno para refinamento
            # - Baixa entropia (≤threshold): bloco grande para preservação
            if ent[cur] > thr:
                # Região com alta variabilidade → bloco mínimo
                length = self.min_block_len
            else:
                # Região com baixa variabilidade → bloco duplo
                length = self.min_block_len * 2

            # Cria bloco respeitando limites da string
            end_pos = min(self.L, cur + length)
            blocks.append((cur, end_pos))

            # Avança para próxima região
            cur += length

        return blocks

    def _evaluate_population_parallel(self, pop: Population) -> list[tuple[str, int]]:
        """
        Avalia a população em paralelo usando threads para acelerar computação.

        A avaliação de fitness é frequentemente o gargalo computacional em
        algoritmos genéticos, especialmente para o CSP onde cada avaliação
        requer O(n*L) comparações. Esta função paraleliza o processo para
        aproveitar múltiplos cores do processador.

        ESTRATÉGIA DE PARALELIZAÇÃO:

        1. **DETECÇÃO AUTOMÁTICA**:
           - Se internal_workers ≤ 1: usa avaliação sequencial
           - Caso contrário: usa ThreadPoolExecutor para paralelização

        2. **DISTRIBUIÇÃO DE TRABALHO**:
           - Cada thread avalia um indivíduo independentemente
           - Função max_distance é thread-safe (apenas leitura)
           - Resultados coletados automaticamente pelo ThreadPoolExecutor

        3. **ORDENAÇÃO INTEGRADA**:
           - Retorna resultados já ordenados por fitness
           - Evita passo adicional de ordenação na função chamadora
           - Melhor cache locality para acessos subsequentes

        VANTAGENS DA PARALELIZAÇÃO:
        - **Speedup Linear**: Próximo ao número de cores disponíveis
        - **Overhead Baixo**: ThreadPoolExecutor reutiliza threads
        - **Escalabilidade**: Adapta-se automaticamente ao hardware
        - **Fallback Robusto**: Degrada graciosamente para modo sequencial

        CONSIDERAÇÕES DE PERFORMANCE:
        - Threads são eficazes pois avaliação é CPU-bound
        - Contém de contexto é mínimo (função pura)
        - Balanceamento automático de carga entre threads

        Args:
            pop: População de strings a serem avaliadas

        Returns:
            list[tuple[str, int]]: Lista de tuplas (string, fitness) ordenadas
                                  por fitness crescente (melhor primeiro)

        Note:
            O número de threads é controlado por self.internal_workers,
            permitindo ajuste baseado no hardware disponível.
        """
        # MODO SEQUENCIAL: Para casos simples ou debug
        if self.internal_workers <= 1:
            # Avaliação sequencial simples
            logger.debug(
                f"[PARALLEL-LOG] BLF-GA usando avaliação SEQUENCIAL (internal_workers={self.internal_workers})"
            )
            return [(s, max_distance(s, self.strings)) for s in pop]

        # Log temporário: paralelismo interno
        import threading

        thread_id = threading.get_ident()
        logger.debug(
            f"[PARALLEL-LOG] BLF-GA usando avaliação PARALELA com {self.internal_workers} workers internos (Thread principal: {thread_id})"
        )

        # FUNÇÃO DE AVALIAÇÃO INDIVIDUAL
        # Encapsula lógica para uso com ThreadPoolExecutor
        def evaluate_string(s: str) -> tuple[str, int]:
            """Avalia fitness de uma única string."""
            return (s, max_distance(s, self.strings))

        # PARALELIZAÇÃO COM THREADPOOLEXECUTOR
        # Usa context manager para limpeza automática de recursos
        with ThreadPoolExecutor(max_workers=self.internal_workers) as executor:
            # map() distribui evaluate_string para cada string em pop
            # Cada thread processa uma string independentemente
            # list() coleta resultados mantendo ordem original
            results = list(executor.map(evaluate_string, pop))

        logger.debug(
            f"[PARALLEL-LOG] BLF-GA avaliação paralela CONCLUÍDA - {len(pop)} indivíduos processados"
        )

        # ORDENAÇÃO POR FITNESS
        # Ordena do melhor (menor distância) para o pior (maior distância)
        # key=lambda x: x[1] usa o fitness (segundo elemento da tupla)
        return sorted(results, key=lambda x: x[1])

    def _sort_population_parallel(self, pop: Population) -> Population:
        """
        Ordena a população por fitness usando paralelização eficiente.

        Esta função é um wrapper conveniente que combina avaliação paralela
        com extração da população ordenada. Simplifica o código nas funções
        chamadoras e centraliza a lógica de ordenação populacional.

        FUNCIONAMENTO:
        1. Chama _evaluate_population_parallel para avaliação paralela
        2. Extrai apenas as strings da lista ordenada (descarta fitness)
        3. Retorna população ordenada pronta para uso

        VANTAGENS DO WRAPPER:
        - **Simplicidade**: Interface limpa para funções chamadoras
        - **Reutilização**: Aproveita avaliação paralela existente
        - **Consistência**: Garante ordenação uniforme em todo o código
        - **Manutenibilidade**: Centraliza lógica de ordenação

        Args:
            pop: População a ser ordenada

        Returns:
            Population: População ordenada por fitness (melhor → pior)

        Note:
            Esta função é equivalente a sorted(pop, key=fitness) mas
            aproveitando paralelização e cache de avaliações.
        """
        # 1. AVALIAÇÃO PARALELA
        # Obtém lista de tuplas (string, fitness) ordenadas
        evaluated = self._evaluate_population_parallel(pop)

        # 2. EXTRAÇÃO DA POPULAÇÃO
        # Extrai apenas as strings da lista ordenada
        # [s for s, _ in evaluated] usa list comprehension para eficiência
        return [s for s, _ in evaluated]
