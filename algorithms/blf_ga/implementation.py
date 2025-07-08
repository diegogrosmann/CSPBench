"""
Implementação da heurística BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) para CSP.

O BLF-GA é uma metaheurística híbrida que combina três estratégias principais:
1. **Blockwise Learning**: Divide as strings em blocos e aprende padrões locais
2. **Genetic Algorithm**: Usa operações genéticas para evolução global
3. **Fusion**: Combina conhecimento local e global para melhor convergência

ALGORITMO GERAL:
├── Inicialização
│   ├── Criar população inicial (consenso + variações + aleatórias)
│   └── Dividir strings em blocos adaptativos
├── Loop Principal (até critério de parada)
│   ├── Aprendizado por Blocos
│   │   ├── Para cada bloco, encontrar padrão consenso
│   │   └── Criar repositório de conhecimento local
│   ├── Evolução Genética
│   │   ├── Seleção por torneio
│   │   ├── Crossover (one-point, uniform, blend-blocks)
│   │   ├── Mutação (multi, inversion, transposition)
│   │   └── Elitismo adaptativo
│   ├── Mecanismos Adaptativos
│   │   ├── Redivisão de blocos baseada em entropia
│   │   ├── Imigrantes aleatórios para diversidade
│   │   ├── Mutação adaptativa baseada em convergência
│   │   └── Refinamento local dos melhores indivíduos
│   └── Critérios de Parada
│       ├── Solução ótima encontrada
│       ├── Limite de tempo/gerações
│       └── Early stopping por estagnação
└── Retorno da melhor solução encontrada

Classes:
    BLFGA: Implementa o algoritmo BLF-GA.

Funções auxiliares:
    hamming_dist(a, b): Wrapper para distância de Hamming.
    max_distance(center, strings): Wrapper para distância máxima.
"""

import logging
import random
import time
from collections import Counter
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from src.utils.distance import hamming_distance, max_distance

from .config import BLF_GA_DEFAULTS
from .ops import genetic_ops

logger = logging.getLogger(__name__)

String = str
Population = list[String]


def hamming_dist(a: String, b: String) -> int:
    """Wrapper para manter compatibilidade."""
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
        # Temporariamente desabilitado para debug
        self.internal_workers = 1  # int(os.environ.get('INTERNAL_WORKERS', '1'))

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
        Para pop_size (mode='multiplier'):
            - Se param for float >= 1, retorna int(param * ref_value)
            - Se param for int, retorna param
        Para initial_blocks (mode='proportion'):
            - Se param for float entre 0 e 1, retorna int(param * ref_value)
            - Se param for int, retorna param

        Args:
            param: Parâmetro a ser resolvido
            ref_value: Valor de referência para cálculo
            mode: Modo de resolução ('multiplier' ou 'proportion')
            min_value: Valor mínimo a ser garantido (opcional)
        """
        if mode == "multiplier":
            if isinstance(param, float) and param >= 1:
                result = int(param * ref_value)
            else:
                result = int(param)
        else:
            if isinstance(param, float) and 0 < param <= 1:
                result = int(param * ref_value)
            else:
                result = int(param)

        # Garante valor mínimo se especificado
        if min_value is not None:
            result = max(min_value, result)
        else:
            result = max(1, result)

        return result

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Define um callback para relatar o progresso."""
        self.progress_callback = callback

    def update_params(self, **params):
        """
        Atualiza parâmetros do algoritmo dinamicamente.

        Args:
            **params: Parâmetros a serem atualizados
        """
        # Atualizar parâmetros básicos
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)

        # Reprocessar parâmetros dinâmicos se necessário
        if "pop_size" in params:
            min_pop_size = params.get("min_pop_size", self.min_pop_size)
            self.pop_size = self._resolve_dynamic_param(
                params["pop_size"], self.n, mode="multiplier", min_value=min_pop_size
            )

        if "initial_blocks" in params:
            self.initial_blocks = self._resolve_dynamic_param(
                params["initial_blocks"], self.L, mode="proportion"
            )
            # Reinicializar blocos se necessário
            self.blocks = self._initial_blocking()

        if "no_improve_patience" in params:
            patience_param = params["no_improve_patience"]
            if isinstance(patience_param, float) and 0 < patience_param < 1:
                self.no_improve_patience = max(1, int(patience_param * self.max_gens))
            else:
                self.no_improve_patience = int(patience_param)

        if "seed" in params:
            self.rng = random.Random(params["seed"])

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

            # --- MECANISMO 2: MUTAÇÃO ADAPTATIVA ---
            # Ajusta taxa de mutação baseada na diversidade populacional e convergência
            diversity = genetic_ops.mean_hamming_distance(pop)

            # Se diversidade baixa, aumenta mutação temporariamente
            if diversity < self.diversity_threshold * self.L:
                self.mut_prob = mut_prob_backup * self.mutation_adapt_factor
                mut_adapt_timer = self.mutation_adapt_duration

            # Se fitness estagnado por N gerações, aumenta mutação
            if len(self.history) > self.mutation_adapt_N and all(
                self.history[-i] == self.history[-1]
                for i in range(1, self.mutation_adapt_N + 1)
            ):
                self.mut_prob = mut_prob_backup * self.mutation_adapt_factor
                mut_adapt_timer = self.mutation_adapt_duration

            # Decrementa timer da mutação adaptativa
            if mut_adapt_timer > 0:
                mut_adapt_timer -= 1
                if mut_adapt_timer == 0:
                    self.mut_prob = mut_prob_backup

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
            else:
                elites = pop[:k]

            # --- MECANISMO 6: REFINAMENTO LOCAL ---
            # Aplica busca local intensiva nos melhores indivíduos
            if self.refine_elites == "all":
                pop[:k] = self._refine_elites(elites)  # Refina todos os elites
            else:
                if elites:
                    pop[0] = self._refine_elites([elites[0]])[
                        0
                    ]  # Refina apenas o melhor

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
                for i in range(n_restart):
                    pop[-(i + 1)] = "".join(
                        self.rng.choice(self.alphabet) for _ in range(self.L)
                    )
                no_improve = 0

            # --- CRITÉRIO DE PARADA 1: EARLY STOPPING ---
            # Encerra se não houver melhoria por X gerações
            if self.no_improve_patience and no_improve >= self.no_improve_patience:
                if self.progress_callback:
                    self.progress_callback(
                        f"Encerrando por early stopping após {no_improve} gerações sem melhoria."
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
                self.blocks = self._adaptive_blocking(pop)

        # === FASE 3: FINALIZAÇÃO ===
        return best, best_val, self.history

    def _init_population(self) -> Population:
        """
        Cria população inicial usando estratégia híbrida inteligente.

        ESTRATÉGIA DE INICIALIZAÇÃO:
        1. **Consenso Global**: Cria uma string consenso baseada na moda de cada posição
        2. **Variações Inteligentes**: Cria 1/3 da população substituindo blocos do consenso
           por blocos correspondentes das strings originais
        3. **Diversidade Aleatória**: Preenche o restante com strings completamente aleatórias

        Esta estratégia garante:
        - Qualidade inicial alta (consenso)
        - Diversidade estrutural (variações por blocos)
        - Exploração ampla (aleatoriedade)

        Returns:
            Population: Lista de strings representando a população inicial
        """
        # 1. Cria consenso: para cada posição, pega o símbolo mais frequente
        consensus = "".join(
            Counter(pos).most_common(1)[0][0] for pos in zip(*self.strings)
        )
        pop = [consensus]

        # 2. Cria variações inteligentes do consenso
        for _ in range(self.pop_size // 3):
            s = list(consensus)
            # Para cada bloco, com 50% de chance substitui por bloco de string original
            for l, r in self.blocks:
                if self.rng.random() < 0.5:
                    src = self.rng.choice(self.strings)
                    s[l:r] = src[l:r]
            pop.append("".join(s))

        # 3. Preenche restante com strings aleatórias para diversidade
        while len(pop) < self.pop_size:
            rand_s = "".join(self.rng.choice(self.alphabet) for _ in range(self.L))
            pop.append(rand_s)
        return pop

    def _initial_blocking(self) -> list[tuple[int, int]]:
        """
        Cria divisão inicial em blocos de tamanho uniforme.

        Returns:
            list[tuple[int, int]]: Lista de tuplas (início, fim) definindo os blocos
        """
        size = max(self.min_block_len, self.L // self.initial_blocks)
        blocks = [(i, min(i + size, self.L)) for i in range(0, self.L, size)]
        return blocks

    def _learn_blocks(self, pop: Population) -> list[String]:
        """
        Extrai conhecimento local de cada bloco através de aprendizado por consenso.

        FUNCIONAMENTO:
        1. Para cada bloco definido, coleta todos os segmentos correspondentes da população
        2. Calcula consenso posição por posição dentro do bloco
        3. Constrói repositório de padrões aprendidos

        Este é o coração do "Learning" no BLF-GA:
        - Identifica padrões locais eficazes
        - Cria conhecimento reutilizável
        - Guia a evolução com informação estruturada

        Args:
            pop: População atual

        Returns:
            list[String]: Repositório de padrões consenso por bloco
        """
        repo = []
        for l, r in self.blocks:
            # Coleta segmentos do bloco atual de todos os indivíduos
            symbols = [ind[l:r] for ind in pop if len(ind) >= r]
            # Cria consenso posição por posição dentro do bloco
            block = "".join(
                Counter(chars).most_common(1)[0][0] for chars in zip(*symbols)
            )
            repo.append(block)
        return repo

    def _tournament_selection(self, pop: Population, k: int = 2) -> String:
        """
        Seleciona indivíduo através de torneio.

        FUNCIONAMENTO:
        1. Escolhe k indivíduos aleatoriamente da população
        2. Retorna o melhor dentre eles (menor distância máxima)

        Args:
            pop: População para seleção
            k: Tamanho do torneio

        Returns:
            String: Indivíduo selecionado
        """
        if k is None:
            k = self.tournament_k if self.tournament_k is not None else 2
        tournament_pool = self.rng.sample(pop, k)
        winner = min(tournament_pool, key=lambda s: max_distance(s, self.strings))
        return winner

    def _apply_crossover(self, p1, p2):
        """
        Aplica operação de crossover entre dois pais.

        TIPOS DE CROSSOVER:
        - one_point: Ponto único de corte
        - uniform: Troca uniforme gene a gene
        - blend_blocks: Mistura blocos inteiros (específico do BLF-GA)

        Args:
            p1, p2: Pais para crossover

        Returns:
            tuple: (filho1, filho2)
        """
        if self.crossover_type == "one_point":
            c1, c2 = genetic_ops.crossover_one_point(p1, p2, self.rng)
        elif self.crossover_type == "uniform":
            c1, c2 = genetic_ops.crossover_uniform(p1, p2, self.rng)
        elif self.crossover_type == "blend_blocks":
            c1, c2 = genetic_ops.crossover_blend_blocks(p1, p2, self.blocks, self.rng)
        else:
            c1, c2 = p1, p2
        return c1, c2

    def _apply_mutation(self, ind):
        """
        Aplica operação de mutação a um indivíduo.

        TIPOS DE MUTAÇÃO:
        - multi: Muda N posições aleatórias
        - inversion: Inverte um segmento
        - transposition: Transpõe dois segmentos

        Args:
            ind: Indivíduo para mutar

        Returns:
            String: Indivíduo mutado
        """
        if self.mutation_type == "multi":
            return genetic_ops.mutate_multi(
                ind, self.alphabet, self.rng, n=self.mutation_multi_n
            )
        elif self.mutation_type == "inversion":
            return genetic_ops.mutate_inversion(ind, self.rng)
        elif self.mutation_type == "transposition":
            return genetic_ops.mutate_transposition(ind, self.rng)
        else:
            return ind

    def _apply_refinement(self, ind):
        """
        Aplica refinamento local (busca local) a um indivíduo.

        TIPOS DE REFINAMENTO:
        - greedy: Busca gulosa posição por posição
        - swap: Troca de posições
        - insertion: Inserção de segmentos
        - 2opt: Otimização 2-opt

        Args:
            ind: Indivíduo para refinar

        Returns:
            String: Indivíduo refinado
        """
        if self.refinement_type == "greedy":
            return genetic_ops.refine_greedy(ind, self.strings)
        elif self.refinement_type == "swap":
            return genetic_ops.refine_swap(ind, self.strings)
        elif self.refinement_type == "insertion":
            return genetic_ops.refine_insertion(ind, self.strings)
        elif self.refinement_type == "2opt":
            return genetic_ops.refine_2opt(ind, self.strings)
        else:
            return ind

    def _next_generation(
        self, pop: Population, repo: list[String]
    ) -> Population:  # pylint: disable=unused-argument
        """
        Gera nova população através de operações genéticas clássicas.

        ALGORITMO GENÉTICO PADRÃO:
        1. **Avaliação**: Calcula fitness de todos os indivíduos
        2. **Elitismo**: Preserva os melhores indivíduos automaticamente
        3. **Reprodução**: Gera filhos até completar população através de:
           - Seleção por torneio
           - Crossover (com probabilidade cross_prob)
           - Mutação (com probabilidade mut_prob)

        Args:
            pop: População atual
            repo: Repositório de conhecimento (não usado nesta implementação)

        Returns:
            Population: Nova população gerada
        """
        # 1. Avaliar e ordenar população usando paralelismo
        evaluated_pop = self._evaluate_population_parallel(pop)
        pop_sorted = [s for s, _ in evaluated_pop]

        # 2. Elitismo - preserva melhores indivíduos
        elite_n = max(1, int(self.elite_rate * self.pop_size))
        new_pop = pop_sorted[:elite_n]

        # 3. Reprodução - gera filhos até completar população
        while len(new_pop) < self.pop_size:
            p1 = self._tournament_selection(pop_sorted, k=self.tournament_k)
            p2 = self._tournament_selection(pop_sorted, k=self.tournament_k)
            if self.rng.random() < self.cross_prob:
                c1, c2 = self._apply_crossover(p1, p2)
            else:
                c1, c2 = p1, p2
            c1 = self._apply_mutation(c1)
            c2 = self._apply_mutation(c2)
            new_pop.append(c1)
            if len(new_pop) < self.pop_size:
                new_pop.append(c2)

        # 4. Aplica niching se habilitado
        if self.niching:
            new_pop = self._apply_niching(new_pop)

        return new_pop[: self.pop_size]

    def _refine_elites(self, pop: Population) -> Population:
        # Paralelização opcional pode ser feita aqui
        def refine(ind):
            return self._apply_refinement(ind)

        with ThreadPoolExecutor() as executor:
            refined = list(executor.map(refine, pop))
        return refined

    def _apply_niching(self, pop: Population) -> Population:
        """
        Aplica niching para preservar diversidade local na população.

        O niching mantém distância mínima entre indivíduos para evitar
        convergência prematura e preservar múltiplas regiões promissoras.

        Args:
            pop: População para aplicar niching

        Returns:
            Population: População com niching aplicado
        """
        if not self.niching or len(pop) <= 1:
            return pop

        # Ordena população por fitness
        pop_with_fitness = [(ind, max_distance(ind, self.strings)) for ind in pop]
        pop_with_fitness.sort(key=lambda x: x[1])

        # Aplica niching mantendo diversidade
        niched_pop = []
        for ind, fitness in pop_with_fitness:
            # Verifica se está muito próximo de algum indivíduo já aceito
            too_close = False
            for accepted_ind in niched_pop:
                if hamming_dist(ind, accepted_ind) < self.niching_radius:
                    too_close = True
                    break

            # Aceita se não está muito próximo ou se a população ainda é pequena
            if not too_close or len(niched_pop) < self.pop_size // 2:
                niched_pop.append(ind)

            # Para quando tiver indivíduos suficientes
            if len(niched_pop) >= self.pop_size:
                break

        # Preenche com indivíduos aleatórios se necessário
        while len(niched_pop) < self.pop_size:
            rand_s = "".join(self.rng.choice(self.alphabet) for _ in range(self.L))
            niched_pop.append(rand_s)

        return niched_pop

    def _adaptive_blocking(self, pop: Population) -> list[tuple[int, int]]:
        """
        Redivide strings em blocos adaptativos baseado na entropia de cada posição.

        FUNCIONAMENTO:
        1. **Cálculo de Entropia**: Para cada posição, calcula entropia baseada na
           distribuição de símbolos na população atual
        2. **Threshold Adaptativo**: Define limiar baseado na entropia máxima
        3. **Divisão Inteligente**:
           - Posições com alta entropia (mais diversas) → blocos menores
           - Posições com baixa entropia (mais consenso) → blocos maiores

        INTUIÇÃO:
        - Alta entropia = muita variação = precisa de blocos menores para capturar detalhes
        - Baixa entropia = consenso claro = pode usar blocos maiores

        Args:
            pop: População atual para análise

        Returns:
            list[tuple[int, int]]: Nova divisão em blocos
        """
        # 1. Calcula entropia para cada posição
        ent = np.zeros(self.L)
        for pos in range(self.L):
            # Considere apenas strings que são longas o suficiente
            valid_chars = [ind[pos] for ind in pop if pos < len(ind)]
            if valid_chars:  # Se temos pelo menos um caractere válido
                cnt = Counter(valid_chars)
                probs = np.array(list(cnt.values())) / len(valid_chars)
                ent[pos] = -np.sum(probs * np.log2(probs))

        # 2. Define threshold baseado na entropia máxima
        thr = 0.7 * ent.max() if ent.max() > 0 else 0.0

        # 3. Cria blocos adaptativos
        blocks = []
        cur = 0
        while cur < self.L:
            # Alta entropia → bloco menor, baixa entropia → bloco maior
            length = self.min_block_len if ent[cur] > thr else self.min_block_len * 2
            blocks.append((cur, min(self.L, cur + length)))
            cur += length
        return blocks

    def _evaluate_population_parallel(self, pop: Population) -> list[tuple[str, int]]:
        """
        Avalia a população em paralelo usando threads.

        Args:
            pop: População a ser avaliada

        Returns:
            Lista de tuplas (string, fitness) ordenadas por fitness
        """
        if self.internal_workers <= 1:
            # Se apenas 1 worker, usar avaliação sequencial
            return [(s, max_distance(s, self.strings)) for s in pop]

        # Função para avaliar uma string
        def evaluate_string(s: str) -> tuple[str, int]:
            return (s, max_distance(s, self.strings))

        # Usar ThreadPoolExecutor para paralelizar avaliações
        with ThreadPoolExecutor(max_workers=self.internal_workers) as executor:
            results = list(executor.map(evaluate_string, pop))

        # Ordenar por fitness (distância máxima)
        return sorted(results, key=lambda x: x[1])

    def _sort_population_parallel(self, pop: Population) -> Population:
        """
        Ordena a população por fitness usando paralelismo.

        Args:
            pop: População a ser ordenada

        Returns:
            População ordenada por fitness
        """
        evaluated = self._evaluate_population_parallel(pop)
        return [s for s, _ in evaluated]
