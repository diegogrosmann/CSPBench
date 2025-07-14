"""
H³-CSP: Hybrid Hierarchical Hamming Search para CSP.

Este módulo implementa a classe wrapper que integra o algoritmo H³-CSP
ao framework CSP-BLFGA. O H³-CSP é um algoritmo híbrido de três camadas
que combina decomposição hierárquica, técnicas especializadas por bloco
e refinamento global para resolver o Closest String Problem.

Arquitetura do H³-CSP:
1. B-Splitter: Divide as strings em blocos contíguos (~√L)
2. Smart-Core: Seleciona técnica ótima por bloco baseada na dificuldade
3. Global Refine: Combina blocos e aplica hill-climbing global

Classes:
    H3CSPAlgorithm: Wrapper para integração do H³-CSP ao framework CSP.
"""

from collections.abc import Callable

from src.domain.algorithms import CSPAlgorithm, register_algorithm

from .config import H3_CSP_DEFAULTS
from .implementation import H3CSP


@register_algorithm
class H3CSPAlgorithm(CSPAlgorithm):
    """
    H³-CSP: Hybrid Hierarchical Hamming Search para o Closest String Problem.

    O H³-CSP é um algoritmo híbrido que utiliza uma abordagem de três camadas
    para resolver o CSP de forma eficiente. Ele divide o problema em blocos
    menores, aplica técnicas especializadas em cada bloco e depois combina
    os resultados com refinamento global.

    Características:
    - Decomposição hierárquica baseada na regra √L
    - Seleção adaptativa de técnicas por bloco
    - Refinamento global por busca local (hill-climbing)
    - Determinístico e sem suporte a paralelismo interno

    Parâmetros principais:
    - beam_width: Largura do beam search (padrão: 32)
    - k_candidates: Número de candidatos por bloco (padrão: 5)
    - local_iters: Iterações de refinamento local (padrão: 3)
    - max_time: Tempo máximo de execução em segundos (padrão: 300)
    - seed: Semente para reprodutibilidade (padrão: None)

    Args:
        strings (list[str]): Lista de strings de entrada para o CSP.
        alphabet (str): Alfabeto utilizado nas strings.
        **params: Parâmetros específicos do algoritmo (ver config.py).

    Attributes:
        name (str): Nome do algoritmo ("H³-CSP").
        default_params (dict): Parâmetros padrão do algoritmo.
        supports_internal_parallel (bool): False - não suporta paralelismo interno.
        is_deterministic (bool): True - algoritmo determinístico.
        h3_csp_instance (H3CSP): Instância da implementação do algoritmo.

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> alg = H3CSPAlgorithm(strings, "ACGT", beam_width=16)
        >>> center, distance, metadata = alg.run()
        >>> print(f"Centro: {center}, Distância: {distance}")
    """

    name = "H³-CSP"
    default_params = H3_CSP_DEFAULTS
    supports_internal_parallel = False  # H³-CSP não suporta paralelismo interno
    is_deterministic = True  # H³-CSP é determinístico

    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Inicializa o wrapper H3CSPAlgorithm.

        Configura os parâmetros e cria uma instância da implementação
        H3CSP que será utilizada para executar o algoritmo.

        Args:
            strings (list[str]): Lista de strings de entrada para o CSP.
                                Todas as strings devem ter o mesmo comprimento.
            alphabet (str): Alfabeto utilizado nas strings (ex: "ACGT" para DNA).
            **params: Parâmetros específicos do algoritmo que sobrescreverão
                     os valores padrão definidos em H3_CSP_DEFAULTS.

        Raises:
            ValueError: Se as strings tiverem comprimentos diferentes.
            ValueError: Se o alfabeto estiver vazio.
        """
        super().__init__(strings, alphabet, **params)
        self.h3_csp_instance = H3CSP(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define um callback para monitorar o progresso da execução.

        O callback será chamado durante a execução do algoritmo com
        mensagens descritivas sobre o progresso atual, permitindo
        interfaces gráficas ou logs de progresso.

        Args:
            callback (Callable[[str], None]): Função que será chamada
                                             com mensagens de progresso.
                                             Deve aceitar uma string como
                                             parâmetro e não retornar nada.

        Example:
            >>> def progress_handler(msg):
            ...     print(f"Progresso: {msg}")
            >>> alg.set_progress_callback(progress_handler)
        """
        super().set_progress_callback(callback)
        self.h3_csp_instance.set_progress_callback(callback)

    def run(self) -> tuple[str, int, dict]:
        """
        Executa o algoritmo H³-CSP e retorna o resultado.

        O algoritmo seguirá as seguintes etapas:
        1. Divisão das strings em blocos (~√L blocos)
        2. Análise de cada bloco para determinar a técnica ótima
        3. Geração de candidatos por bloco usando a técnica selecionada
        4. Fusão dos melhores candidatos de cada bloco
        5. Refinamento global usando hill-climbing

        Returns:
            tuple[str, int, dict]: Tupla contendo:
                - str: String central encontrada (solução do CSP)
                - int: Distância máxima da string central para as strings originais
                - dict: Metadados da execução contendo:
                    - iteracoes: Número de iterações realizadas
                    - algoritmo: Nome do algoritmo ("H³-CSP")
                    - parametros_usados: Parâmetros utilizados na execução
                    - centro_encontrado: String central encontrada

        Raises:
            TimeoutError: Se o tempo máximo de execução for excedido.
            RuntimeError: Se ocorrer um erro durante a execução do algoritmo.

        Example:
            >>> alg = H3CSPAlgorithm(["ACGT", "AGCT", "ATCT"], "ACGT")
            >>> center, distance, metadata = alg.run()
            >>> print(f"Solução: {center} com distância {distance}")
        """
        # Salvar estado inicial no histórico se habilitado
        if self.save_history:
            self._save_history_entry(
                0,
                phase="initialization",
                parameters=self.params,
                message="Iniciando algoritmo H³-CSP",
            )

        center, dist = self.h3_csp_instance.run()

        metadata = {
            "iteracoes": getattr(self.h3_csp_instance, "iterations", 1),
            "algoritmo": "H³-CSP",
            "parametros_usados": self.params,
            "centro_encontrado": center,
        }

        # Salvar estado final no histórico se habilitado
        if self.save_history:
            self._save_history_entry(
                metadata["iteracoes"],
                phase="completion",
                best_fitness=dist,
                best_solution=center,
                message="Algoritmo H³-CSP finalizado",
            )

            # Adicionar histórico aos metadados
            metadata["history"] = self.get_history()

        return center, dist, metadata
