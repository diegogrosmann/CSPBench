"""
H³-CSP: Hybrid Hierarchical Hamming Search for CSP.

This module implements the wrapper class that integrates the H³-CSP algorithm
into the CSP-BLFGA framework. H³-CSP is a three-layer hybrid algorithm
that combines hierarchical decomposition, specialized block techniques
and global refinement to solve the Closest String Problem.

H³-CSP Architecture:
1. B-Splitter: Divides strings into contiguous blocks (~√L)
2. Smart-Core: Selects optimal technique per block based on difficulty
3. Global Refine: Combines blocks and applies global hill-climbing

Classes:
    H3CSPAlgorithm: Wrapper for H³-CSP integration into CSP framework.
"""

from collections.abc import Callable

from src.domain.algorithms import CSPAlgorithm, register_algorithm

from .config import H3_CSP_DEFAULTS
from .implementation import H3CSP


@register_algorithm
class H3CSPAlgorithm(CSPAlgorithm):
    """
    H³-CSP: Hybrid Hierarchical Hamming Search for the Closest String Problem.

    H³-CSP is a hybrid algorithm that uses a three-layer approach
    to solve CSP efficiently. It divides the problem into smaller
    blocks, applies specialized techniques to each block, and then combines
    the results with global refinement.

    Features:
    - Hierarchical decomposition based on √L rule
    - Adaptive technique selection per block
    - Global refinement through local search (hill-climbing)
    - Deterministic and no internal parallelism support

    Main parameters:
    - beam_width: Beam search width (default: 32)
    - k_candidates: Number of candidates per block (default: 5)
    - local_iters: Local refinement iterations (default: 3)
    - max_time: Maximum execution time in seconds (default: 300)
    - seed: Seed for reproducibility (default: None)

    Args:
        strings (list[str]): List of input strings for CSP.
        alphabet (str): Alphabet used in strings.
        **params: Algorithm-specific parameters (see config.py).

    Attributes:
        name (str): Algorithm name ("H³-CSP").
        default_params (dict): Algorithm default parameters.
        supports_internal_parallel (bool): False - does not support internal parallelism.
        is_deterministic (bool): True - deterministic algorithm.
        h3_csp_instance (H3CSP): Instance of algorithm implementation.

    Example:
        >>> strings = ["ACGT", "AGCT", "ATCT"]
        >>> alg = H3CSPAlgorithm(strings, "ACGT", beam_width=16)
        >>> center, distance, metadata = alg.run()
        >>> print(f"Center: {center}, Distance: {distance}")
    """

    name = "H³-CSP"
    default_params = H3_CSP_DEFAULTS
    supports_internal_parallel = False  # H³-CSP does not support internal parallelism
    is_deterministic = True  # H³-CSP is deterministic

    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Initialize the H3CSPAlgorithm wrapper.

        Configure parameters and create an instance of the H3CSP implementation
        that will be used to execute the algorithm.

        Args:
            strings (list[str]): List of input strings for CSP.
                                All strings must have the same length.
            alphabet (str): Alphabet used in strings (e.g. "ACGT" for DNA).
            **params: Algorithm-specific parameters that will override
                     the default values defined in H3_CSP_DEFAULTS.

        Raises:
            ValueError: If strings have different lengths.
            ValueError: If alphabet is empty.
        """
        super().__init__(strings, alphabet, **params)
        self.h3_csp_instance = H3CSP(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define a callback to monitor execution progress.

        The callback will be called during algorithm execution with
        descriptive messages about current progress, enabling
        graphical interfaces or progress logs.
    """

    name = "H³-CSP"
    default_params = H3_CSP_DEFAULTS
    supports_internal_parallel = False  # H³-CSP does not support internal parallelism
    is_deterministic = True  # H³-CSP is deterministic

    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Initialize the H3CSPAlgorithm wrapper.

        Configure parameters and create an instance of the H3CSP
        implementation that will be used to execute the algorithm.

        Args:
            strings (list[str]): List of input strings for CSP.
                                All strings must have the same length.
            alphabet (str): Alphabet used in strings (e.g., "ACGT" for DNA).
            **params: Algorithm-specific parameters that will override
                     the default values defined in H3_CSP_DEFAULTS.

        Raises:
            ValueError: If strings have different lengths.
            ValueError: If alphabet is empty.
        """
        super().__init__(strings, alphabet, **params)
        self.h3_csp_instance = H3CSP(self.strings, self.alphabet, **self.params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define a callback to monitor execution progress.

        The callback will be called during algorithm execution with
        descriptive messages about current progress, allowing
        graphical interfaces or progress logs.

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
