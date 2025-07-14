"""
Domínio: Algoritmos CSP

Este módulo contém a lógica central de algoritmos para o Closest String Problem (CSP),
incluindo interfaces abstratas, registry de algoritmos e operadores genéticos.
Livre de dependências externas conforme arquitetura hexagonal.
"""

import random
from abc import ABC, abstractmethod
from typing import Any, Callable, Optional

# =============================================================================
# REGISTRY DE ALGORITMOS
# =============================================================================

global_registry: dict[str, type] = {}


def register_algorithm(cls: type) -> type:
    """
    Decorador para registrar uma classe de algoritmo no registry global.

    Args:
        cls: Classe do algoritmo a ser registrada

    Returns:
        type: A própria classe, permitindo uso como decorador
    """
    algorithm_name = getattr(cls, "name", cls.__name__)
    global_registry[algorithm_name] = cls
    return cls


# =============================================================================
# INTERFACES DE ALGORITMOS
# =============================================================================


class CSPAlgorithm(ABC):
    """Interface base abstrata para todos os algoritmos CSP."""

    # Atributos de classe obrigatórios
    name: str
    default_params: dict
    is_deterministic: bool = False
    supports_internal_parallel: bool = False

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Inicializa o algoritmo com as strings e alfabeto.

        Args:
            strings: Lista de strings do dataset
            alphabet: Alfabeto utilizado
            **params: Parâmetros específicos do algoritmo
        """
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}
        self.progress_callback: Optional[Callable[[str], None]] = None
        self.warning_callback: Optional[Callable[[str], None]] = None

        # Configurações de histórico
        self.save_history = params.get("save_history", False)
        self.history_frequency = params.get(
            "history_frequency", 1
        )  # A cada N iterações
        self.history = []

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """Define callback para relatar progresso do algoritmo."""
        self.progress_callback = callback

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """Define callback para relatar warnings do algoritmo."""
        self.warning_callback = callback

    def _report_progress(self, message: str) -> None:
        """Relata progresso se callback estiver definido."""
        if self.progress_callback:
            self.progress_callback(message)

    def _report_warning(self, message: str) -> None:
        """Relata warning se callback estiver definido."""
        if self.warning_callback:
            self.warning_callback(message)

    def _save_history_entry(self, iteration: int, **data) -> None:
        """
        Salva uma entrada no histórico se habilitado.

        Args:
            iteration: Número da iteração atual
            **data: Dados do estado atual (fitness, melhor solução, etc.)
        """
        if self.save_history and (iteration % self.history_frequency == 0):
            entry = {"iteration": iteration, "timestamp": self._get_timestamp(), **data}
            self.history.append(entry)

    def _get_timestamp(self) -> float:
        """Retorna timestamp atual para histórico."""
        import time

        return time.time()

    def get_history(self) -> list[dict]:
        """Retorna o histórico de execução."""
        return self.history.copy()

    def clear_history(self) -> None:
        """Limpa o histórico atual."""
        self.history.clear()

    @abstractmethod
    def run(self) -> tuple[str, int, dict[str, Any]]:
        """
        Executa o algoritmo e retorna resultado estruturado.

        Returns:
            tuple: (string_central, distancia_maxima, metadata)
                metadata deve conter:
                - history: Lista de estados durante execução (se habilitado)
                - iterations: Número total de iterações
                - convergence_data: Dados de convergência (se aplicável)
                - outras informações específicas do algoritmo
        """

    def set_params(self, **params) -> None:
        """Define novos parâmetros para o algoritmo."""
        self.params.update(params)

    def get_metadata(self) -> dict[str, Any]:
        """Retorna metadados do algoritmo."""
        return {
            "name": self.name,
            "params": self.params.copy(),
            "is_deterministic": self.is_deterministic,
            "supports_internal_parallel": self.supports_internal_parallel,
            "input_size": len(self.strings),
            "string_length": len(self.strings[0]) if self.strings else 0,
            "alphabet_size": len(self.alphabet),
        }


class Algorithm(ABC):
    """Interface legacy para compatibilidade com código existente."""

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str):
        """Inicializa algoritmo legacy."""
        self.strings = strings
        self.alphabet = alphabet

    @abstractmethod
    def solve(self) -> str:
        """Resolve o problema CSP retornando string central."""


# =============================================================================
# OPERADORES GENÉTICOS
# =============================================================================

String = str
Population = list[String]


def mean_hamming_distance(pop: Population) -> float:
    """
    Calcula a distância de Hamming média entre todos os pares de indivíduos.

    Args:
        pop: População de strings

    Returns:
        float: Distância média entre pares
    """
    if len(pop) < 2:
        return 0.0

    total_distance = 0
    total_pairs = 0

    for i in range(len(pop)):
        for j in range(i + 1, len(pop)):
            distance = sum(c1 != c2 for c1, c2 in zip(pop[i], pop[j]))
            total_distance += distance
            total_pairs += 1

    return total_distance / total_pairs if total_pairs > 0 else 0.0


def mutate_multi(ind: str, alphabet: str, rng: random.Random, n: int = 2) -> str:
    """
    Realiza mutação multi-ponto alterando até n posições aleatórias.

    Args:
        ind: String original
        alphabet: Conjunto de símbolos válidos
        rng: Gerador de números aleatórios
        n: Número de mutações a aplicar

    Returns:
        str: String mutada
    """
    chars = list(ind)
    L = len(chars)

    for _ in range(n):
        pos = rng.randint(0, L - 1)
        old = chars[pos]

        # Escolhe símbolo diferente do atual
        available = [c for c in alphabet if c != old]
        if available:
            chars[pos] = rng.choice(available)

    return "".join(chars)


def mutate_inversion(ind: str, rng: random.Random) -> str:
    """
    Realiza mutação por inversão de um segmento aleatório.

    Args:
        ind: String original
        rng: Gerador de números aleatórios

    Returns:
        str: String com segmento invertido
    """
    chars = list(ind)
    L = len(chars)

    if L < 2:
        return ind

    # Seleciona dois pontos aleatórios
    i, j = sorted(rng.sample(range(L), 2))

    # Inverte o segmento entre i e j
    chars[i : j + 1] = chars[i : j + 1][::-1]

    return "".join(chars)


def crossover_one_point(
    parent1: str, parent2: str, rng: random.Random
) -> tuple[str, str]:
    """
    Realiza crossover de um ponto entre dois pais.

    Args:
        parent1: Primeiro pai
        parent2: Segundo pai
        rng: Gerador de números aleatórios

    Returns:
        tuple: Dois filhos resultantes
    """
    L = len(parent1)
    if L <= 1:
        return parent1, parent2

    # Seleciona ponto de corte
    cut_point = rng.randint(1, L - 1)

    # Cria filhos trocando sufixos
    child1 = parent1[:cut_point] + parent2[cut_point:]
    child2 = parent2[:cut_point] + parent1[cut_point:]

    return child1, child2


def crossover_uniform(
    parent1: str, parent2: str, rng: random.Random, rate: float = 0.5
) -> tuple[str, str]:
    """
    Realiza crossover uniforme entre dois pais.

    Args:
        parent1: Primeiro pai
        parent2: Segundo pai
        rng: Gerador de números aleatórios
        rate: Taxa de troca por posição

    Returns:
        tuple: Dois filhos resultantes
    """
    chars1 = list(parent1)
    chars2 = list(parent2)

    for i in range(len(chars1)):
        if rng.random() < rate:
            chars1[i], chars2[i] = chars2[i], chars1[i]

    return "".join(chars1), "".join(chars2)


def refine_greedy(individual: str, strings: list[str]) -> str:
    """
    Refinamento local guloso posição por posição.

    Args:
        individual: String a ser refinada
        strings: Strings de referência

    Returns:
        str: String refinada
    """
    chars = list(individual)
    alphabet = set("".join(strings))

    for pos in range(len(chars)):
        current_char = chars[pos]
        best_char = current_char
        best_distance = _max_distance_to_strings(chars, strings)

        # Testa cada símbolo do alfabeto
        for symbol in alphabet:
            if symbol != current_char:
                chars[pos] = symbol
                distance = _max_distance_to_strings(chars, strings)

                if distance < best_distance:
                    best_distance = distance
                    best_char = symbol

        chars[pos] = best_char

    return "".join(chars)


def _max_distance_to_strings(chars: list[str], strings: list[str]) -> int:
    """
    Calcula distância máxima de um candidato para conjunto de strings.

    Args:
        chars: Lista de caracteres formando string candidata
        strings: Strings de referência

    Returns:
        int: Distância máxima
    """
    candidate = "".join(chars)
    max_dist = 0

    for string in strings:
        dist = sum(c1 != c2 for c1, c2 in zip(candidate, string))
        max_dist = max(max_dist, dist)

    return max_dist


# =============================================================================
# ALGORITMOS ESPECÍFICOS
# =============================================================================

# Algoritmos específicos foram movidos para a pasta algorithms/
# para serem carregados dinamicamente como plug-ins
