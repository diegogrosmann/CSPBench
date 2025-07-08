# Guia do Desenvolvedor: Adicionando Novos Algoritmos

## 🎯 Visão Geral

Este guia explica como adicionar novos algoritmos ao CSP-BLFGA de forma padronizada e extensível. O sistema foi projetado para facilitar a adição de novos algoritmos sem modificar o código central.

## 🏗️ Arquitetura de Algoritmos

### Sistema de Registro Automático

O CSP-BLFGA usa um sistema de registro automático baseado em decoradores:

```python
from algorithms.base import CSPAlgorithm, register_algorithm

@register_algorithm  # ← Este decorador registra automaticamente
class MeuAlgoritmo(CSPAlgorithm):
    name = "MeuAlgoritmo"  # Nome que aparece na interface
    # ... implementação
```

### Interface Padronizada

Todos os algoritmos devem implementar a interface `CSPAlgorithm`:

```python
class IAlgorithm(Protocol):
    def run(self) -> Tuple[str, int, Dict[str, Any]]: ...
    def set_progress_callback(self, callback: Callable[[str], None]) -> None: ...
    def set_warning_callback(self, callback: Callable[[str], None]) -> None: ...
```

## 📁 Estrutura de Diretórios

### Organização Recomendada

```
algorithms/
├── base.py                    # Classes base e registry
├── meu_algoritmo/            # Seu novo algoritmo
│   ├── __init__.py
│   ├── algorithm.py          # Classe principal
│   ├── implementation.py     # Lógica do algoritmo
│   ├── config.py            # Configurações
│   ├── README.md            # Documentação
│   ├── TECHNICAL_DOCUMENTATION.md
│   └── tests/
│       ├── __init__.py
│       ├── test_algorithm.py
│       └── test_implementation.py
└── ...
```

### Exemplo de `__init__.py`

```python
# algorithms/meu_algoritmo/__init__.py
"""
Meu Algoritmo para CSP.

Este módulo implementa um algoritmo personalizado para o Closest String Problem.
"""

from .algorithm import MeuAlgoritmo

__all__ = ["MeuAlgoritmo"]
```

## 🧬 Implementação Passo a Passo

### 1. Classe Principal (`algorithm.py`)

```python
# algorithms/meu_algoritmo/algorithm.py
from typing import Dict, Any, Tuple, List
from algorithms.base import CSPAlgorithm, register_algorithm
from .config import MEU_ALGORITMO_DEFAULTS
from .implementation import MeuAlgoritmoCore

@register_algorithm
class MeuAlgoritmo(CSPAlgorithm):
    """
    Implementação do Meu Algoritmo para CSP.
    
    Este algoritmo utiliza [descreva a abordagem] para encontrar
    uma solução aproximada do Closest String Problem.
    
    Características:
    - Complexidade: O(n * m * k) onde n=strings, m=length, k=iterations
    - Determinístico: [True/False]
    - Paralelizável: [True/False]
    - Memória: O(n * m)
    
    Referências:
        [1] Autor, A. "Título do Paper". Conferência, Ano.
        [2] Implementação baseada em: [fonte]
    """
    
    # Metadados obrigatórios
    name = "MeuAlgoritmo"
    default_params = MEU_ALGORITMO_DEFAULTS
    is_deterministic = False  # True se for determinístico
    supports_internal_parallel = True  # True se suporta paralelismo interno
    
    def __init__(self, strings: List[str], alphabet: str, **params):
        """
        Inicializa o algoritmo.
        
        Args:
            strings: Lista de strings do dataset
            alphabet: Alfabeto utilizado
            **params: Parâmetros específicos do algoritmo
        """
        super().__init__(strings, alphabet, **params)
        
        # Validar parâmetros específicos
        self._validate_algorithm_params()
        
        # Inicializar o core do algoritmo
        self.core = MeuAlgoritmoCore(
            strings=self.strings,
            alphabet=self.alphabet,
            **self.params
        )
        
        # Estado interno
        self.iteration_count = 0
        self.best_solution = None
        self.convergence_history = []
    
    def _validate_algorithm_params(self) -> None:
        """Valida parâmetros específicos do algoritmo."""
        # Exemplo de validações
        if self.params.get('max_iterations', 0) <= 0:
            raise ValueError("max_iterations deve ser positivo")
        
        if not 0 <= self.params.get('learning_rate', 0) <= 1:
            raise ValueError("learning_rate deve estar entre 0 e 1")
    
    def run(self) -> Tuple[str, int, Dict[str, Any]]:
        """
        Executa o algoritmo principal.
        
        Returns:
            Tuple contendo:
            - str: String centro encontrada
            - int: Distância máxima
            - Dict[str, Any]: Metadados detalhados
        """
        import time
        
        start_time = time.time()
        
        try:
            # Inicialização
            self._report_progress("Inicializando algoritmo...")
            self.core.initialize()
            
            # Loop principal
            self._report_progress("Iniciando iterações...")
            while not self._should_stop():
                # Executar uma iteração
                self.core.iterate()
                self.iteration_count += 1
                
                # Reportar progresso
                if self.iteration_count % 10 == 0:
                    current_best = self.core.get_current_best_distance()
                    self._report_progress(
                        f"Iteração {self.iteration_count}: "
                        f"melhor distância = {current_best}"
                    )
                
                # Verificar convergência
                if self._check_convergence():
                    self._report_progress("Convergência detectada!")
                    break
            
            # Obter resultado final
            center = self.core.get_best_solution()
            distance = self.core.get_best_distance()
            
            # Calcular tempo de execução
            execution_time = time.time() - start_time
            
            # Coletar metadados
            metadata = self._collect_metadata(execution_time)
            
            return center, distance, metadata
            
        except Exception as e:
            self._report_warning(f"Erro durante execução: {e}")
            raise
    
    def _should_stop(self) -> bool:
        """Verifica se o algoritmo deve parar."""
        max_iterations = self.params.get('max_iterations', 1000)
        return self.iteration_count >= max_iterations
    
    def _check_convergence(self) -> bool:
        """Verifica se o algoritmo convergiu."""
        # Implementar lógica de convergência específica
        convergence_threshold = self.params.get('convergence_threshold', 0.001)
        window_size = self.params.get('convergence_window', 10)
        
        # Exemplo: verificar se melhoria nos últimos N iterações < threshold
        if len(self.convergence_history) >= window_size:
            recent_improvements = self.convergence_history[-window_size:]
            avg_improvement = sum(recent_improvements) / window_size
            return avg_improvement < convergence_threshold
        
        return False
    
    def _collect_metadata(self, execution_time: float) -> Dict[str, Any]:
        """Coleta metadados detalhados da execução."""
        return {
            # Métricas básicas
            'execution_time': execution_time,
            'iterations': self.iteration_count,
            'algorithm_name': self.name,
            'parameters_used': self.params.copy(),
            
            # Métricas de convergência
            'convergence_history': self.convergence_history.copy(),
            'converged': self._check_convergence(),
            
            # Métricas específicas do algoritmo
            'custom_metric_1': self.core.get_custom_metric_1(),
            'custom_metric_2': self.core.get_custom_metric_2(),
            
            # Informações do dataset
            'dataset_size': len(self.strings),
            'sequence_length': len(self.strings[0]) if self.strings else 0,
            'alphabet_size': len(self.alphabet),
            
            # Informações de qualidade da solução
            'solution_quality': self.core.assess_solution_quality(),
            'exploration_vs_exploitation': self.core.get_exploration_ratio(),
        }
```

### 2. Implementação do Core (`implementation.py`)

```python
# algorithms/meu_algoritmo/implementation.py
from typing import List, Dict, Any, Optional
import random
import numpy as np
from src.utils.distance import hamming_distance

class MeuAlgoritmoCore:
    """
    Implementação do núcleo do algoritmo.
    
    Separada da classe principal para facilitar testes
    e manutenção do código.
    """
    
    def __init__(self, strings: List[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = params
        
        # Estado interno
        self.current_solution = None
        self.best_solution = None
        self.best_distance = float('inf')
        
        # Estruturas de dados específicas
        self.population = []
        self.fitness_scores = []
        self.generation = 0
        
        # Métricas personalizadas
        self.custom_metrics = {}
    
    def initialize(self) -> None:
        """Inicializa o algoritmo."""
        # Inicializar população/solução inicial
        self._initialize_population()
        
        # Calcular fitness inicial
        self._evaluate_population()
        
        # Definir melhor solução inicial
        self._update_best_solution()
    
    def iterate(self) -> None:
        """Executa uma iteração do algoritmo."""
        # Exemplo de iteração (personalizar para seu algoritmo)
        
        # 1. Seleção
        parents = self._selection()
        
        # 2. Crossover
        offspring = self._crossover(parents)
        
        # 3. Mutação
        mutated_offspring = self._mutation(offspring)
        
        # 4. Avaliação
        self._evaluate_individuals(mutated_offspring)
        
        # 5. Substituição
        self._replacement(mutated_offspring)
        
        # 6. Atualizar melhor solução
        self._update_best_solution()
        
        # 7. Atualizar métricas
        self._update_metrics()
        
        self.generation += 1
    
    def _initialize_population(self) -> None:
        """Inicializa a população."""
        population_size = self.params.get('population_size', 50)
        sequence_length = len(self.strings[0]) if self.strings else 0
        
        self.population = []
        for _ in range(population_size):
            # Criar indivíduo aleatório
            individual = ''.join(
                random.choice(self.alphabet) 
                for _ in range(sequence_length)
            )
            self.population.append(individual)
    
    def _evaluate_population(self) -> None:
        """Avalia toda a população."""
        self.fitness_scores = []
        for individual in self.population:
            fitness = self._calculate_fitness(individual)
            self.fitness_scores.append(fitness)
    
    def _calculate_fitness(self, solution: str) -> float:
        """Calcula fitness de uma solução."""
        # Para CSP, fitness é tipicamente o negativo da distância máxima
        max_distance = max(
            hamming_distance(solution, seq) 
            for seq in self.strings
        )
        return -max_distance  # Negativo para maximização
    
    def _selection(self) -> List[str]:
        """Seleciona pais para reprodução."""
        # Implementar seleção (torneio, roleta, etc.)
        tournament_size = self.params.get('tournament_size', 3)
        parents = []
        
        for _ in range(len(self.population)):
            # Seleção por torneio
            tournament = random.sample(
                list(zip(self.population, self.fitness_scores)),
                tournament_size
            )
            winner = max(tournament, key=lambda x: x[1])
            parents.append(winner[0])
        
        return parents
    
    def _crossover(self, parents: List[str]) -> List[str]:
        """Aplica crossover nos pais."""
        crossover_rate = self.params.get('crossover_rate', 0.8)
        offspring = []
        
        for i in range(0, len(parents), 2):
            parent1 = parents[i]
            parent2 = parents[i + 1] if i + 1 < len(parents) else parents[0]
            
            if random.random() < crossover_rate:
                # Crossover de um ponto
                point = random.randint(1, len(parent1) - 1)
                child1 = parent1[:point] + parent2[point:]
                child2 = parent2[:point] + parent1[point:]
                offspring.extend([child1, child2])
            else:
                offspring.extend([parent1, parent2])
        
        return offspring
    
    def _mutation(self, individuals: List[str]) -> List[str]:
        """Aplica mutação nos indivíduos."""
        mutation_rate = self.params.get('mutation_rate', 0.01)
        mutated = []
        
        for individual in individuals:
            if random.random() < mutation_rate:
                # Mutação de um ponto
                position = random.randint(0, len(individual) - 1)
                new_char = random.choice(self.alphabet)
                mutated_individual = (
                    individual[:position] + 
                    new_char + 
                    individual[position + 1:]
                )
                mutated.append(mutated_individual)
            else:
                mutated.append(individual)
        
        return mutated
    
    def _evaluate_individuals(self, individuals: List[str]) -> List[float]:
        """Avalia lista de indivíduos."""
        return [self._calculate_fitness(ind) for ind in individuals]
    
    def _replacement(self, offspring: List[str]) -> None:
        """Substitui população atual pela descendência."""
        # Estratégia elitista: manter os melhores
        offspring_fitness = self._evaluate_individuals(offspring)
        
        # Combinar populações
        combined_pop = self.population + offspring
        combined_fitness = self.fitness_scores + offspring_fitness
        
        # Ordenar por fitness
        combined = list(zip(combined_pop, combined_fitness))
        combined.sort(key=lambda x: x[1], reverse=True)
        
        # Manter apenas os melhores
        population_size = self.params.get('population_size', 50)
        self.population = [ind for ind, _ in combined[:population_size]]
        self.fitness_scores = [fit for _, fit in combined[:population_size]]
    
    def _update_best_solution(self) -> None:
        """Atualiza a melhor solução encontrada."""
        if self.fitness_scores:
            best_idx = max(range(len(self.fitness_scores)), 
                          key=lambda i: self.fitness_scores[i])
            current_best = self.population[best_idx]
            current_best_distance = -self.fitness_scores[best_idx]
            
            if current_best_distance < self.best_distance:
                self.best_solution = current_best
                self.best_distance = current_best_distance
    
    def _update_metrics(self) -> None:
        """Atualiza métricas personalizadas."""
        # Diversidade da população
        diversity = self._calculate_population_diversity()
        
        # Convergência
        convergence_rate = self._calculate_convergence_rate()
        
        self.custom_metrics.update({
            'diversity': diversity,
            'convergence_rate': convergence_rate,
            'generation': self.generation
        })
    
    def _calculate_population_diversity(self) -> float:
        """Calcula diversidade da população."""
        if len(self.population) < 2:
            return 0.0
        
        total_distance = 0
        count = 0
        
        for i in range(len(self.population)):
            for j in range(i + 1, len(self.population)):
                total_distance += hamming_distance(
                    self.population[i], 
                    self.population[j]
                )
                count += 1
        
        return total_distance / count if count > 0 else 0.0
    
    def _calculate_convergence_rate(self) -> float:
        """Calcula taxa de convergência."""
        # Implementar lógica específica
        return 0.0  # Placeholder
    
    # Métodos públicos para acesso aos resultados
    def get_best_solution(self) -> str:
        """Retorna a melhor solução encontrada."""
        return self.best_solution or ""
    
    def get_best_distance(self) -> int:
        """Retorna a melhor distância encontrada."""
        return int(self.best_distance) if self.best_distance != float('inf') else 0
    
    def get_current_best_distance(self) -> int:
        """Retorna a melhor distância da geração atual."""
        if self.fitness_scores:
            return int(-max(self.fitness_scores))
        return 0
    
    def get_custom_metric_1(self) -> float:
        """Retorna métrica personalizada 1."""
        return self.custom_metrics.get('diversity', 0.0)
    
    def get_custom_metric_2(self) -> float:
        """Retorna métrica personalizada 2."""
        return self.custom_metrics.get('convergence_rate', 0.0)
    
    def assess_solution_quality(self) -> str:
        """Avalia qualidade da solução."""
        if self.best_distance == float('inf'):
            return "no_solution"
        elif self.best_distance <= 5:
            return "excellent"
        elif self.best_distance <= 10:
            return "good"
        elif self.best_distance <= 20:
            return "fair"
        else:
            return "poor"
    
    def get_exploration_ratio(self) -> float:
        """Retorna razão exploração vs exploração."""
        return self.custom_metrics.get('diversity', 0.0) / 100.0
```

### 3. Configuração (`config.py`)

```python
# algorithms/meu_algoritmo/config.py
"""
Configurações padrão para Meu Algoritmo.
"""

MEU_ALGORITMO_DEFAULTS = {
    # Parâmetros básicos
    'population_size': 50,
    'max_iterations': 1000,
    'max_generations': 500,
    
    # Parâmetros genéticos
    'crossover_rate': 0.8,
    'mutation_rate': 0.01,
    'tournament_size': 3,
    
    # Parâmetros de convergência
    'convergence_threshold': 0.001,
    'convergence_window': 20,
    'max_stagnation': 50,
    
    # Parâmetros específicos do algoritmo
    'learning_rate': 0.1,
    'exploration_factor': 0.3,
    'local_search_intensity': 0.2,
    
    # Parâmetros de performance
    'use_parallel': True,
    'num_threads': 4,
    'memory_limit_mb': 512,
    
    # Parâmetros de qualidade
    'elitism_rate': 0.1,
    'diversity_threshold': 0.5,
    'adaptive_parameters': True,
    
    # Parâmetros de debug
    'verbose': False,
    'save_history': False,
    'checkpoint_interval': 100,
}

# Validação de parâmetros
PARAMETER_CONSTRAINTS = {
    'population_size': {'min': 10, 'max': 1000, 'type': int},
    'max_iterations': {'min': 1, 'max': 10000, 'type': int},
    'crossover_rate': {'min': 0.0, 'max': 1.0, 'type': float},
    'mutation_rate': {'min': 0.0, 'max': 1.0, 'type': float},
    'learning_rate': {'min': 0.0, 'max': 1.0, 'type': float},
}

def validate_parameters(params: dict) -> dict:
    """
    Valida e corrige parâmetros do algoritmo.
    
    Args:
        params: Parâmetros a serem validados
        
    Returns:
        Parâmetros validados e corrigidos
        
    Raises:
        ValueError: Se parâmetros são inválidos
    """
    validated = params.copy()
    
    for param, constraints in PARAMETER_CONSTRAINTS.items():
        if param in validated:
            value = validated[param]
            
            # Verificar tipo
            if not isinstance(value, constraints['type']):
                try:
                    validated[param] = constraints['type'](value)
                except (ValueError, TypeError):
                    raise ValueError(f"Parâmetro '{param}' deve ser do tipo {constraints['type'].__name__}")
            
            # Verificar limites
            if 'min' in constraints and validated[param] < constraints['min']:
                raise ValueError(f"Parâmetro '{param}' deve ser >= {constraints['min']}")
            
            if 'max' in constraints and validated[param] > constraints['max']:
                raise ValueError(f"Parâmetro '{param}' deve ser <= {constraints['max']}")
    
    return validated
```

### 4. Testes (`tests/test_algorithm.py`)

```python
# algorithms/meu_algoritmo/tests/test_algorithm.py
import pytest
from unittest.mock import Mock, patch
from algorithms.meu_algoritmo.algorithm import MeuAlgoritmo

class TestMeuAlgoritmo:
    def setup_method(self):
        """Configuração para cada teste."""
        self.sequences = [
            "ACGTACGTACGT",
            "ACGTACGTACGA",
            "ACGTACGTACGC",
            "ACGTACGTACGG"
        ]
        self.alphabet = "ACGT"
        self.params = {
            'population_size': 10,
            'max_iterations': 50,
            'crossover_rate': 0.8,
            'mutation_rate': 0.01
        }
    
    def test_initialization(self):
        """Testa inicialização do algoritmo."""
        algo = MeuAlgoritmo(self.sequences, self.alphabet, **self.params)
        
        assert algo.name == "MeuAlgoritmo"
        assert algo.strings == self.sequences
        assert algo.alphabet == self.alphabet
        assert algo.params['population_size'] == 10
    
    def test_invalid_parameters(self):
        """Testa validação de parâmetros inválidos."""
        with pytest.raises(ValueError, match="max_iterations deve ser positivo"):
            MeuAlgoritmo(self.sequences, self.alphabet, max_iterations=0)
        
        with pytest.raises(ValueError, match="learning_rate deve estar entre 0 e 1"):
            MeuAlgoritmo(self.sequences, self.alphabet, learning_rate=1.5)
    
    def test_run_basic(self):
        """Testa execução básica do algoritmo."""
        algo = MeuAlgoritmo(self.sequences, self.alphabet, **self.params)
        center, distance, metadata = algo.run()
        
        # Verificar tipos de retorno
        assert isinstance(center, str)
        assert isinstance(distance, int)
        assert isinstance(metadata, dict)
        
        # Verificar valores básicos
        assert len(center) == len(self.sequences[0])
        assert distance >= 0
        assert all(c in self.alphabet for c in center)
    
    def test_progress_callbacks(self):
        """Testa callbacks de progresso."""
        progress_messages = []
        warning_messages = []
        
        def progress_callback(msg):
            progress_messages.append(msg)
        
        def warning_callback(msg):
            warning_messages.append(msg)
        
        algo = MeuAlgoritmo(self.sequences, self.alphabet, **self.params)
        algo.set_progress_callback(progress_callback)
        algo.set_warning_callback(warning_callback)
        
        algo.run()
        
        # Verificar se callbacks foram chamados
        assert len(progress_messages) > 0
        assert any("Inicializando" in msg for msg in progress_messages)
    
    def test_metadata_completeness(self):
        """Testa completude dos metadados."""
        algo = MeuAlgoritmo(self.sequences, self.alphabet, **self.params)
        _, _, metadata = algo.run()
        
        # Verificar metadados obrigatórios
        required_keys = [
            'execution_time', 'iterations', 'algorithm_name',
            'parameters_used', 'dataset_size', 'sequence_length'
        ]
        
        for key in required_keys:
            assert key in metadata, f"Metadado '{key}' não encontrado"
        
        # Verificar tipos
        assert isinstance(metadata['execution_time'], (int, float))
        assert isinstance(metadata['iterations'], int)
        assert isinstance(metadata['algorithm_name'], str)
        assert isinstance(metadata['parameters_used'], dict)
    
    def test_deterministic_behavior(self):
        """Testa comportamento determinístico (se aplicável)."""
        if MeuAlgoritmo.is_deterministic:
            algo1 = MeuAlgoritmo(self.sequences, self.alphabet, **self.params)
            algo2 = MeuAlgoritmo(self.sequences, self.alphabet, **self.params)
            
            result1 = algo1.run()
            result2 = algo2.run()
            
            assert result1[0] == result2[0]  # Mesmo centro
            assert result1[1] == result2[1]  # Mesma distância
    
    def test_empty_sequences(self):
        """Testa comportamento com sequências vazias."""
        with pytest.raises(ValueError):
            MeuAlgoritmo([], self.alphabet, **self.params)
    
    def test_single_sequence(self):
        """Testa comportamento com uma única sequência."""
        single_seq = [self.sequences[0]]
        algo = MeuAlgoritmo(single_seq, self.alphabet, **self.params)
        center, distance, metadata = algo.run()
        
        assert center == single_seq[0]
        assert distance == 0
    
    def test_convergence_detection(self):
        """Testa detecção de convergência."""
        # Usar parâmetros que forçam convergência rápida
        fast_params = {
            **self.params,
            'convergence_threshold': 0.1,
            'convergence_window': 5
        }
        
        algo = MeuAlgoritmo(self.sequences, self.alphabet, **fast_params)
        _, _, metadata = algo.run()
        
        # Verificar se convergência foi detectada
        assert 'converged' in metadata
        assert isinstance(metadata['converged'], bool)
    
    @patch('algorithms.meu_algoritmo.implementation.random.choice')
    def test_mocked_randomness(self, mock_choice):
        """Testa com randomness mockada."""
        mock_choice.return_value = 'A'
        
        algo = MeuAlgoritmo(self.sequences, self.alphabet, **self.params)
        center, distance, metadata = algo.run()
        
        # Verificar que o mock foi usado
        assert mock_choice.called
        assert 'A' in center
```

### 5. Documentação (`README.md`)

```markdown
# Meu Algoritmo

## Visão Geral

Este algoritmo implementa uma solução personalizada para o Closest String Problem baseada em [descreva a abordagem].

## Características

- **Complexidade**: O(n * m * k)
- **Memória**: O(n * m)
- **Determinístico**: Não
- **Paralelizável**: Sim
- **Adequado para**: Instâncias de tamanho médio a grande

## Parâmetros

### Parâmetros Principais

- `population_size` (int, padrão: 50): Tamanho da população
- `max_iterations` (int, padrão: 1000): Número máximo de iterações
- `crossover_rate` (float, padrão: 0.8): Taxa de crossover
- `mutation_rate` (float, padrão: 0.01): Taxa de mutação

### Parâmetros Avançados

- `learning_rate` (float, padrão: 0.1): Taxa de aprendizado
- `convergence_threshold` (float, padrão: 0.001): Threshold de convergência
- `use_parallel` (bool, padrão: True): Usar processamento paralelo

## Exemplo de Uso

```python
from algorithms.meu_algoritmo.algorithm import MeuAlgoritmo

# Dados de exemplo
sequences = ["ACGTACGT", "ACGTACGA", "ACGTACGC"]
alphabet = "ACGT"

# Criar instância
algo = MeuAlgoritmo(sequences, alphabet, 
                   population_size=100,
                   max_iterations=500)

# Executar
center, distance, metadata = algo.run()

print(f"Centro: {center}")
print(f"Distância: {distance}")
print(f"Iterações: {metadata['iterations']}")
```

## Referências

1. [Referência do algoritmo original]
2. [Implementação baseada em]
3. [Papers relacionados]
```

## 🧪 Testando o Algoritmo

### Testes Unitários

```bash
# Executar testes do algoritmo
pytest algorithms/meu_algoritmo/tests/ -v

# Executar com cobertura
pytest algorithms/meu_algoritmo/tests/ --cov=algorithms.meu_algoritmo --cov-report=html
```

### Testes de Integração

```bash
# Testar integração com o sistema
python main.py --silent --dataset synthetic --algorithms MeuAlgoritmo --num-execs 3
```

### Benchmarking

```python
import time
from algorithms.meu_algoritmo.algorithm import MeuAlgoritmo

def benchmark_algorithm():
    sequences = ["ACGTACGT"] * 100  # Dataset grande
    alphabet = "ACGT"
    
    algo = MeuAlgoritmo(sequences, alphabet)
    
    start = time.time()
    center, distance, metadata = algo.run()
    end = time.time()
    
    print(f"Tempo: {end - start:.2f}s")
    print(f"Distância: {distance}")
    print(f"Iterações: {metadata['iterations']}")
```

## 📊 Métricas e Avaliação

### Métricas Coletadas

O algoritmo coleta automaticamente:

- **Performance**: Tempo de execução, iterações
- **Qualidade**: Distância final, taxa de convergência
- **Diversidade**: Diversidade populacional
- **Recursos**: Uso de memória, CPU

### Comparação com Outros Algoritmos

```python
from src.core.comparison import compare_algorithms

results = compare_algorithms(
    algorithms=['Baseline', 'BLF-GA', 'MeuAlgoritmo'],
    datasets=['synthetic_small', 'synthetic_large'],
    num_executions=10
)

print(results.summary())
```

## 🔧 Otimização e Tuning

### Otimização Automática

```python
from src.optimization.optimizer import AlgorithmOptimizer

optimizer = AlgorithmOptimizer(
    algorithm='MeuAlgoritmo',
    dataset='synthetic_medium',
    parameter_space={
        'population_size': (20, 200),
        'mutation_rate': (0.001, 0.1),
        'learning_rate': (0.01, 0.5)
    }
)

best_params = optimizer.optimize(n_trials=100)
print(f"Melhores parâmetros: {best_params}")
```

### Análise de Sensibilidade

```python
from src.optimization.sensitivity import sensitivity_analysis

results = sensitivity_analysis(
    algorithm='MeuAlgoritmo',
    dataset='synthetic_medium',
    parameter='population_size',
    values=[10, 25, 50, 100, 200]
)

results.plot()
```

## 📈 Monitoramento e Debugging

### Logs Detalhados

```python
import logging

# Configurar logging detalhado
logging.basicConfig(level=logging.DEBUG)

algo = MeuAlgoritmo(sequences, alphabet, verbose=True)
result = algo.run()
```

### Visualização de Convergência

```python
def plot_convergence(metadata):
    import matplotlib.pyplot as plt
    
    history = metadata.get('convergence_history', [])
    plt.plot(history)
    plt.xlabel('Iteração')
    plt.ylabel('Melhor Distância')
    plt.title('Convergência do Algoritmo')
    plt.show()
```

## 🚀 Próximos Passos

1. **Implementar seu algoritmo** seguindo este guia
2. **Adicionar testes** abrangentes
3. **Documentar** adequadamente
4. **Otimizar** parâmetros
5. **Comparar** com algoritmos existentes
6. **Publicar** resultados

## 💡 Dicas Importantes

- **Sempre validar** parâmetros de entrada
- **Implementar callbacks** de progresso
- **Coletar métricas** detalhadas
- **Tratar erros** adequadamente
- **Documentar** comportamento
- **Testar** exaustivamente

---

Com este guia, você deve ser capaz de implementar e integrar novos algoritmos ao CSP-BLFGA de forma padronizada e robusta. 🎯
