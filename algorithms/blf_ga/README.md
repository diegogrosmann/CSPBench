# BLF-GA: Blockwise Learning Fusion + Genetic Algorithm

The **BLF-GA** is an advanced hybrid metaheuristic that combines blockwise learning with global genetic algorithm, offering a sophisticated approach to solve the Closest String Problem.

## 🧬 Overview

### **Hybrid Architecture**
BLF-GA operates on multiple layers, combining:
- **Blockwise Learning**: Local optimization of string segments
- **Genetic Evolution**: Global search through evolutionary population
- **Adaptive Fusion**: Intelligent combination of local and global knowledge
- **Elite Refinement**: Intensive local search on best individuals

### **Algorithm Phases**
1. **Initialization**: Initial population with controlled diversity
2. **Block Division**: Adaptive string segmentation
3. **Local Learning**: Block-wise optimization using consensus and local search
4. **Global Evolution**: Genetic operators (selection, crossover, mutation)
5. **Knowledge Fusion**: Combination of local and global learning
6. **Elite Refinement**: Local search on best individuals
7. **Dynamic Redivision**: Block reconfiguration based on evolution

## 🏗️ Technical Components

### **Block System**
- **B-Splitter**: Intelligent division into contiguous blocks
- **Adaptive Size**: Blocks adjusted based on progress
- **Dynamic Redivision**: Periodic reconfiguration to escape local optima
- **Block Learning**: Specialized local optimization

### **Algoritmo Genético**
- **População Diversa**: Inicialização garantindo diversidade genética
- **Seleção por Torneio**: Pressão seletiva balanceada
- **Crossover Especializado**: Operadores adaptados para strings
- **Mutação Inteligente**: Taxa adaptativa baseada na diversidade
- **Elitismo Controlado**: Preservação dos melhores indivíduos

### **Mecanismos Adaptativos**
- **Controle de Diversidade**: Monitoramento e manutenção da diversidade populacional
- **Taxas Dinâmicas**: Ajuste automático de crossover e mutação
- **Critérios de Convergência**: Múltiplos critérios de parada
- **Intensificação/Diversificação**: Balanceamento automático

## ⚙️ Parâmetros Principais

### **População e Evolução**
```python
"population_size": 100,        # Tamanho da população
"max_generations": 300,        # Máximo de gerações
"elite_rate": 0.1,            # Elitism rate (10%)
"tournament_size": 3,          # Tournament size
```

### **Genetic Operators**
```python
"crossover_prob": 0.8,         # Crossover probability
"mutation_prob": 0.1,          # Base mutation probability
"adaptive_mutation": True,      # Adaptive mutation
"local_search_prob": 0.3,      # Local search probability
```

### **Block System**
```python
"initial_blocks": 4,           # Initial number of blocks
"min_block_length": 3,         # Minimum block size
"redivision_frequency": 50,    # Redivision frequency
"block_learning_rate": 0.1,    # Learning rate per block
```

### **Control and Convergence**
```python
"max_time": 300,              # Maximum time (seconds)
"convergence_generations": 50, # Generations without improvement to stop
"diversity_threshold": 0.1,    # Minimum diversity threshold
"seed": None,                 # Seed for reproducibility
```

## 🎯 Strategies and Heuristics

### **Blockwise Learning**
- **Local Consensus**: Optimal consensus generation per block
- **Local Exhaustive Search**: Complete exploration for small blocks
- **Hill Climbing**: Local refinement by symbol substitution
- **Block Cache**: Reuse of solutions for similar blocks

### **Knowledge Fusion**
- **Voting Scheme**: Weighted combination of solutions
- **Block Replacement**: Quality-based block substitution
- **Hybrid Offspring**: Generation of hybrid descendants
- **Knowledge Transfer**: Transfer between generations

### **Dynamic Adaptation**
- **Population Diversity Control**: Genetic diversity maintenance
- **Operator Rate Adaptation**: Automatic rate adjustment
- **Block Size Adaptation**: Performance-based resizing
- **Search Strategy Switching**: Strategy alternation

## 💻 Usage Example

### **Basic Usage**
```python
from algorithms.blf_ga.algorithm import BLFGAAlgorithm

# Example dataset
strings = ["ACGTACGTACGT", "AGGTACGTAAGT", "ACGTAAGTTCGT"]
alphabet = "ACGT"

# Configure algorithm
algorithm = BLFGAAlgorithm(
    strings, alphabet,
    population_size=100,
    max_generations=200,
    crossover_prob=0.8,
    mutation_prob=0.1
)

# Execute
center, distance, metadata = algorithm.run()

print(f"Center found: {center}")
print(f"Distance: {distance}")
print(f"Generations: {metadata['generations_run']}")
print(f"Time: {metadata['execution_time']:.2f}s")
```

### **Advanced Configuration**
```python
# Configuration for large instances
algorithm = BLFGAAlgorithm(
    strings, alphabet,
    population_size=200,
    max_generations=500,
    initial_blocks=8,
    redivision_frequency=100,
    elite_rate=0.15,
    adaptive_mutation=True,
    max_time=600
)
```

### **Via Framework**
```bash
# Basic execution
python main.py --algorithms BLF-GA --dataset synthetic

# With custom parameters
python main.py --algorithms BLF-GA --dataset synthetic --workers 4

# Optimized execution for large instances
python main.py --algorithms BLF-GA --dataset file --timeout 600
```

### **YAML Configuration**
```yaml
algorithms: ["BLF-GA"]
algorithm_params:
  "BLF-GA":
    population_size: 150
    max_generations: 400
    crossover_prob: 0.85
    mutation_prob: 0.12
    initial_blocks: 6
    max_time: 450
```

## 📈 Performance and Characteristics

### **Computational Complexity**
- **Time**: O(G × P × L × n) where:
  - G: number of generations
  - P: population size
  - L: string length
  - n: number of strings
- **Space**: O(P × L) for population + O(B × L) for blocks

### **Parallelization**
- ✅ **Internal Support**: `supports_internal_parallel = True`
- ✅ **Parallel Evaluation**: Population evaluated in parallel
- ✅ **Parallel Blocks**: Simultaneous block processing
- ✅ **Auto-configuration**: Workers adjusted automatically

### **Scalability**
- **Small Instances** (n≤20, L≤50): ~10-30s
- **Medium Instances** (n≤100, L≤200): ~1-5 min
- **Large Instances** (n≤500, L≤1000): ~10-30 min

## 🎯 Use Cases

### **✅ Ideal For**
- **Medium/Large Instances**: n > 20, L > 100
- **Data with Local Patterns**: Structured biological sequences
- **Quality Priority**: When near-optimal solution is essential
- **Available Computational Resources**: Multi-core systems
- **Batch Execution**: Multiple runs with statistics

### **❌ Limitations**
- **Configuration Complexity**: Many parameters to tune
- **Execution Time**: Can be slow for small instances
- **Variability**: Results may vary between runs
- **Resources**: Consumes more CPU and memory than simple algorithms

## 🔬 Collected Metadata

```python
{
    "generations_run": 245,
    "convergence_generation": 201,
    "final_population_diversity": 0.85,
    "best_fitness_evolution": [25, 22, 18, 15, 13],
    "block_divisions": [4, 6, 4, 8],
    "redivisions_count": 4,
    "elite_preserved": 10,
    "crossovers_performed": 1960,
    "mutations_performed": 196,
    "local_searches_performed": 588,
    "execution_time": 124.5,
    "memory_peak": 256.3,
    "convergence_curve": [...],
    "diversity_curve": [...]
}
```

## 🧪 Recommended Configurations

### **For Small Instances**
```python
{
    "population_size": 50,
    "max_generations": 100,
    "initial_blocks": 2,
    "max_time": 60
}
```

### **For Medium Instances**
```python
{
    "population_size": 100,
    "max_generations": 300,
    "initial_blocks": 4,
    "max_time": 300
}
```

### **For Large Instances**
```python
{
    "population_size": 200,
    "max_generations": 500,
    "initial_blocks": 8,
    "max_time": 900
}
```

### **For Fast Execution**
```python
{
    "population_size": 30,
    "max_generations": 50,
    "convergence_generations": 20,
    "max_time": 30
}
```

## 🎨 Analysis and Visualizations

### **Convergence Curve**
```python
import matplotlib.pyplot as plt

def plot_convergence(metadata):
    generations = range(len(metadata['convergence_curve']))
    fitness = metadata['convergence_curve']
    
    plt.figure(figsize=(10, 6))
    plt.plot(generations, fitness, 'b-', linewidth=2)
    plt.xlabel('Generation')
    plt.ylabel('Best Fitness')
    plt.title('BLF-GA Convergence')
    plt.grid(True, alpha=0.3)
    plt.show()
```

### **Population Diversity**
```python
def plot_diversity(metadata):
    generations = range(len(metadata['diversity_curve']))
    diversity = metadata['diversity_curve']
    
    plt.figure(figsize=(10, 6))
    plt.plot(generations, diversity, 'r-', linewidth=2)
    plt.xlabel('Generation')
    plt.ylabel('Diversity')
    plt.title('Population Diversity')
    plt.grid(True, alpha=0.3)
    plt.show()
```

## 🔗 CSPBench Integration

### **Framework Features**
- **Automatic Registration**: `@register_algorithm`
- **Internal Parallelism**: Automatic worker configuration
- **Progress Callbacks**: Real-time reports
- **Timeouts**: Maximum time control
- **Monitoring**: Curses interface compatible

### **Hyperparameter Optimization**
```yaml
# Configuration for Optuna
optimization_config:
  param_space:
    "BLF-GA":
      population_size: ["int", 50, 300]
      max_generations: ["int", 100, 500]
      crossover_prob: ["uniform", 0.6, 0.9]
      mutation_prob: ["uniform", 0.05, 0.2]
      initial_blocks: ["int", 2, 10]
```

### **Sensitivity Analysis**
```yaml
# Configuration for SALib
sensitivity_config:
  param_space:
    "BLF-GA": [
      "population_size",
      "max_generations", 
      "crossover_prob",
      "mutation_prob",
      "elite_rate"
    ]
```

## 🚀 Extensions and Improvements

### **Future Versions**
1. **Multi-Objective BLF-GA**: Simultaneous optimization of multiple criteria
2. **Distributed BLF-GA**: Cluster/grid execution
3. **Adaptive BLF-GA**: Self-adjusting algorithm
4. **Quantum BLF-GA**: Quantum computing inspiration

### **Implementable Improvements**
- **Smart Cache**: Evaluation reuse
- **Specialized Operators**: CSP-specific crossover and mutation
- **Online Learning**: History-based adaptation
- **Advanced Parallelization**: GPU computing

---

*BLF-GA: When excellence in quality and algorithmic sophistication are priorities for solving CSP.*

## Documentation

- See code for detailed docstrings (Google style).
- Automatic integration with CSP framework via `@register_algorithm` decorator.
