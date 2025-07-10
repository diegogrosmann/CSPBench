# BLF-GA: Blockwise Learning Fusion + Genetic Algorithm

O **BLF-GA** é uma metaheurística híbrida avançada que combina aprendizado por blocos (Blockwise Learning) com algoritmo genético global, oferecendo uma abordagem sofisticada para resolver o Closest String Problem.

## 🧬 Visão Geral

### **Arquitetura Híbrida**
O BLF-GA opera em múltiplas camadas, combinando:
- **Aprendizado por Blocos**: Otimização local de segmentos das strings
- **Evolução Genética**: Busca global através de população evolucionária
- **Fusão Adaptativa**: Combinação inteligente de conhecimento local e global
- **Refinamento Elite**: Busca local intensiva nos melhores indivíduos

### **Fases do Algoritmo**
1. **Inicialização**: População inicial com diversidade controlada
2. **Divisão em Blocos**: Segmentação adaptativa das strings
3. **Aprendizado Local**: Otimização por bloco usando consenso e busca local
4. **Evolução Global**: Operadores genéticos (seleção, crossover, mutação)
5. **Fusão de Conhecimento**: Combinação de aprendizado local e global
6. **Refinamento Elite**: Busca local nos melhores indivíduos
7. **Redivisão Dinâmica**: Reconfiguração de blocos baseada na evolução

## 🏗️ Componentes Técnicos

### **Sistema de Blocos**
- **B-Splitter**: Divisão inteligente em blocos contíguos
- **Tamanho Adaptativo**: Blocos ajustados baseado no progresso
- **Redivisão Dinâmica**: Reconfiguração periódica para escape de ótimos locais
- **Aprendizado por Bloco**: Otimização local especializada

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
"elite_rate": 0.1,            # Taxa de elitismo (10%)
"tournament_size": 3,          # Tamanho do torneio
```

### **Operadores Genéticos**
```python
"crossover_prob": 0.8,         # Probabilidade de crossover
"mutation_prob": 0.1,          # Probabilidade de mutação base
"adaptive_mutation": True,      # Mutação adaptativa
"local_search_prob": 0.3,      # Probabilidade de busca local
```

### **Sistema de Blocos**
```python
"initial_blocks": 4,           # Número inicial de blocos
"min_block_length": 3,         # Tamanho mínimo de bloco
"redivision_frequency": 50,    # Frequência de redivisão
"block_learning_rate": 0.1,    # Taxa de aprendizado por bloco
```

### **Controle e Convergência**
```python
"max_time": 300,              # Tempo máximo (segundos)
"convergence_generations": 50, # Gerações sem melhoria para parar
"diversity_threshold": 0.1,    # Limiar de diversidade mínima
"seed": None,                 # Semente para reprodutibilidade
```

## 🎯 Estratégias e Heurísticas

### **Blockwise Learning**
- **Consenso Local**: Geração de consenso ótimo por bloco
- **Busca Exaustiva Local**: Exploração completa para blocos pequenos
- **Hill Climbing**: Refinamento local por substituição de símbolos
- **Cache de Blocos**: Reutilização de soluções de blocos similares

### **Fusão de Conhecimento**
- **Voting Scheme**: Combinação ponderada de soluções
- **Block Replacement**: Substituição de blocos baseada em qualidade
- **Hybrid Offspring**: Geração de descendentes híbridos
- **Knowledge Transfer**: Transferência entre gerações

### **Adaptação Dinâmica**
- **Population Diversity Control**: Manutenção de diversidade genética
- **Operator Rate Adaptation**: Ajuste automático de taxas
- **Block Size Adaptation**: Redimensionamento baseado em performance
- **Search Strategy Switching**: Alternância entre estratégias

## 💻 Exemplo de Uso

### **Uso Básico**
```python
from algorithms.blf_ga.algorithm import BLFGAAlgorithm

# Dataset de exemplo
strings = ["ACGTACGTACGT", "AGGTACGTAAGT", "ACGTAAGTTCGT"]
alphabet = "ACGT"

# Configurar algoritmo
algorithm = BLFGAAlgorithm(
    strings, alphabet,
    population_size=100,
    max_generations=200,
    crossover_prob=0.8,
    mutation_prob=0.1
)

# Executar
center, distance, metadata = algorithm.run()

print(f"Centro encontrado: {center}")
print(f"Distância: {distance}")
print(f"Gerações: {metadata['generations_run']}")
print(f"Tempo: {metadata['execution_time']:.2f}s")
```

### **Configuração Avançada**
```python
# Configuração para instâncias grandes
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
# Execução básica
python main.py --algorithms BLF-GA --dataset synthetic

# Com parâmetros customizados
python main.py --algorithms BLF-GA --dataset synthetic --workers 4

# Execução otimizada para instâncias grandes
python main.py --algorithms BLF-GA --dataset file --timeout 600
```

### **Configuração YAML**
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

## 📈 Performance e Características

### **Complexidade Computacional**
- **Temporal**: O(G × P × L × n) onde:
  - G: número de gerações
  - P: tamanho da população
  - L: comprimento das strings
  - n: número de strings
- **Espacial**: O(P × L) para população + O(B × L) para blocos

### **Paralelização**
- ✅ **Suporte Interno**: `supports_internal_parallel = True`
- ✅ **Avaliação Paralela**: População avaliada em paralelo
- ✅ **Blocos Paralelos**: Processamento simultâneo de blocos
- ✅ **Auto-configuração**: Workers ajustados automaticamente

### **Escalabilidade**
- **Instâncias Pequenas** (n≤20, L≤50): ~10-30s
- **Instâncias Médias** (n≤100, L≤200): ~1-5 min
- **Instâncias Grandes** (n≤500, L≤1000): ~10-30 min

## 🎯 Casos de Uso

### **✅ Ideal Para**
- **Instâncias Médias/Grandes**: n > 20, L > 100
- **Dados com Padrões Locais**: Sequências biológicas estruturadas
- **Qualidade Prioritária**: Quando solução próxima do ótimo é essencial
- **Recursos Computacionais Disponíveis**: Sistemas multi-core
- **Execução em Lote**: Múltiplas execuções com estatísticas

### **❌ Limitações**
- **Complexidade de Configuração**: Muitos parâmetros para ajustar
- **Tempo de Execução**: Pode ser lento para instâncias pequenas
- **Variabilidade**: Resultados podem variar entre execuções
- **Recursos**: Consome mais CPU e memória que algoritmos simples

## 🔬 Metadados Coletados

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

## 🧪 Configurações Recomendadas

### **Para Instâncias Pequenas**
```python
{
    "population_size": 50,
    "max_generations": 100,
    "initial_blocks": 2,
    "max_time": 60
}
```

### **Para Instâncias Médias**
```python
{
    "population_size": 100,
    "max_generations": 300,
    "initial_blocks": 4,
    "max_time": 300
}
```

### **Para Instâncias Grandes**
```python
{
    "population_size": 200,
    "max_generations": 500,
    "initial_blocks": 8,
    "max_time": 900
}
```

### **Para Execução Rápida**
```python
{
    "population_size": 30,
    "max_generations": 50,
    "convergence_generations": 20,
    "max_time": 30
}
```

## 🎨 Análise e Visualizações

### **Curva de Convergência**
```python
import matplotlib.pyplot as plt

def plot_convergence(metadata):
    generations = range(len(metadata['convergence_curve']))
    fitness = metadata['convergence_curve']
    
    plt.figure(figsize=(10, 6))
    plt.plot(generations, fitness, 'b-', linewidth=2)
    plt.xlabel('Geração')
    plt.ylabel('Melhor Fitness')
    plt.title('Convergência do BLF-GA')
    plt.grid(True, alpha=0.3)
    plt.show()
```

### **Diversidade Populacional**
```python
def plot_diversity(metadata):
    generations = range(len(metadata['diversity_curve']))
    diversity = metadata['diversity_curve']
    
    plt.figure(figsize=(10, 6))
    plt.plot(generations, diversity, 'r-', linewidth=2)
    plt.xlabel('Geração')
    plt.ylabel('Diversidade')
    plt.title('Diversidade Populacional')
    plt.grid(True, alpha=0.3)
    plt.show()
```

## 🔗 Integração com CSPBench

### **Recursos do Framework**
- **Registro Automático**: `@register_algorithm`
- **Paralelismo Interno**: Configuração automática de workers
- **Callbacks de Progresso**: Relatórios em tempo real
- **Timeouts**: Controle de tempo máximo
- **Monitoramento**: Interface curses compatível

### **Otimização de Hiperparâmetros**
```yaml
# Configuração para Optuna
optimization_config:
  param_space:
    "BLF-GA":
      population_size: ["int", 50, 300]
      max_generations: ["int", 100, 500]
      crossover_prob: ["uniform", 0.6, 0.9]
      mutation_prob: ["uniform", 0.05, 0.2]
      initial_blocks: ["int", 2, 10]
```

### **Análise de Sensibilidade**
```yaml
# Configuração para SALib
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

## 🚀 Extensões e Melhorias

### **Versões Futuras**
1. **BLF-GA Multi-Objetivo**: Otimização simultânea de múltiplos critérios
2. **BLF-GA Distribuído**: Execução em cluster/grid
3. **BLF-GA Adaptativo**: Algoritmo que se auto-ajusta
4. **BLF-GA Quântico**: Inspiração em computação quântica

### **Melhorias Implementáveis**
- **Cache Inteligente**: Reutilização de avaliações
- **Operadores Especializados**: Crossover e mutação específicos para CSP
- **Aprendizado Online**: Adaptação baseada em histórico
- **Paralelização Avançada**: GPU computing

---

*BLF-GA: Quando excelência em qualidade e sofisticação algorítmica são prioridades para resolver o CSP.*

## Documentação

- Consulte o código para docstrings detalhadas (Google style).
- Integração automática com o framework CSP via decorador `@register_algorithm`.
