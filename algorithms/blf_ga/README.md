# BLF-GA: Blockwise Learning Fusion + Genetic Algorithm

## Visão Geral
O BLF-GA é uma metaheurística híbrida que combina **aprendizado por blocos** com **algoritmo genético** para resolver o Closest String Problem. A abordagem divide o problema em segmentos menores, aplica aprendizado local em cada bloco e usa evolução genética para otimização global.

## Funcionamento

### Arquitetura Híbrida
O algoritmo opera em **três camadas principais**:
1. **B-Splitter**: Divisão adaptativa em blocos
2. **Block Learning**: Aprendizado local por segmento  
3. **Genetic Evolution**: Evolução populacional global

### Fluxo de Execução Detalhado

#### 1. Inicialização
- **População Inicial**: Cria população diversificada incluindo consenso ganancioso
- **Blocking Inicial**: Divide strings em √L blocos aproximadamente iguais
- **Avaliação Base**: Calcula fitness inicial (distância máxima) para cada indivíduo

#### 2. Ciclo Evolutivo Principal
```
Para cada geração g = 1 até max_gens:
  ├── Block Learning Phase
  ├── Genetic Operations  
  ├── Population Update
  ├── Adaptive Blocking
  └── Convergence Check
```

#### 3. Block Learning Phase
- **Análise por Bloco**: Para cada bloco [i,j), identifica padrões locais
- **Consensus Learning**: Gera consenso ótimo para o bloco atual
- **Pattern Mining**: Extrai subsequências promissoras da população
- **Repository Update**: Armazena melhores blocos encontrados

#### 4. Genetic Operations
- **Seleção por Torneio**: Escolhe pais via competição (k=3)
- **Crossover Blocks**: Troca blocos inteiros entre pais
- **Mutação Adaptativa**: Modifica posições com taxa baseada na diversidade
- **Elitismo**: Preserva melhores indivíduos (5% da população)

#### 5. Adaptive Blocking
- **Entropy Analysis**: Calcula entropia por posição
- **Dynamic Segmentation**: Redefine blocos com base na variabilidade
- **Block Refinement**: Ajusta tamanho de blocos conforme necessário

## Heurísticas Principais

### Learning Fusion
**Princípio**: Combinar conhecimento local (blocos) com busca global (GA)
- **Bloco Consensus**: Cada bloco aprende independentemente seu melhor padrão
- **Cross-Block Exchange**: Blocos bons são propagados entre indivíduos
- **Global Optimization**: GA refina a combinação de blocos

### Adaptive Blocking Strategy
**Critério de Divisão**: Balancear granularidade vs. eficiência
- **Inicial**: √L blocos para começar com granularidade moderada
- **Entropy-Based**: Redividir baseado na variabilidade observada
- **Minimum Size**: Manter blocos ≥ 3 posições para aprendizado efetivo

### Evolutionary Pressure
**Multi-Level Selection**: Pressão seletiva em diferentes granularidades
- **Individual Level**: Seleção baseada em fitness total
- **Block Level**: Blocos bons se espalham independentemente
- **Population Level**: Rediversificação periódica para evitar convergência prematura

### Elite Refinement
**Local Search Intensification**: Melhoria iterativa dos melhores
- **Hillclimbing**: Tenta melhorar posição por posição
- **Neighborhood Search**: Explora vizinhança de Hamming
- **Greedy Refinement**: Aplica mudanças que reduzem distância máxima

## Parâmetros de Configuração

### População e Evolução
- **pop_size**: Tamanho da população (padrão: 100)
- **max_gens**: Gerações máximas (padrão: 100)  
- **elite_rate**: Taxa de elitismo (padrão: 5%)

### Operadores Genéticos
- **cross_prob**: Probabilidade de crossover (padrão: 90%)
- **mut_prob**: Probabilidade de mutação (padrão: 80%)
- **tournament_size**: Tamanho do torneio (fixo: 3)

### Blocking e Aprendizado
- **initial_blocks**: Blocos iniciais (padrão: 10)
- **min_block_len**: Tamanho mínimo de bloco (padrão: 3)
- **rediv_freq**: Frequência de redivisão (padrão: 10 gerações)

### Controle de Execução
- **max_time**: Tempo limite em segundos (padrão: 600s)
- **seed**: Semente para reprodutibilidade

## Características Avançadas

### Convergence Management
- **Diversity Tracking**: Monitora diversidade populacional
- **Rediversification**: Introduz novos indivíduos quando necessário
- **Early Stopping**: Para quando não há melhoria por N gerações

### Memory Management  
- **Block Repository**: Cache de melhores blocos descobertos
- **Elite History**: Mantém histórico de melhores soluções
- **Pattern Database**: Armazena padrões recorrentes promissores

### Parallel-Ready Design
- **Block Independence**: Blocos podem ser processados paralelamente
- **Population Chunks**: População pode ser dividida para processamento paralelo
- **Fitness Caching**: Evita recálculos desnecessários

## Performance e Escalabilidade

### Complexidade Temporal
- **Por Geração**: O(pop_size × n × L + block_operations)
- **Total**: O(max_gens × pop_size × n × L)
- **Block Learning**: O(num_blocks × L/num_blocks × search_complexity)

### Cenários de Uso Ideal
- **Instâncias Médias a Grandes**: n > 20, L > 100
- **Dados com Estrutura Local**: Padrões regionais conservados
- **Tempo Disponível**: Quando se pode aguardar otimização iterativa
- **Qualidade Prioritária**: Quando resultado próximo do ótimo é crucial

### Vantagens Competitivas
- **Hibridização Efetiva**: Combina busca local e global
- **Adaptabilidade**: Ajusta estratégia durante execução
- **Robustez**: Funciona bem em diferentes tipos de instância
- **Escalabilidade**: Performance se mantém em problemas maiores

### Limitações
- **Complexidade Paramétrica**: Muitos parâmetros para ajustar
- **Tempo de Execução**: Pode ser lento para instâncias muito pequenas
- **Convergência**: Pode estagnar em ótimos locais
- **Memória**: Usa mais memória que algoritmos simples

O BLF-GA representa uma abordagem sofisticada que combina o melhor dos mundos do aprendizado local e otimização global, sendo especialmente efetivo para instâncias complexas onde métodos mais simples falham.
