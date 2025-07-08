# BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) — Documentação Técnica

## Visão Geral
O **BLF-GA** é uma metaheurística híbrida avançada para o Closest String Problem (CSP), combinando aprendizado por blocos, algoritmos genéticos e mecanismos adaptativos. Ele é projetado para instâncias médias e grandes, onde métodos exatos são inviáveis.

## Arquitetura do Algoritmo

### 1. Blockwise Learning (Aprendizado por Blocos)
- Divide as strings em blocos adaptativos (baseado em entropia)
- Aprende padrões locais (consenso por bloco)
- Cria repositório de conhecimento local

### 2. Genetic Algorithm (Algoritmo Genético)
- População inicial: consenso, variações e aleatórios
- Seleção por torneio
- Crossover: one-point, uniforme, blend-blocks
- Mutação: multi, inversão, transposição
- Elitismo adaptativo

### 3. Fusion (Fusão)
- Combina conhecimento local e global
- Refinamento local dos melhores indivíduos (hill-climbing)
- Mecanismos adaptativos: redivisão de blocos, imigrantes, mutação adaptativa

## Fluxo Geral
1. Inicialização da população e blocos
2. Loop principal:
   - Aprendizado por blocos
   - Evolução genética
   - Mecanismos adaptativos
   - Critérios de parada (ótimo, tempo, estagnação)
3. Retorno da melhor solução

## Complexidade
- **Tempo:** O(G × pop_size × L × |Σ|), onde G = número de gerações
- **Espaço:** O(pop_size × L)

## Pontos Fortes
- Altamente adaptativo e robusto
- Suporte a paralelismo interno
- Diversidade de operadores genéticos
- Refinamento local intensivo

## Limitações
- Estocástico (resultados podem variar)
- Muitos parâmetros para ajuste fino
- Overhead para instâncias pequenas

## Exemplo de Uso
```python
from algorithms.blf_ga.algorithm import BLFGAAlgorithm
alg = BLFGAAlgorithm(strings, alphabet, pop_size=100, generations=200)
center, dist, meta = alg.run()
```

## Estruturas de Dados
- `strings`: lista de strings de entrada
- `alphabet`: string com todos os símbolos possíveis
- `pop_size`: tamanho da população

## Funções/Classes Principais
- `BLFGA`: Implementação principal do algoritmo
- `BLFGAAlgorithm`: Wrapper de integração ao framework
- `genetic_ops`: Operadores genéticos (crossover, mutação, etc)

## Referências
- Goldberg, D. E. (1989). Genetic Algorithms in Search, Optimization, and Machine Learning.
- Eiben, A. E., & Smith, J. E. (2015). Introduction to Evolutionary Computing.

---
*Recomendado para instâncias médias/grandes e quando se busca soluções de alta qualidade com flexibilidade adaptativa.*
