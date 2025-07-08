# CSC (Consensus String Clustering) — Documentação Técnica

## Visão Geral
O **CSC** (Consensus String Clustering) é um algoritmo heurístico para o Closest String Problem (CSP) que combina técnicas de consenso, divisão em blocos e agrupamento (clustering) para encontrar soluções aproximadas de alta qualidade.

## Arquitetura do Algoritmo

### 1. Consenso por Maioria
- Calcula a string consenso para cada posição (baseline)
- Serve como ponto de partida para a heurística

### 2. Divisão em Blocos
- Divide as strings em n blocos aproximadamente iguais
- Permite recombinação de blocos de diferentes consensos

### 3. Clustering (DBSCAN)
- Converte strings em arrays numéricos
- Aplica DBSCAN para agrupar strings similares
- Cada cluster gera um consenso local

### 4. Recombinação e Busca
- Gera candidatos recombinando blocos de diferentes consensos
- Avalia todos os candidatos e seleciona o de menor distância máxima

## Fluxo Geral
1. Conversão das strings para arrays numéricos
2. Cálculo automático de parâmetros (raio d, n_blocks)
3. Clustering das strings (DBSCAN)
4. Cálculo de consensos locais por cluster
5. Geração de candidatos por recombinação de blocos
6. Seleção do melhor candidato

## Complexidade
- **Tempo:** O(n² × L) para clustering + O(B^n × L) para recombinação (B = n_blocks)
- **Espaço:** O(n × L)

## Pontos Fortes
- Explora estrutura local e global dos dados
- Pode encontrar soluções melhores que o consenso simples
- Determinístico

## Limitações
- Sensível à escolha dos parâmetros (d, n_blocks)
- Overhead para instâncias pequenas
- Não garante ótimo global

## Exemplo de Uso
```python
from algorithms.csc.algorithm import CSCAlgorithm
alg = CSCAlgorithm(strings, alphabet, d=3, n_blocks=4)
center, dist, meta = alg.run()
```

## Estruturas de Dados
- `strings`: lista de strings de entrada
- `alphabet`: string com todos os símbolos possíveis
- `d`: raio para clustering (DBSCAN)
- `n_blocks`: número de blocos para recombinação

## Funções/Classes Principais
- `heuristic_closest_string`: Heurística principal do CSC
- `consensus_string`: Consenso por maioria
- `split_blocks`, `recombine_blocks`: Manipulação de blocos

## Referências
- Ester, M. et al. (1996). A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise (DBSCAN).
- Gusfield, D. (1997). Algorithms on Strings, Trees, and Sequences.

---
*Recomendado para instâncias médias e quando se deseja explorar padrões locais e globais simultaneamente.*
