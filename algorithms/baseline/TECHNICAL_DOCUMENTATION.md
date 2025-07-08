# Baseline (Consenso Ganancioso) — Documentação Técnica

## Visão Geral

O algoritmo Baseline implementa uma estratégia de consenso ganancioso (greedy consensus) para o Closest String Problem (CSP). É uma solução simples, eficiente e determinística, servindo como referência para comparação com algoritmos mais sofisticados.

## Fundamentação Teórica

### Closest String Problem (CSP)
Dado um conjunto de strings S = {s₁, s₂, ..., sₙ} sobre um alfabeto Σ, o objetivo é encontrar uma string c que minimize a maior distância de Hamming para qualquer string do conjunto.

### Estratégia Greedy
O algoritmo constrói a string consenso posição por posição. Para cada posição, escolhe o símbolo do alfabeto que minimiza a maior distância de Hamming parcial até aquele ponto.

## Algoritmo

1. Para cada posição i (0 ≤ i < L):
   - Para cada símbolo do alfabeto:
     - Calcula a maior distância de Hamming parcial se esse símbolo for escolhido na posição i.
   - Seleciona o símbolo que minimiza essa distância máxima.
2. Repete até completar todas as posições.
3. Retorna a string consenso construída.

## Complexidade
- **Tempo:** O(L × |Σ| × n)
- **Espaço:** O(L)

## Pontos Fortes
- Extremamente rápido e simples
- Determinístico
- Serve como baseline para avaliação de heurísticas

## Limitações
- Não garante solução ótima global
- Pode ser sensível a outliers
- Não explora múltiplas soluções

## Exemplo de Uso
```python
from algorithms.baseline.algorithm import BaselineAlg
alg = BaselineAlg(strings, alphabet)
center, dist, meta = alg.run()
```

## Estruturas de Dados
- `strings`: lista de strings de entrada
- `alphabet`: string com todos os símbolos possíveis

## Funções Principais
- `greedy_consensus(strings, alphabet)`: Gera a string consenso
- `max_distance(center, strings)`: Calcula a maior distância de Hamming

## Referências
- Gusfield, D. (1997). Algorithms on Strings, Trees, and Sequences. Cambridge University Press.

---
*Este algoritmo é recomendado para instâncias pequenas ou como referência de comparação.*
