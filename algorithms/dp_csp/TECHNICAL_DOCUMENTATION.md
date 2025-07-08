# DP-CSP (Dynamic Programming Closest String Problem) — Documentação Técnica

## Visão Geral
O **DP-CSP** é um algoritmo exato baseado em programação dinâmica para o Closest String Problem (CSP). Ele resolve o problema de forma determinística, encontrando a menor distância máxima possível entre a string centro e todas as strings do conjunto.

## Arquitetura do Algoritmo

### 1. Decisão por Programação Dinâmica
- Para um valor de raio d, verifica se existe uma string centro c tal que max_distance(c, S) ≤ d
- Usa estados representando o número de erros restantes permitidos para cada string
- Avança posição por posição, testando todos os caracteres possíveis
- Reconstrói a solução se viável

### 2. Busca Binária Incremental
- Executa a DP para valores crescentes de d até encontrar o menor valor viável
- Garante a solução ótima global

## Fluxo Geral
1. Inicializa d (pode usar upper bound do baseline)
2. Para cada d crescente:
   - Executa _dp_decision para verificar viabilidade
   - Se viável, reconstrói a string centro
   - Se não, incrementa d
3. Retorna a menor string centro encontrada

## Complexidade
- **Tempo:** O(L × n × |Σ| × d^n) (exponencial no número de strings)
- **Espaço:** O(L × d^n)

## Pontos Fortes
- Garante solução ótima global
- Determinístico
- Útil para instâncias pequenas e médias

## Limitações
- Inviável para n > 8 ou L muito grande (explosão combinatória)
- Alto consumo de memória

## Exemplo de Uso
```python
from algorithms.dp_csp.algorithm import DPCSPAlgorithm
alg = DPCSPAlgorithm(strings, alphabet, max_d=5)
center, dist, meta = alg.run()
```

## Estruturas de Dados
- `strings`: lista de strings de entrada
- `alphabet`: string com todos os símbolos possíveis
- `max_d`: limite superior para a distância máxima

## Funções/Classes Principais
- `exact_dp_closest_string`: Busca incremental do menor raio possível
- `_dp_decision`: DP decisório para existência de centro

## Referências
- Ma, B., & Sun, F. (2009). More efficient algorithms for closest string and substring problems. SIAM J. Comput.
- Gusfield, D. (1997). Algorithms on Strings, Trees, and Sequences.

---
*Recomendado apenas para instâncias pequenas ou para validação de heurísticas.*
