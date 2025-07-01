# CSC: Consensus String Clustering

O CSC resolve o Closest String Problem via clusterização, consenso local e recombinação inteligente de blocos.

## Estratégia

1. **Clusterização**: Agrupa strings similares (DBSCAN).
2. **Consenso Local**: Gera consenso por cluster.
3. **Recombinação**: Combina blocos dos consensos.
4. **Busca Local**: Refina a melhor solução.

## Heurísticas

- Parâmetros automáticos baseados nos dados.
- Robustez a outliers e ruído.
- Exploração combinatorial controlada.

## Parâmetros Principais

- `min_d`, `d_factor`, `min_blocks`, `max_blocks`, `n_div`, `l_div` (veja `config.py`)

## Uso Ideal

- Dados com estrutura natural de agrupamento.
- Instâncias médias (n = 20-200).
- Alfabetos pequenos (DNA/RNA).

## Limitações

- Complexidade combinatorial pode crescer rápido.
- Performance depende da existência de clusters naturais.

## Exemplo de Uso

```python
from algorithms.csc.algorithm import CSCAlgorithm
alg = CSCAlgorithm(strings, alphabet, min_blocks=3)
center, dist = alg.run()
```

## Documentação

- Consulte o código para docstrings detalhadas (Google style).
- Integração automática com o framework CSP via decorador `@register_algorithm`.
