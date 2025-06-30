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

## Parâmetros

- Raio de clusterização, número de blocos, critérios automáticos.

## Uso Ideal

- Dados com estrutura natural de agrupamento.
- Instâncias médias (n = 20-200).
- Alfabetos pequenos (DNA/RNA).

## Limitações

- Complexidade combinatorial pode crescer rápido.
- Performance depende da existência de clusters naturais.
