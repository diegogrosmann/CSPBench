# Baseline Algorithm

O algoritmo Baseline implementa uma solução simples e eficiente para o Closest String Problem usando consenso ganancioso. Serve como referência para comparação com métodos mais sofisticados.

## Estratégia

- Para cada posição, escolhe o símbolo mais frequente entre todas as strings.
- Constrói a string consenso e calcula a distância máxima de Hamming em relação às entradas.

## Características

- **Determinístico** e **instantâneo**.
- Complexidade linear: O(n × L × |Σ|).
- Não requer callback de progresso.
- Documentação detalhada no estilo Google disponível no código-fonte.

## Parâmetros

- Não possui parâmetros ajustáveis além do dataset.

## Uso

Ideal como baseline para benchmarking e validação de instâncias simples ou bem estruturadas.

### Exemplo de Uso

```python
from algorithms.baseline.algorithm import BaselineAlg
alg = BaselineAlg(strings, alphabet)
center, dist = alg.run()
```

## Limitações

- Não considera dependências entre posições.
- Pode não encontrar o ótimo global em casos de empate ou alta diversidade.

## Cenários de Uso Ideal
- Benchmark para comparação com algoritmos mais complexos
- Solução rápida quando tempo é limitado
- Pré-processamento para outros algoritmos
- Dados com forte consenso por posição

## Documentação

- Consulte o código para docstrings detalhadas (Google style).
- Integração automática com o framework CSP via decorador `@register_algorithm`.
