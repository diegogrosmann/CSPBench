# H³-CSP: Hybrid Hierarchical Hamming Search

O H³-CSP é um algoritmo híbrido de três camadas que combina decomposição hierárquica, técnicas especializadas por bloco e refinamento global para resolver o Closest String Problem.

## Arquitetura

1. **B-Splitter**: Divide as strings em blocos contíguos (~√L).
2. **Smart-Core**: Seleciona técnica ótima por bloco (busca exaustiva ou beam search).
3. **Global Refine**: Combina blocos e aplica hill-climbing global.

## Heurísticas

- Decomposição adaptativa (√L Rule).
- Seleção inteligente de técnica por bloco.
- Geração de candidatos diversos.
- Refinamento global por busca local.

## Parâmetros Principais

- `auto_blocks`, `min_block_size`, `max_blocks`, `exhaustive_limit`, `beam_width`, `k_candidates`, `local_search_iters`, `max_time`, `diversity_threshold`, `seed` (veja `config.py`)

## Uso Ideal

- Instâncias de tamanho médio (L entre 50-500).
- Dados com padrões locais.
- Quando se busca equilíbrio entre qualidade e eficiência.

## Limitações

- Mais complexo que algoritmos simples.
- Overhead para instâncias pequenas.
- Vários parâmetros para ajuste fino.

## Exemplo de Uso

```python
from algorithms.h3_csp.algorithm import H3CSPAlgorithm
alg = H3CSPAlgorithm(strings, alphabet, beam_width=16)
center, dist = alg.run()
```

## Documentação

- Consulte o código para docstrings detalhadas (Google style).
- Integração automática com o framework CSP via decorador `@register_algorithm`.
