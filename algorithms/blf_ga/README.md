# BLF-GA: Blockwise Learning Fusion + Genetic Algorithm

O BLF-GA é uma metaheurística híbrida que combina aprendizado por blocos e algoritmo genético para resolver o Closest String Problem.

## Estratégia

- Divide as strings em blocos (B-Splitter).
- Aplica aprendizado local em cada bloco (Block Learning).
- Usa evolução genética para otimização global (Genetic Evolution).
- Reporta progresso via callback.

## Heurísticas

- **Learning Fusion**: Combina conhecimento local e busca global.
- **Adaptive Blocking**: Redefine blocos dinamicamente.
- **Elite Refinement**: Busca local intensiva nos melhores indivíduos.

## Parâmetros Principais

- `pop_size`: Tamanho da população
- `initial_blocks`: Número inicial de blocos
- `min_block_len`: Tamanho mínimo de bloco
- `cross_prob`: Probabilidade de crossover
- `mut_prob`: Probabilidade de mutação
- `elite_rate`: Taxa de elitismo
- `rediv_freq`: Frequência de redivisão dos blocos
- `max_gens`: Número máximo de gerações
- `max_time`: Tempo máximo de execução (segundos)
- `seed`: Semente aleatória

Veja `config.py` para todos os parâmetros e valores padrão.

## Uso Ideal

- Instâncias médias a grandes (n > 20, L > 100).
- Dados com padrões locais.
- Quando qualidade próxima do ótimo é prioritária.

## Limitações

- Muitos parâmetros.
- Pode ser lento para instâncias pequenas.
- Pode estagnar em ótimos locais.

## Exemplo de Uso

```python
from algorithms.blf_ga.algorithm import BLFGAAlgorithm
alg = BLFGAAlgorithm(strings, alphabet, pop_size=100, max_gens=50)
center, dist = alg.run()
```

## Documentação

- Consulte o código para docstrings detalhadas (Google style).
- Integração automática com o framework CSP via decorador `@register_algorithm`.
