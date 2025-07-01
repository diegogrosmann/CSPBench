# DP-CSP: Dynamic Programming for Closest String Problem

O DP-CSP é um algoritmo exato baseado em programação dinâmica para encontrar a solução ótima do Closest String Problem.

## Estratégia

- Busca incremental pelo menor raio d viável.
- Mantém estados representando orçamentos de erro por string.
- Reconstrói a solução ótima via tabela de transições.

## Garantias

- Sempre encontra o ótimo global.
- Determinístico e completo.

## Limitações

- Complexidade exponencial em n (número de strings).
- Prático apenas para n ≤ 10 e d pequeno.
- Uso intensivo de memória.

## Uso Ideal

- Instâncias pequenas.
- Validação de heurísticas.
- Benchmark de qualidade.

## Parâmetros Principais

- `max_d`: Limite superior para raio
- `warn_threshold`: Threshold para alertas de complexidade
- `max_time`: Timeout em segundos

Veja `config.py` para todos os parâmetros e valores padrão.

## Exemplo de Uso

```python
from algorithms.dp_csp.algorithm import DPCSPAlgorithm
alg = DPCSPAlgorithm(strings, alphabet, max_d=5)
center, dist = alg.run()
```

## Documentação

- Consulte o código para docstrings detalhadas (Google style).
- Integração automática com o framework CSP via decorador `@register_algorithm`.

## Características e Limitações

### Garantias Teóricas
- **Optimalidade**: Sempre encontra a solução de raio mínimo
- **Completude**: Se existe solução ≤ max_d, será encontrada
- **Determinismo**: Sempre produz o mesmo resultado

### Complexidade Exponencial
**Escalabilidade**: Viável apenas para instâncias pequenas
- **Prático**: n ≤ 10, d ≤ 8
- **Limite**: (d+1)ⁿ cresce exponencialmente
- **Memória**: Pode esgotar RAM rapidamente

### Cenários de Uso

#### Ideal Para:
- **Instâncias Pequenas**: n ≤ 10 strings
- **Raio Pequeno**: d ≤ 8 esperado
- **Verificação de Optimalidade**: Validar outros algoritmos
- **Benchmark**: Comparação com heurísticas

#### Inadequado Para:
- **Instâncias Grandes**: n > 15
- **Raio Grande**: d > 10
- **Tempo Limitado**: Pode ser muito lento
- **Memória Limitada**: Uso exponencial de RAM

### Otimizações Possíveis

#### Técnicas Avançadas
- **Branching**: Escolha inteligente de símbolos por posição
- **Memoização**: Cache de subproblemas recorrentes
- **Poda por Cotas**: Eliminate estados claramente subótimos
- **Paralelização**: Estados independentes podem ser processados em paralelo

#### Variações do Algoritmo
- **Beam DP**: Manter apenas k melhores estados por nível
- **A* Search**: Usar heurística para guiar busca
- **Iterative Deepening**: Combinar busca incremental com limites de memória

## Performance Típica

### Instâncias Favoráveis
- **Strings Similares**: Quando raio ótimo é pequeno
- **Alfabeto Pequeno**: DNA (4 símbolos) vs. proteínas (20 símbolos)
- **Dados Estruturados**: Padrões que facilitam poda

### Tempo de Execução
- **n=5, d=3**: Milissegundos
- **n=8, d=5**: Segundos
- **n=10, d=7**: Minutos
- **n=12, d=8**: Horas ou inviável

### Uso de Memória
- **Estados Simultâneos**: Até (d+1)ⁿ estados por nível
- **Parent Table**: L × estados_visitados
- **Total**: Pode atingir gigabytes para n ≥ 12

O DP-CSP é fundamental como **padrão-ouro** para validação de outros algoritmos, garantindo que sempre encontramos a resposta correta quando computacionalmente viável, servindo como referência absoluta para problemas pequenos.
