# H³-CSP: Hybrid Hierarchical Hamming Search

O H³-CSP é um algoritmo híbrido de três camadas que combina decomposição hierárquica, técnicas especializadas por bloco e refinamento global.

## Arquitetura

1. **B-Splitter**: Divide as strings em blocos contíguos (~√L).
2. **Smart-Core**: Seleciona técnica ótima por bloco (busca exaustiva ou beam search).
3. **Global Refine**: Combina blocos e aplica hill-climbing global.

## Heurísticas

- Decomposição adaptativa (√L Rule).
- Seleção inteligente de técnica por bloco.
- Geração de candidatos diversos.
- Refinamento global por busca local.

## Parâmetros

- Número de blocos, limites para busca exaustiva, largura do beam, candidatos por bloco, iterações de refinamento.

## Uso Ideal

- Instâncias de tamanho médio (L entre 50-500).
- Dados com padrões locais.
- Quando se busca equilíbrio entre qualidade e eficiência.

## Limitações

- Mais complexo que algoritmos simples.
- Overhead para instâncias pequenas.
- Vários parâmetros para ajuste fino.
B = ⌈√L⌉  // número de blocos
block_size = ⌈L/B⌉

Para i = 0 até B-1:
  ├── início = i × block_size
  ├── fim = min((i+1) × block_size, L)
  └── bloco[i] = substring[início:fim]
```

#### Fase 2: Processamento Inteligente por Bloco
```
Para cada bloco b:
  ├── Analisa características (tamanho, diversidade)
  ├── IF tamanho_pequeno AND |Σ|^tamanho < limite:
  │     Usa busca_exaustiva(b)
  ├── ELSE:
  │     Usa beam_search(b, beam_width)
  ├── Gera k melhores candidatos para bloco b
  └── Adiciona candidatos ao repositório
```

#### Fase 3: Fusão e Refinamento Global
```
candidatos_globais = []
Para cada combinação de blocos:
  ├── Concatena blocos para formar candidato completo
  ├── Avalia qualidade (distância máxima)
  └── Adiciona aos candidatos globais

melhor = min(candidatos_globais, key=qualidade)
melhor_refinado = hill_climbing(melhor)
```

## Heurísticas Principais

### Decomposição Adaptativa
**Princípio**: Tamanho de bloco balanceia precisão vs. eficiência
- **√L Rule**: Blocos de tamanho √L maximizam trade-off teórico
- **Contiguidade**: Blocos adjacentes preservam correlações posicionais
- **Flexibilidade**: Último bloco pode ter tamanho diferente

### Seleção Inteligente de Técnica
**Estratégia**: Técnica ótima varia com características do bloco

#### Busca Exaustiva
- **Quando**: |Σ|^m < 10,000 (m = tamanho do bloco)
- **Vantagem**: Encontra ótimo local garantido
- **Custo**: Exponencial, mas controlado pelo limite

#### Beam Search
- **Quando**: Blocos grandes ou alfabeto grande
- **Estratégia**: Mantém beam_width melhores candidatos por nível
- **Balanceamento**: Explora largura suficiente sem explodir computação

### Geração de Candidatos Diversificados
**Objetivo**: k candidatos diversos e promissores por bloco
- **Quality-Based**: Ordena por distância de Hamming
- **Diversity Injection**: Inclui consenso local sempre
- **Fallback Strategy**: Se busca falha, usa candidatos existentes

### Hill-Climbing Global
**Refinement**: Melhoria iterativa da melhor solução
- **Position-wise**: Testa cada posição independentemente
- **Alphabet Scanning**: Para cada posição, testa todos os símbolos
- **Greedy Acceptance**: Aceita mudanças que melhoram fitness
- **Local Optimum**: Para quando nenhuma melhoria é possível

## Parâmetros de Configuração

### Decomposição
- **auto_blocks**: Usa √L como número de blocos (padrão: true)
- **min_block_size**: Tamanho mínimo de bloco (padrão: 2)
- **max_blocks**: Limite superior de blocos (padrão: L/2)

### Técnicas por Bloco
- **exhaustive_limit**: Limite para busca exaustiva (padrão: 10,000)
- **beam_width**: Largura do beam search (padrão: 32)
- **k_candidates**: Candidatos por bloco (padrão: 5)

### Refinamento Global
- **local_search_iters**: Iterações de hill-climbing (padrão: 3)
- **max_time**: Tempo limite total (padrão: 300s)

### Controle de Qualidade
- **diversity_threshold**: Mínima diversidade entre candidatos
- **fallback_enabled**: Usar estratégias de fallback (padrão: true)

## Características Avançadas

### Análise Adaptativa de Blocos
**Métricas de Decisão**:
- **Diversidade**: Entropia posicional dentro do bloco
- **Tamanho**: Número de posições no bloco
- **Alfabeto Efetivo**: Símbolos únicos observados no bloco
- **Complexidade**: |Σ|^m estimada

### Estratégias de Fallback
**Robustez**: Garantias mínimas mesmo quando técnicas principais falham
- **Consensus Fallback**: Se busca exaustiva falha, usa consenso
- **Random Sampling**: Se beam search falha, amostra aleatoriamente
- **Existing Blocks**: Reutiliza blocos de strings originais

### Otimizações de Performance
- **Early Termination**: Para beam search quando convergir
- **Candidate Caching**: Evita recomputações de blocos
- **Memory Management**: Controla uso de memória em instâncias grandes

## Performance e Escalabilidade

### Complexidade Temporal
**Por Bloco**:
- Exaustiva: O(|Σ|^m × n) onde m é tamanho do bloco
- Beam: O(m × beam_width × |Σ| × n)

**Total**: O(B × max(busca_por_bloco) + fusão + refinamento)

### Complexidade Espacial
- **Candidatos**: O(B × k × tamanho_médio_bloco)
- **Beam**: O(beam_width × m) temporariamente
- **Total**: Linear no tamanho do problema

### Cenários de Uso Ideal

#### Instâncias Favoráveis
- **Tamanho Médio**: L entre 50-500 (onde √L é balanceado)
- **Estrutura Local**: Padrões conservados em regiões específicas
- **Alfabeto Moderado**: |Σ| entre 4-20
- **Tempo Disponível**: Quando se pode aguardar processamento mais sofisticado

#### Vantagens Competitivas
- **Adaptabilidade**: Escolhe técnica ótima automaticamente
- **Escalabilidade**: Cresce suavemente com tamanho do problema
- **Robustez**: Múltiplas estratégias de fallback
- **Qualidade**: Combina precisão local com otimização global

### Limitações
- **Complexidade de Implementação**: Mais complexo que algoritmos simples
- **Overhead**: Pode ser excessivo para instâncias muito pequenas
- **Parâmetros**: Vários parâmetros para potencial ajuste
- **Dependência de Estrutura**: Performance varia com características dos dados

### Comparação com Outras Abordagens
- **vs. Baseline**: Melhor qualidade, maior custo computacional
- **vs. BLF-GA**: Menor tempo, pode ter qualidade ligeiramente inferior
- **vs. CSC**: Mais sistemático, menos dependente de estrutura de clusters
- **vs. DP-CSP**: Escalável, mas sem garantia de otimalidade

O H³-CSP representa uma abordagem equilibrada que combina a sistematicidade de métodos exatos com a praticidade de heurísticas, sendo especialmente efetivo para instâncias de tamanho médio onde tanto qualidade quanto eficiência são importantes.
