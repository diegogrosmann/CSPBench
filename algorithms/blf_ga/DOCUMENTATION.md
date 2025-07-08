# Documentação do Algoritmo BLF-GA

## Resumo Executivo

O **BLF-GA** (Blockwise Learning Fusion + Genetic Algorithm) é uma metaheurística híbrida desenvolvida para resolver o **Closest String Problem (CSP)**. O algoritmo combina três estratégias principais:

1. **Blockwise Learning** - Aprendizado local por divisão em blocos
2. **Genetic Algorithm** - Evolução global através de operadores genéticos
3. **Fusion** - Fusão inteligente de conhecimento local e global

## Arquitetura do Algoritmo

### Visão Geral
```
BLF-GA
├── Inicialização Inteligente
│   ├── Consenso das strings de entrada
│   ├── Variações baseadas em blocos
│   └── Diversidade aleatória
├── Loop Principal de Evolução
│   ├── Mecanismos Adaptativos
│   │   ├── Imigrantes aleatórios
│   │   ├── Mutação adaptativa
│   │   └── Controle de diversidade
│   ├── Aprendizado por Blocos
│   │   ├── Extração de padrões locais
│   │   └── Criação de repositório
│   ├── Evolução Genética
│   │   ├── Seleção por torneio
│   │   ├── Crossover multi-tipo
│   │   ├── Mutação multi-tipo
│   │   └── Elitismo adaptativo
│   ├── Refinamento Local
│   │   ├── Busca local nos elites
│   │   └── Melhoria de qualidade
│   └── Redivisão Adaptativa
│       ├── Análise de entropia
│       └── Ajuste da estrutura de blocos
└── Critérios de Parada
    ├── Solução ótima encontrada
    ├── Limite de tempo/gerações
    ├── Early stopping
    └── Restart automático
```

### Componentes Principais

#### 1. Blockwise Learning (Aprendizado por Blocos)
- **Objetivo**: Capturar padrões locais eficazes
- **Método**: Divide strings em blocos e aprende consenso
- **Vantagens**: 
  - Reduz complexidade do espaço de busca
  - Explora estrutura local dos dados
  - Permite adaptação dinâmica

#### 2. Genetic Algorithm (Algoritmo Genético)
- **Operadores**: Seleção, crossover, mutação, elitismo
- **Tipos de Crossover**: one-point, uniform, blend-blocks
- **Tipos de Mutação**: multi, inversion, transposition
- **Seleção**: Torneio com tamanho configurável

#### 3. Fusion (Fusão)
- **Mecanismos Adaptativos**: 
  - Mutação baseada em diversidade/convergência
  - Imigrantes aleatórios para diversidade
  - Refinamento local intensivo
- **Controles Dinâmicos**:
  - Redivisão de blocos por entropia
  - Elitismo adaptativo
  - Restart automático

## Parâmetros Principais

### População
- `pop_size`: Tamanho (int fixo ou float para multiplicador)
- `min_pop_size`: Tamanho mínimo garantido (20)
- `elite_rate`: Taxa de elitismo (0.05 = 5%)

### Blocos
- `initial_blocks`: Número inicial (int fixo ou float para proporção)
- `min_block_len`: Tamanho mínimo (1)
- `rediv_freq`: Frequência de redivisão (10 gerações)

### Operadores Genéticos
- `cross_prob`: Probabilidade de crossover (0.9)
- `mut_prob`: Probabilidade de mutação (0.1)
- `tournament_k`: Tamanho do torneio (2)

### Mecanismos Adaptativos
- `immigrant_freq`: Frequência de imigrantes (10 gerações)
- `immigrant_ratio`: Proporção de imigrantes (0.2)
- `diversity_threshold`: Limiar de diversidade (0.4)
- `mutation_adapt_factor`: Fator de adaptação da mutação (2.0)

### Critérios de Parada
- `max_gens`: Máximo de gerações (100)
- `max_time`: Tempo máximo em segundos (1200)
- `no_improve_patience`: Early stopping (0.2 = 20% de max_gens)

## Recomendações de Uso

### Por Tamanho de Instância

| Tamanho | n (strings) | L (comprimento) | pop_size | max_gens | Observações |
|---------|-------------|-----------------|----------|----------|-------------|
| Pequeno | < 10 | < 50 | 50-100 | 50-100 | Rápido, boa qualidade |
| Médio | 10-50 | 50-200 | 100-200 | 100-200 | Balanceado |
| Grande | > 50 | > 200 | 200-500 | 200-500 | Foco em qualidade |

### Por Características dos Dados

| Característica | Recomendação | Justificativa |
|----------------|--------------|---------------|
| Dados com padrões locais | `initial_blocks=0.1-0.3` | Melhor captura de padrões |
| Dados aleatórios | `immigrant_ratio=0.3-0.5` | Maior diversidade necessária |
| Alfabeto pequeno | `mut_prob=0.05-0.1` | Menos mutação necessária |
| Alfabeto grande | `mut_prob=0.1-0.2` | Mais exploração necessária |

## Exemplo de Uso

```python
from algorithms.blf_ga.implementation import BLFGA

# Dados de exemplo
strings = ['ATCGATCG', 'ATCGATGG', 'ATCGATCC']
alphabet = 'ATCG'

# Configuração padrão
blfga = BLFGA(strings, alphabet)
best_string, best_fitness, history = blfga.run()

# Configuração personalizada
blfga_custom = BLFGA(
    strings=strings,
    alphabet=alphabet,
    pop_size=100,          # População fixa de 100
    max_gens=50,           # Máximo 50 gerações
    cross_prob=0.95,       # Alta recombinação
    mut_prob=0.05,         # Baixa mutação
    elite_rate=0.1,        # 10% de elites
    immigrant_ratio=0.3,   # 30% de imigrantes
    seed=42                # Reprodutibilidade
)
```

## Vantagens do BLF-GA

1. **Híbrido Inteligente**: Combina busca local e global eficientemente
2. **Adaptativo**: Ajusta parâmetros dinamicamente durante execução
3. **Robusto**: Múltiplos mecanismos anti-estagnação
4. **Flexível**: Muitos parâmetros configuráveis
5. **Eficiente**: Paralelização e otimizações de performance

## Limitações

1. **Complexidade**: Muitos parâmetros para ajustar
2. **Overhead**: Pode ser lento para instâncias muito pequenas
3. **Memória**: Uso moderado a alto de memória
4. **Tuning**: Requer ajuste fino para dados específicos

## Status da Implementação

- ✅ Algoritmo principal implementado
- ✅ Todos os mecanismos adaptativos funcionais
- ✅ Paralelização básica implementada
- ✅ Documentação completa
- ✅ Garantia de população mínima
- ✅ Testes de validação passando
- ✅ Integração com framework CSP

## Contribuições Principais

1. **Garantia de População Mínima**: Evita populações muito pequenas em instâncias menores
2. **Documentação Completa**: Código totalmente documentado com explicações didáticas
3. **Flexibilidade de Configuração**: Suporte a valores fixos e proporcionais
4. **Mecanismos Adaptativos**: Controles dinâmicos para melhor performance
