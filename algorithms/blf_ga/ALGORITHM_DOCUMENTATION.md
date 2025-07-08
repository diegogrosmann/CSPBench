# Documentação do Algoritmo BLF-GA

## Visão Geral

O BLF-GA (Blockwise Learning Fusion + Genetic Algorithm) é uma metaheurística híbrida desenvolvida para resolver o Closest String Problem (CSP). Ele combina três estratégias principais:

1. **Blockwise Learning**: Aprendizado local por blocos
2. **Genetic Algorithm**: Algoritmo genético para evolução global
3. **Fusion**: Mecanismos de fusão entre conhecimento local e global

## Arquitetura do Algoritmo

```
┌─────────────────────────────────────────────────────────────────┐
│                        BLF-GA ARCHITECTURE                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌─────────────────┐    ┌─────────────────┐    ┌─────────────┐  │
│  │   BLOCKWISE     │    │    GENETIC      │    │   FUSION    │  │
│  │   LEARNING      │◄──►│   ALGORITHM     │◄──►│ MECHANISMS  │  │
│  │                 │    │                 │    │             │  │
│  │ • Block Division│    │ • Population    │    │ • Adaptive  │  │
│  │ • Local Consensus│   │ • Selection     │    │ • Refinement│  │
│  │ • Pattern Repo  │    │ • Crossover     │    │ • Diversity │  │
│  │ • Entropy-based │    │ • Mutation      │    │ • Restart   │  │
│  └─────────────────┘    └─────────────────┘    └─────────────┘  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Fluxo de Execução

### 1. Inicialização
- **População Inicial**: Criada usando estratégia híbrida
  - Consenso global das strings de entrada
  - Variações inteligentes substituindo blocos
  - Strings aleatórias para diversidade
- **Divisão em Blocos**: Blocos uniformes iniciais
- **Parâmetros**: Configuração dos operadores e critérios

### 2. Loop Principal

```python
for geração in range(max_gerações):
    # 1. MECANISMOS ADAPTATIVOS
    if geração % immigrant_freq == 0:
        injetar_imigrantes_aleatórios()
    
    diversidade = calcular_diversidade_populacional()
    if diversidade < threshold:
        aumentar_mutação_temporariamente()
    
    # 2. APRENDIZADO POR BLOCOS
    repositório = aprender_padrões_por_bloco(população)
    
    # 3. EVOLUÇÃO GENÉTICA
    nova_população = []
    for i in range(tamanho_população):
        pai1 = seleção_torneio(população)
        pai2 = seleção_torneio(população)
        
        if random() < prob_crossover:
            filho1, filho2 = crossover(pai1, pai2)
        else:
            filho1, filho2 = pai1, pai2
        
        filho1 = mutação(filho1)
        filho2 = mutação(filho2)
        
        nova_população.extend([filho1, filho2])
    
    # 4. ELITISMO E REFINAMENTO
    elites = melhores_indivíduos(população)
    elites_refinados = refinamento_local(elites)
    nova_população[:len(elites)] = elites_refinados
    
    # 5. REDIVISÃO ADAPTATIVA
    if geração % rediv_freq == 0:
        blocos = redivisão_baseada_em_entropia(população)
    
    # 6. CRITÉRIOS DE PARADA
    if solução_ótima_encontrada() or timeout() or early_stopping():
        break
```

### 3. Finalização
- Retorna melhor solução encontrada
- Histórico de fitness por geração
- Estatísticas de execução

## Componentes Principais

### Blockwise Learning
- **Divisão em Blocos**: Strings são divididas em segmentos menores
- **Consenso Local**: Para cada bloco, calcula-se o padrão mais comum
- **Repositório de Padrões**: Armazena conhecimento aprendido
- **Adaptação Dinâmica**: Blocos são redefinidos baseado em entropia

### Genetic Algorithm
- **Seleção**: Torneio binário ou k-ário
- **Crossover**: One-point, uniform, ou blend-blocks
- **Mutação**: Multi-point, inversion, ou transposition
- **Elitismo**: Preserva os melhores indivíduos

### Fusion Mechanisms
- **Mutação Adaptativa**: Ajusta taxa baseada em diversidade
- **Imigrantes Aleatórios**: Injeta diversidade periodicamente
- **Refinamento Local**: Busca intensiva nos melhores
- **Restart**: Reinicia parcialmente se estagnado

## Parâmetros Importantes

### População
- `pop_size`: Tamanho da população (absoluto ou proporcional)
- `min_pop_size`: Tamanho mínimo garantido
- `elite_rate`: Porcentagem de elites preservados

### Blocos
- `initial_blocks`: Número inicial de blocos
- `min_block_len`: Tamanho mínimo de bloco
- `rediv_freq`: Frequência de redivisão

### Operadores
- `cross_prob`: Probabilidade de crossover
- `mut_prob`: Probabilidade de mutação
- `crossover_type`: Tipo de crossover
- `mutation_type`: Tipo de mutação

### Adaptação
- `immigrant_freq`: Frequência de imigrantes
- `diversity_threshold`: Limite para diversidade
- `mutation_adapt_factor`: Fator de adaptação da mutação

### Critérios
- `max_gens`: Gerações máximas
- `max_time`: Tempo máximo
- `no_improve_patience`: Paciência para early stopping

## Vantagens

1. **Aprendizado Local**: Captura padrões específicos em regiões das strings
2. **Evolução Global**: Mantém diversidade e exploração ampla
3. **Adaptabilidade**: Ajusta comportamento dinamicamente
4. **Robustez**: Múltiplos mecanismos de escape de ótimos locais
5. **Eficiência**: Paralelização de avaliações

## Limitações

1. **Complexidade**: Muitos parâmetros para ajustar
2. **Overhead**: Pode ser lento para instâncias muito pequenas
3. **Memória**: Requer armazenamento de população e estruturas auxiliares
4. **Sintonização**: Diferentes tipos de instâncias podem precisar configurações específicas

## Casos de Uso Recomendados

- **Instâncias Médias/Grandes**: n > 20, L > 100
- **Dados com Padrões Locais**: Quando existe estrutura regional
- **Qualidade vs Tempo**: Quando qualidade é prioritária sobre velocidade
- **Exploração Extensiva**: Quando é importante evitar ótimos locais

## Exemplos de Configuração

### Instância Pequena (n<10)
```python
config = {
    'pop_size': 50,
    'max_gens': 50,
    'initial_blocks': 0.3,
    'cross_prob': 0.8,
    'mut_prob': 0.15
}
```

### Instância Média (n=10-50)
```python
config = {
    'pop_size': 100,
    'max_gens': 100,
    'initial_blocks': 0.2,
    'cross_prob': 0.9,
    'mut_prob': 0.1
}
```

### Instância Grande (n>50)
```python
config = {
    'pop_size': 200,
    'max_gens': 200,
    'initial_blocks': 0.1,
    'cross_prob': 0.95,
    'mut_prob': 0.05
}
```

## Referências

O BLF-GA foi desenvolvido especificamente para o Closest String Problem, combinando técnicas de:
- Algoritmos Genéticos clássicos
- Aprendizado de máquina não supervisionado
- Otimização por enxame de partículas
- Busca local intensiva
