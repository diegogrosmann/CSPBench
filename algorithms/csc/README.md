# CSC: Consensus String Clustering

## Visão Geral
O algoritmo CSC (Consensus String Clustering) resolve o Closest String Problem através de uma abordagem em três etapas: **clusterização por similaridade**, **consenso local por cluster** e **recombinação inteligente de blocos**. A estratégia central é identificar grupos naturais nas strings e combinar suas características mais promissoras.

## Funcionamento

### Estratégia Principal
O CSC baseia-se na hipótese de que strings similares devem ser agrupadas e processadas juntas, permitindo descoberta de padrões locais que podem ser recombinados para formar uma solução global superior.

### Fluxo de Execução Detalhado

#### 1. Pré-processamento e Análise
- **Cálculo de Distâncias**: Computa matriz de distâncias de Hamming entre todos os pares
- **Estatísticas Descritivas**: Analisa média, mediana, mínimo e máximo das distâncias
- **Parâmetros Automáticos**: Define automaticamente raio de clusterização e número de blocos

#### 2. Clusterização DBSCAN
```
Para cada string s:
  ├── Identificar vizinhos dentro do raio d
  ├── Formar cluster se min_samples satisfeito
  ├── Expandir cluster recursivamente
  └── Marcar outliers como ruído
```

#### 3. Consenso Local por Cluster
- **Análise Posicional**: Para cada posição em cada cluster
- **Votação Majoritária**: Escolhe símbolo mais frequente por posição
- **Resolução de Empates**: Critério determinístico para desempate
- **Validação Local**: Verifica qualidade do consenso gerado

#### 4. Recombinação de Blocos
- **Divisão em Blocos**: Segmenta cada consenso em n_blocks blocos
- **Exploração Combinatorial**: Testa todas as combinações possíveis de blocos
- **Avaliação de Candidatos**: Calcula distância máxima para cada combinação
- **Seleção do Melhor**: Escolhe candidato com menor raio

#### 5. Busca Local (Hill-Climbing)
- **Refinamento Posicional**: Tenta melhorar cada posição individualmente
- **Exploração de Vizinhança**: Testa símbolos alternativos
- **Critério de Melhoria**: Aceita mudanças que reduzem distância máxima
- **Convergência**: Para quando nenhuma melhoria é possível

## Heurísticas Principais

### Clusterização Baseada em Densidade
**Princípio**: Strings próximas compartilham características estruturais
- **DBSCAN**: Identifica clusters de densidade arbitrária
- **Tolerância a Ruído**: Outliers não formam clusters próprios
- **Flexibilidade**: Não assume número fixo de clusters

### Consenso Hierárquico
**Estratégia**: Consenso local antes de combinação global
- **Redução de Ruído**: Cada cluster filtra suas variações internas
- **Preservação de Padrões**: Mantém características dominantes de cada grupo
- **Robustez**: Reduz impacto de strings atípicas

### Recombinação Inteligente
**Abordagem**: Combina melhores partes de diferentes soluções
- **Granularidade Controlada**: n_blocks define resolução da recombinação
- **Exploração Exaustiva**: Testa todas as combinações viáveis
- **Avaliação Rigorosa**: Cada candidato é testado contra todas as strings

### Busca Local Intensiva
**Refinamento**: Melhoria iterativa da melhor solução encontrada
- **Hill-Climbing Simples**: Mudanças que sempre melhoram
- **Busca Posicional**: Explora alternativas para cada posição
- **Convergência Local**: Garante ótimo local na vizinhança

## Parâmetros de Configuração

### Clusterização Automática
- **d (raio)**: Calculado como 80% da distância média (padrão: auto)
- **min_samples**: Mínimo 2 strings por cluster
- **distance_metric**: Hamming distance customizada

### Recombinação
- **n_blocks**: Número de blocos para divisão (padrão: auto, baseado em log(n))
- **max_combinations**: Limitado pelo produto cartesiano dos blocos

### Critérios Automáticos
```
d = max(2, floor(média_distâncias × 0.8))
n_blocks = max(2, min(4, n//6, L//25))
```

## Características Avançadas

### Parâmetros Adaptativos
- **Auto-tuning**: Parâmetros ajustados baseados nas características dos dados
- **Análise Estatística**: Usa distribuição de distâncias para calibração
- **Robustez**: Funciona bem mesmo com parâmetros subótimos

### Tratamento de Casos Especiais
- **Clusters Únicos**: Quando todas as strings formam um cluster
- **Sem Clusters**: Fallback para consenso global
- **Strings Idênticas**: Detecção e tratamento eficiente

### Validação Integrada
- **Verificação de Qualidade**: Cada etapa valida seus resultados
- **Fallback Strategies**: Alternativas quando etapas falham
- **Logging Detalhado**: Rastreamento completo do processo

## Performance e Escalabilidade

### Complexidade Temporal
- **Clusterização**: O(n² × L) para cálculo de distâncias + O(n log n) para DBSCAN
- **Consenso**: O(k × L × |Σ|) onde k é número de clusters
- **Recombinação**: O(B^k × n × L) onde B é blocos por cluster
- **Busca Local**: O(L × |Σ| × n) no pior caso

### Complexidade Espacial
- **Matriz de Distâncias**: O(n²) se materializadas (pode ser otimizada)
- **Clusters**: O(n) para armazenamento
- **Candidatos**: O(B^k × L) para todas as combinações

### Cenários de Uso Ideal
- **Dados com Estrutura Natural**: Quando strings formam grupos óbvios
- **Instâncias Médias**: n = 20-200, onde clusterização é efetiva
- **Alfabetos Pequenos**: DNA/RNA onde recombinação é mais promissora
- **Qualidade vs. Tempo**: Quando se busca boa solução em tempo moderado

### Vantagens
- **Parâmetros Automáticos**: Pouca necessidade de ajuste manual
- **Robustez a Outliers**: DBSCAN trata naturalmente ruído
- **Exploração Sistemática**: Recombinação garante boa cobertura do espaço
- **Interpretabilidade**: Clusters oferecem insights sobre os dados

### Limitações
- **Complexidade Combinatorial**: Recombinação pode explodir com muitos clusters
- **Dependência de Clusters**: Performance degrada se não há estrutura clara
- **Memória**: Pode usar muita memória para instâncias grandes
- **Parâmetros Sensíveis**: Embora automáticos, podem não ser ótimos para todos os casos

O CSC é especialmente efetivo para dados que possuem estrutura natural de agrupamento, oferecendo uma abordagem sistemática que combina análise de dados com otimização combinatorial.
