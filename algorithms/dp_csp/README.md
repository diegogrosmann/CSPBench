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

## Parâmetros

- Limite superior para raio.
- Threshold para alertas de complexidade.
### Fluxo de Execução Detalhado

#### 1. Pré-computação
- **Matriz de Discrepâncias**: Para cada posição pos e símbolo σ:
  ```
  δ[pos][σ] = (δ₁, δ₂, ..., δₙ)
  onde δᵢ = 1 se strings[i][pos] ≠ σ, senão 0
  ```
- **Estimativa de Complexidade**: Calcula (d+1)ⁿ para alertar sobre uso de memória

#### 2. Inicialização do DP
- **Estado Inicial**: `(d, d, ..., d)` - orçamento completo para todas as strings
- **Fronteira**: Conjunto de estados alcançáveis
- **Tabela Parent**: Rastreia caminho para reconstrução da solução

#### 3. Propagação por Posições
```
Para cada posição pos = 0 até L-1:
  Para cada estado (r₁, ..., rₙ) na fronteira atual:
    Para cada símbolo σ no alfabeto:
      ├── Calcula novo_estado = (r₁-δ₁, ..., rₙ-δₙ)
      ├── Se algum rᵢ < 0: descarta (inviável)
      ├── Senão: adiciona à nova fronteira
      └── Registra transição na tabela parent
```

#### 4. Verificação de Viabilidade
- **Fronteira Vazia**: Se fronteira fica vazia, não há solução para este d
- **Fronteira Não-Vazia**: Existe pelo menos uma solução

#### 5. Reconstrução da Solução
```
estado_final = qualquer estado na fronteira final
pos = L
Enquanto pos > 0:
  ├── Consulta parent[pos][estado_final]
  ├── Recupera (estado_anterior, símbolo_usado)
  ├── Adiciona símbolo_usado ao resultado
  └── pos--, estado_final = estado_anterior
```

## Heurísticas e Optimizações

### Busca Incremental por Raio
**Princípio**: Testa valores crescentes de d até encontrar solução
- **Guarantee**: Primeiro d viável é necessariamente ótimo
- **Eficiência**: Evita testar raios desnecessariamente grandes
- **Early Stopping**: Para assim que encontra primeira solução

### Poda de Estados Inviáveis
**Critério**: Estados com orçamento negativo são descartados
- **Válido**: `min(r₁, ..., rₙ) ≥ 0`
- **Inválido**: Qualquer rᵢ < 0 indica impossibilidade
- **Benefit**: Reduz drasticamente o espaço de estados

### Representação Eficiente de Estados
**Estrutura**: Tuplas de inteiros como chaves de dicionário
- **Hashable**: Permite uso como chaves em dicionários Python
- **Compact**: Representação mínima de memória
- **Fast Access**: Operações O(1) para consulta/inserção

### Reconstrução Otimizada
**Parent Table**: Armazena apenas transições necessárias
- **Memória**: O(L × |Estados|) em vez de armazenar soluções completas
- **Reconstrução**: O(L) para recuperar solução
- **Flexibilidade**: Permite múltiplas soluções ótimas

## Parâmetros de Configuração

### Controle de Execução
- **max_d**: Limite superior para raio (padrão: distância baseline)
- **warn_threshold**: Alerta quando (d+1)ⁿ > 10^threshold

### Estimativas de Complexidade
- **Estados Potenciais**: (d+1)ⁿ
- **Memória**: O(L × (d+1)ⁿ) no pior caso
- **Tempo**: O(L × |Σ| × (d+1)ⁿ)

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
