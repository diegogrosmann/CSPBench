# Baseline Algorithm

## Visão Geral
O algoritmo Baseline implementa uma solução simples e eficiente para o Closest String Problem usando consenso ganancioso. Serve como ponto de referência para comparar a performance de outros algoritmos mais sofisticados.

## Funcionamento

### Estratégia Principal
O algoritmo utiliza **consenso por maioria**: para cada posição, escolhe o símbolo que aparece com maior frequência naquela posição entre todas as strings de entrada.

### Fluxo de Execução

1. **Análise Posicional**
   - Para cada posição i (0 ≤ i < L):
   - Conta a frequência de cada símbolo nessa posição
   - Seleciona o símbolo mais frequente

2. **Construção do Consenso**
   - Concatena os símbolos mais frequentes de cada posição
   - Forma a string consenso final

3. **Cálculo da Distância**
   - Calcula a distância de Hamming máxima entre o consenso e todas as strings
   - Retorna o par (string_consenso, distância_máxima)

### Heurísticas Utilizadas

#### Princípio da Maioria Local
- **Premissa**: O símbolo mais comum em cada posição provavelmente minimiza os conflitos
- **Vantagem**: Simples e rápido de calcular
- **Limitação**: Não considera dependências entre posições

#### Estratégia Gananciosa
- **Abordagem**: Toma decisões localmente ótimas posição por posição
- **Benefício**: Complexidade linear O(n×L)
- **Trade-off**: Pode não encontrar o ótimo global

## Características

### Complexidade
- **Temporal**: O(n × L × |Σ|) onde n=número de strings, L=comprimento, |Σ|=tamanho do alfabeto
- **Espacial**: O(|Σ|) para contadores de frequência

### Garantias
- **Cota Superior**: Sempre produz uma solução válida
- **Qualidade**: Frequentemente próxima do ótimo para instâncias bem estruturadas
- **Determinismo**: Sempre produz o mesmo resultado para a mesma entrada

### Cenários de Uso Ideal
- **Baseline para comparação** com algoritmos mais complexos
- **Solução rápida** quando tempo é limitado
- **Pré-processamento** para outros algoritmos
- **Dados com forte consenso** por posição

## Limitações

### Casos Problemáticos
- **Empates**: Quando múltiplos símbolos têm a mesma frequência máxima
- **Distribuição Uniforme**: Quando não há consenso claro por posição
- **Correlações Posicionais**: Não considera padrões entre posições adjacentes

### Exemplo de Limitação
```
Strings: ["AB", "BA"]
Consenso: "AA" ou "BB" (dependendo do critério de desempate)
Ótimo real: qualquer string de entrada teria distância 1
Resultado baseline: distância 2
```

## Variações Possíveis

### Critérios de Desempate
- **Lexicográfico**: Escolhe o menor símbolo alfabeticamente
- **Aleatório**: Escolha randômica entre empatados
- **Primeiro Encontrado**: Mantém a ordem de processamento

### Extensões
- **Consenso Ponderado**: Dar pesos diferentes a strings
- **Consenso por Blocos**: Aplicar o princípio a subsequências
- **Múltiplos Consensos**: Gerar várias opções e escolher a melhor

## Performance Típica

### Instâncias Favoráveis
- Dados biológicos com conservação posicional
- Strings com ruído baixo (<20%)
- Alfabetos pequenos (DNA: 4 símbolos)

### Gap Médio do Ótimo
- **Dados sintéticos**: 0-15% acima do ótimo
- **Dados reais**: 5-25% acima do ótimo
- **Pior caso teórico**: Pode ser arbitrariamente distante

O algoritmo Baseline é fundamental como ponto de partida e referência, oferecendo uma solução rápida e interpretável que frequentemente serve como excelente aproximação inicial.
