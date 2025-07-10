# CSC (Consensus String Clustering)

O algoritmo **Consensus String Clustering (CSC)** é uma abordagem híbrida para o Closest String Problem que combina **clusterização de strings** com **recombinação de blocos** para encontrar uma string central de alta qualidade. O algoritmo agrupa strings similares, calcula consensos locais e então recombina segmentos desses consensos para gerar candidatos otimizados.

## 📋 Índice

- [Estratégia Algorítmica](#estratégia-algorítmica)
- [Funcionamento Detalhado](#funcionamento-detalhado)
- [Parâmetros e Configuração](#parâmetros-e-configuração)
- [Casos de Uso](#casos-de-uso)
- [Análise Algorítmica](#análise-algorítmica)
- [Exemplos de Uso](#exemplos-de-uso)
- [Limitações](#limitações)
- [Integração com CSPBench](#integração-com-cspbench)

## 🎯 Estratégia Algorítmica

### Abordagem Principal
O CSC utiliza uma estratégia de **"dividir para conquistar"** combinada com **aprendizado de padrões locais**:

1. **Clusterização por Similaridade**: Agrupa strings com distâncias Hamming próximas usando DBSCAN
2. **Consenso Local**: Calcula a string consenso para cada cluster independentemente
3. **Recombinação de Blocos**: Divide consensos em segmentos e testa todas as combinações possíveis
4. **Busca Local**: Refina o melhor candidato através de otimização posição-a-posição

### Vantagens
- **Explora Estrutura Local**: Aproveita padrões regionais no dataset
- **Híbrido**: Combina clusterização não-supervisionada com busca determinística
- **Escalável**: Performance razoável mesmo com datasets grandes
- **Robusto**: Parâmetros são calculados automaticamente se não especificados

### Filosofia
O CSC assume que strings similares podem compartilhar padrões locais que, quando combinados estrategicamente, podem levar a uma solução global melhor do que consensos únicos.

## ⚙️ Funcionamento Detalhado

### Etapa 1: Preparação e Análise
```
Entrada: [ACGT, AGCT, ATGT, CCGT]
└── Análise de distâncias Hamming
└── Cálculo automático de parâmetros (d, n_blocks)
```

### Etapa 2: Clusterização (DBSCAN)
```
Parâmetro d = 2 (raio de distância)
├── Cluster 1: [ACGT, AGCT, ATGT] (strings similares)
└── Cluster 2: [CCGT] (string isolada)
```

### Etapa 3: Consenso Local
```
Cluster 1: [ACGT, AGCT, ATGT]
├── Posição 0: A,A,A → A (maioria)
├── Posição 1: C,G,T → C (primeiro mais comum)
├── Posição 2: G,C,G → G (maioria)
└── Posição 3: T,T,T → T (maioria)
Resultado: ACGT

Cluster 2: [CCGT]
Resultado: CCGT (consenso trivial)
```

### Etapa 4: Recombinação de Blocos
```
Consensos: [ACGT, CCGT]
Dividindo em n_blocks=2:
├── ACGT → ["AC", "GT"]
└── CCGT → ["CC", "GT"]

Candidatos por recombinação:
├── "AC" + "GT" = "ACGT"
├── "AC" + "GT" = "ACGT" (repetido)
├── "CC" + "GT" = "CCGT"
└── "CC" + "GT" = "CCGT" (repetido)
```

### Etapa 5: Avaliação e Busca Local
```
Melhor candidato: ACGT (menor distância máxima)
└── Busca local: testa melhorias posição-a-posição
└── Resultado final: ACGT
```

## 🔧 Parâmetros e Configuração

### Parâmetros Principais

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `d` | int | Auto | Raio de distância para clusterização DBSCAN |
| `n_blocks` | int | Auto | Número de blocos para recombinação |

### Parâmetros de Cálculo Automático

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `min_d` | int | 2 | Distância mínima para DBSCAN |
| `d_factor` | float | 0.8 | Fator da média das distâncias para calcular d |
| `min_blocks` | int | 2 | Número mínimo de blocos |
| `max_blocks` | int | 4 | Número máximo de blocos |
| `n_div` | int | 6 | Divisor do número de strings para n_blocks |
| `l_div` | int | 25 | Divisor do comprimento das strings para n_blocks |

### Cálculo Automático de Parâmetros
```python
# d (raio DBSCAN)
d = max(min_d, floor(média_distâncias_hamming * d_factor))

# n_blocks (número de blocos)
n_blocks = max(min_blocks, min(max_blocks, n_strings/n_div, L_strings/l_div))
```

## 📊 Casos de Uso

### 🟢 Ideal Para:
- **Datasets com Estrutura Local**: Strings que formam grupos naturais
- **Instâncias Médias**: 10-100 strings, comprimentos 50-500
- **Ruído Moderado**: Datasets com padrões locais preservados
- **Análise Exploratória**: Compreender agrupamentos no dataset

### 🟡 Adequado Para:
- **Datasets Balanceados**: Múltiplos grupos de tamanhos similares
- **Problemas Estruturados**: Sequências com regiões conservadas
- **Análise Comparativa**: Benchmark contra algoritmos mais simples

### 🔴 Limitado Para:
- **Datasets Muito Grandes**: >1000 strings (clusterização custosa)
- **Strings Muito Longas**: >1000 caracteres (muitos blocos)
- **Alto Ruído**: Strings completamente aleatórias
- **Tempo Real**: Execução pode ser lenta para instâncias grandes

## 📈 Análise Algorítmica

### Complexidade Temporal
- **Clusterização**: O(n² × L) onde n = número de strings, L = comprimento
- **Consenso**: O(k × m × L) onde k = clusters, m = strings por cluster
- **Recombinação**: O(k^n_blocks × L) - exponencial no número de blocos
- **Busca Local**: O(iterações × L × |alfabeto|)
- **Total**: O(n² × L + k^n_blocks × L)

### Complexidade Espacial
- **Armazenamento**: O(n × L + k × L + k^n_blocks × L)
- **Pico de Memória**: Durante geração de candidatos

### Performance Esperada
```
n=10,  L=50:   < 1 segundo
n=50,  L=100:  1-5 segundos
n=100, L=200:  5-30 segundos
n=500, L=500:  1-10 minutos
```

## 💡 Exemplos de Uso

### Exemplo 1: Configuração Básica
```python
from algorithms.csc import CSCAlgorithm

strings = ["ACGT", "AGCT", "ATGT", "CCGT"]
algorithm = CSCAlgorithm(strings, alphabet="ACGT")
center, distance, metadata = algorithm.run()

print(f"Centro: {center}")
print(f"Distância: {distance}")
print(f"Clusters encontrados: {metadata['parametros_usados']}")
```

### Exemplo 2: Parâmetros Customizados
```python
# Forçar clusterização mais agressiva
algorithm = CSCAlgorithm(
    strings,
    alphabet="ACGT",
    d=1,  # Raio menor = clusters menores
    n_blocks=3  # Mais blocos = mais combinações
)
center, distance, metadata = algorithm.run()
```

### Exemplo 3: Via Interface do CSPBench
```python
from src.core.interfaces.algorithm_interface import AlgorithmRunner

runner = AlgorithmRunner()
result = runner.run_algorithm(
    algorithm_name="CSC",
    strings=["ACGTACGT", "AGCTACGT", "ATGTACGT"],
    params={"d": 2, "n_blocks": 2}
)
```

### Exemplo 4: Análise de Metadados
```python
center, distance, metadata = algorithm.run()

print("=== Análise CSC ===")
print(f"Centro encontrado: {center}")
print(f"Distância máxima: {distance}")
print(f"Parâmetros utilizados: {metadata['parametros_usados']}")
print(f"Sucesso: {metadata['sucesso']}")

if not metadata['sucesso']:
    print("⚠️ Algoritmo falhou, fallback utilizado")
```

## ⚠️ Limitações

### Limitações Técnicas
1. **Explosão Combinatorial**: k^n_blocks candidatos podem ser muitos
2. **Sensibilidade a Parâmetros**: d e n_blocks afetam drasticamente os resultados
3. **Qualidade de Clusters**: DBSCAN pode falhar com dados esparsos
4. **Overhead de Memória**: Armazena todos os candidatos simultaneamente

### Limitações Práticas
1. **Datasets Desequilibrados**: Clusters de tamanhos muito diferentes
2. **Strings Aleatórias**: Sem estrutura local, clusterização é inútil
3. **Tempo de Execução**: Pode ser lento comparado a heurísticas simples
4. **Determinismo Limitado**: Dependente da implementação do DBSCAN

### Cenários Problemáticos
```python
# Caso 1: Todas as strings são únicas (sem clusters)
strings = ["AAAA", "TTTT", "GGGG", "CCCC"]  # d_max = 4

# Caso 2: Muitos clusters pequenos
strings = ["ACAT", "ACGT", "TCAT", "TGGT", "GCAA", "GCTT"]

# Caso 3: Strings muito longas com muitos blocos
strings = ["A"*1000, "T"*1000, "G"*1000]  # n_blocks pode ser alto
```

## 🔗 Integração com CSPBench

### Registro Automático
O algoritmo é registrado automaticamente no framework via decorador:

```python
@register_algorithm
class CSCAlgorithm(CSPAlgorithm):
    name = "CSC"
    supports_internal_parallel = False
    is_deterministic = True
```

### Configuração via YAML
```yaml
algorithm:
  name: "CSC"
  params:
    d: 3
    n_blocks: 2
```

### Execução via CLI
```bash
python main.py --algorithm CSC --dataset synthetic --d 2 --n_blocks 3
```

### Suporte a Paralelização
- **Paralelismo Interno**: ❌ Não suportado
- **Paralelismo de Runs**: ✅ Múltiplas execuções podem rodar em paralelo
- **Compatibilidade**: ✅ Funciona com batch processing e otimização

### Metadados Retornados
```python
metadata = {
    "iteracoes": 1,
    "parametros_usados": {"d": 2, "n_blocks": 2},
    "centro_encontrado": "ACGT",
    "sucesso": True,
    "fallback_usado": False  # Apenas se sucesso = False
}
```

### Troubleshooting

**Problema**: Nenhum cluster encontrado
```
Solução: Reduzir o parâmetro 'd' ou verificar se strings são muito diferentes
```

**Problema**: Muitos candidatos, execução lenta
```
Solução: Reduzir 'n_blocks' ou aumentar 'd' para clusters maiores
```

**Problema**: Qualidade ruim dos resultados
```
Solução: Ajustar parâmetros manualmente ou usar algoritmo diferente
```

---

**Desenvolvido para CSPBench** - Framework de Experimentação para o Closest String Problem  
📚 Para mais informações, consulte a [documentação principal](../../README.md) do framework.
