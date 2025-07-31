# H³-CSP (Hybrid Hierarchical Hamming Search)

O algoritmo **H³-CSP (Hybrid Hierarchical Hamming Search)** é uma abordagem sofisticada de **três camadas** para o Closest String Problem que combina **decomposição hierárquica**, **seleção adaptativa de técnicas** e **refinamento global**. O algoritmo é especialmente eficaz para instâncias de tamanho médio e oferece um excelente equilíbrio entre qualidade da solução e eficiência computacional.

## 📋 Índice

- [Estratégia Algorítmica](#estratégia-algorítmica)
- [Arquitetura em Três Camadas](#arquitetura-em-três-camadas)
- [Funcionamento Detalhado](#funcionamento-detalhado)
- [Parâmetros e Configuração](#parâmetros-e-configuração)
- [Casos de Uso](#casos-de-uso)
- [Análise Algorítmica](#análise-algorítmica)
- [Exemplos de Uso](#exemplos-de-uso)
- [Limitações](#limitações)
- [Integração com CSPBench](#integração-com-cspbench)

## 🎯 Estratégia Algorítmica

### Filosofia "Dividir para Conquistar Adaptativo"
O H³-CSP se baseia na premissa de que **problemas complexos podem ser decompostos em subproblemas mais simples**, onde cada subproblema pode ser resolvido com a **técnica mais apropriada** para sua dificuldade específica.

### Vantagens Principais
- **Adaptativo**: Seleciona automaticamente a melhor técnica para cada subproblema
- **Hierárquico**: Divisão inteligente em blocos contíguos usando a "regra √L"
- **Híbrido**: Combina busca exaustiva, beam search e hill-climbing
- **Escalável**: Performance consistente em diferentes tamanhos de instância
- **Determinístico**: Resultados reproduzíveis com mesma configuração

### Core Innovation: Smart-Core
O diferencial do H³-CSP é o **Smart-Core**, que analisa a dificuldade de cada bloco (medida pela distância do consenso) e seleciona dinamicamente a técnica de busca mais eficiente:

- **Blocos Fáceis** → Busca Exaustiva (solução ótima)
- **Blocos Médios** → Beam Search Reduzido (eficiência)
- **Blocos Difíceis** → Beam Search Completo (exploração ampla)

## 🏗️ Arquitetura em Três Camadas

### Camada 1: B-Splitter (Block Splitter)
```
Entrada: String de comprimento L
├── Divisão em B ≈ ⌈√L⌉ blocos contíguos
├── Tamanho de bloco ≈ L/B
└── Saída: Lista de blocos [(0,b₁), (b₁,b₂), ..., (bₙ₋₁,L)]
```

**Exemplo**: L=16 → B=4 blocos → [(0,4), (4,8), (8,12), (12,16)]

### Camada 2: Smart-Core (Adaptive Block Processing)
```
Para cada bloco:
├── Calcula consenso local
├── Mede dificuldade d_b = max_distance(consenso, bloco_strings)
├── Seleciona técnica baseada em d_b:
│   ├── d_b ≤ 2: Busca Exaustiva
│   ├── 2 < d_b ≤ 4: Beam Search (largura/2)
│   └── d_b > 4: Beam Search (largura completa)
└── Gera k candidatos para o bloco
```

### Camada 3: Global Refine (Fusion + Local Search)
```
├── Fusão: Concatena melhores candidatos de cada bloco
├── Hill-Climbing: Refinamento posição-a-posição
└── Iteração até convergência ou limite de tempo
```

## ⚙️ Funcionamento Detalhado

### Fase 1: Decomposição Hierárquica (B-Splitter)
```python
def split_in_blocks(L: int) -> list[Block]:
    B = ceil(√L)  # Número de blocos
    block_size = ceil(L / B)  # Tamanho base
    # Distribui posições uniformemente
```

**Intuição**: A regra √L otimiza o trade-off entre granularidade (muitos blocos pequenos) e eficiência (poucos blocos grandes).

### Fase 2: Processamento Adaptativo (Smart-Core)

#### Análise de Dificuldade
```python
consensus = majority_vote(strings_block)
d_b = max_hamming_distance(consensus, strings_block)
```

#### Seleção de Técnica
```python
if d_b <= block_small:          # Padrão: 2
    # Busca Exaustiva
    candidates = all_strings_of_length_m(alphabet, block_length)
    return best_k_candidates(candidates)
    
elif d_b <= block_medium:       # Padrão: 4
    # Beam Search Reduzido
    return beam_search(beam_width // 2)
    
else:                           # d_b > 4
    # Beam Search Completo
    return beam_search(beam_width)
```

#### Busca Exaustiva (Blocos Pequenos)
```python
# Para bloco de tamanho m, testa todos |Σ|^m candidatos
for candidate in product(alphabet, repeat=m):
    distance = max_hamming(candidate, block_strings)
    # Mantém k melhores
```

#### Beam Search (Blocos Grandes)
```python
beam = [""]  # Prefixo vazio
for position in range(block_length):
    new_beam = []
    for prefix in beam:
        for char in alphabet:
            new_prefix = prefix + char
            score = evaluate_partial(new_prefix)
            new_beam.append((score, new_prefix))
    beam = select_best(new_beam, beam_width)
```

### Fase 3: Fusão e Refinamento Global

#### Fusão de Blocos
```python
# Concatena o melhor candidato de cada bloco
center = "".join(best_candidates_per_block)
```

#### Hill-Climbing Local
```python
improved = True
while improved:
    improved = False
    for position in range(len(center)):
        for new_char in alphabet_at_position[position]:
            if try_substitution_improves(position, new_char):
                center[position] = new_char
                improved = True
                break
```

## 🔧 Parâmetros e Configuração

### Parâmetros de Divisão de Blocos

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `auto_blocks` | bool | True | Usa divisão automática por √L |
| `min_block_size` | int | 2 | Tamanho mínimo de bloco |
| `max_blocks` | int | None | Máximo de blocos (None = automático) |

### Parâmetros de Dificuldade

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `block_small` | int | 2 | Limite para busca exaustiva |
| `block_medium` | int | 4 | Limite para beam search reduzido |
| `block_large` | int | 8 | Limite para beam search completo |

### Parâmetros de Busca

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `beam_width` | int | 32 | Largura do beam search |
| `k_candidates` | int | 5 | Número de candidatos por bloco |
| `exhaustive_limit` | int | 10000 | Limite para busca exaustiva |

### Parâmetros de Refinamento

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `local_iters` | int | 3 | Iterações de hill-climbing |
| `max_time` | int | 300 | Timeout em segundos |
| `seed` | int | None | Semente para reprodutibilidade |

## 📊 Casos de Uso

### 🟢 Ideal Para:
- **Instâncias Médias**: L=50-500, n=10-100 strings
- **Dados com Estrutura Local**: Padrões regionais preservados
- **Busca de Alta Qualidade**: Quando precisão é importante
- **Análise Comparativa**: Benchmark contra outros algoritmos

### 🟡 Adequado Para:
- **Datasets Balanceados**: Nenhuma região extremamente difícil
- **Problemas Estruturados**: Sequências biológicas, texto com padrões
- **Experimentação**: Ajuste fino de parâmetros
- **Desenvolvimento**: Prototipagem de novas técnicas

### 🔴 Limitado Para:
- **Instâncias Muito Pequenas**: Overhead de divisão é desnecessário
- **Instâncias Muito Grandes**: L>1000 pode ser lento
- **Tempo Real**: Múltiplas fases podem introduzir latência
- **Dados Completamente Aleatórios**: Estrutura local é inexistente

## 📈 Análise Algorítmica

### Complexidade por Fase

#### Fase 1: B-Splitter
- **Tempo**: O(L) - divisão linear
- **Espaço**: O(B) - lista de blocos

#### Fase 2: Smart-Core
- **Busca Exaustiva**: O(|Σ|^m × k) por bloco pequeno
- **Beam Search**: O(m × |Σ| × beam_width) por bloco grande
- **Total**: O(B × max(|Σ|^m_avg, m_avg × |Σ| × beam_width))

#### Fase 3: Global Refine
- **Fusão**: O(B) - concatenação
- **Hill-Climbing**: O(local_iters × L × |Σ|)

### Complexidade Total
```
Melhor caso: O(L + B × |Σ| × beam_width + local_iters × L × |Σ|)
Pior caso: O(L + B × |Σ|^(L/B) + local_iters × L × |Σ|)
Caso típico: O(L × |Σ| × beam_width)
```

### Performance Estimada
```
L=50,  n=20:   1-5 segundos
L=100, n=50:   5-15 segundos
L=200, n=100:  15-60 segundos
L=500, n=200:  1-5 minutos
```

### Scaling com √L
O brilhantismo do H³-CSP está na regra √L:
- **L=16**: 4 blocos de tamanho 4
- **L=64**: 8 blocos de tamanho 8  
- **L=256**: 16 blocos de tamanho 16
- **L=1024**: 32 blocos de tamanho 32

## 💡 Exemplos de Uso

### Exemplo 1: Configuração Básica
```python
from algorithms.h3_csp import H3CSPAlgorithm

strings = ["ACGTACGT", "AGCTACGT", "ATGTACGT", "ACATACGT"]
algorithm = H3CSPAlgorithm(strings, alphabet="ACGT")
center, distance, metadata = algorithm.run()

print(f"Centro H³-CSP: {center}")
print(f"Distância: {distance}")
print(f"Blocos processados: {len(algorithm.h3_csp_instance.blocks)}")
```

### Exemplo 2: Ajuste Fino de Parâmetros
```python
# Para instâncias maiores, aumentar beam_width
algorithm = H3CSPAlgorithm(
    strings,
    alphabet="ACGT",
    beam_width=64,        # Exploração mais ampla
    k_candidates=8,       # Mais candidatos por bloco
    local_iters=5,        # Refinamento mais intenso
    block_small=1,        # Busca exaustiva mais conservadora
    block_medium=3        # Ajuste dos limiares
)
center, distance, metadata = algorithm.run()
```

### Exemplo 3: Monitoramento de Progresso
```python
def progress_handler(message):
    print(f"[H³-CSP] {message}")

algorithm = H3CSPAlgorithm(strings, "ACGT")
algorithm.set_progress_callback(progress_handler)
center, distance, metadata = algorithm.run()
```

### Exemplo 4: Análise de Performance
```python
import time

def benchmark_h3_csp(strings_list, label):
    print(f"\n=== {label} ===")
    algorithm = H3CSPAlgorithm(strings_list, "ACGT", max_time=60)
    
    start = time.time()
    center, distance, metadata = algorithm.run()
    elapsed = time.time() - start
    
    print(f"Strings: {len(strings_list)}, Comprimento: {len(strings_list[0])}")
    print(f"Blocos: {len(algorithm.h3_csp_instance.blocks)}")
    print(f"Resultado: d={distance} em {elapsed:.2f}s")
    print(f"Parâmetros: {metadata['parametros_usados']}")

# Teste com diferentes tamanhos
small_strings = ["ACGT"] * 5
medium_strings = ["A"*20, "T"*20, "G"*20, "C"*20] * 3
large_strings = ["A"*50] * 10

benchmark_h3_csp(small_strings, "Instância Pequena")
benchmark_h3_csp(medium_strings, "Instância Média") 
benchmark_h3_csp(large_strings, "Instância Grande")
```

### Exemplo 5: Comparação Multi-Algoritmo
```python
from algorithms.baseline import BaselineAlgorithm
from algorithms.csc import CSCAlgorithm

def compare_algorithms(strings):
    algorithms = [
        ("Baseline", BaselineAlgorithm),
        ("CSC", CSCAlgorithm),
        ("H³-CSP", H3CSPAlgorithm)
    ]
    
    results = {}
    for name, AlgClass in algorithms:
        start = time.time()
        center, distance, _ = AlgClass(strings, "ACGT").run()
        elapsed = time.time() - start
        results[name] = {"distance": distance, "time": elapsed}
    
    print("=== Comparação de Algoritmos ===")
    for name, result in results.items():
        print(f"{name:10}: d={result['distance']} em {result['time']:.2f}s")

# Executar comparação
test_strings = ["ACGTACGT", "AGCTACGT", "ATGTACGT", "AAGTACGT"]
compare_algorithms(test_strings)
```

## ⚠️ Limitações

### Limitações Técnicas
1. **Overhead de Divisão**: Para instâncias muito pequenas (L<10), a divisão pode ser contraproducente
2. **Dependência de Estrutura**: Performance degrada com dados completamente aleatórios
3. **Memória para Candidatos**: k_candidates × B pode consumir memória significativa
4. **Determinismo vs Qualidade**: Pode ficar preso em ótimos locais durante hill-climbing

### Limitações Práticas
1. **Configuração de Parâmetros**: Muitos parâmetros para ajustar
2. **Compreensão Complexa**: Algoritmo mais difícil de entender que baselines
3. **Debug Complexo**: Três fases dificultam identificação de problemas
4. **Sensibilidade a Limiares**: block_small e block_medium afetam drasticamente performance

### Cenários Problemáticos
```python
# Caso 1: Strings muito curtas (overhead desnecessário)
strings = ["AC", "AT", "AG", "AA"]  # L=2, apenas 1-2 blocos

# Caso 2: Strings completamente aleatórias
strings = ["AAAA", "TTTT", "GGGG", "CCCC"]  # Sem estrutura local

# Caso 3: Um bloco dominante
strings = ["A"*50 + "TTTTT", "A"*50 + "GGGGG"]  # 1º bloco trivial, 2º difícil

# Caso 4: Limiares mal configurados
# block_small=0 → tudo usa beam search (perde precisão)
# block_small=10 → tudo usa busca exaustiva (muito lento)
```

### Troubleshooting
```python
# Problema: Muito lento
solution = "Reduzir beam_width ou aumentar block_small"

# Problema: Qualidade ruim
solution = "Aumentar k_candidates ou local_iters"

# Problema: Timeout frequente
solution = "Reduzir max_time ou simplificar parâmetros"

# Problema: Uso excessivo de memória
solution = "Reduzir k_candidates ou beam_width"
```

## 🔗 Integração com CSPBench

### Registro Automático
```python
@register_algorithm
class H3CSPAlgorithm(CSPAlgorithm):
    name = "H³-CSP"
    supports_internal_parallel = False
    is_deterministic = True
```

### Configuração via YAML
```yaml
algorithm:
  name: "H³-CSP"
  params:
    beam_width: 32
    k_candidates: 5
    local_iters: 3
    block_small: 2
    block_medium: 4
    max_time: 300
```

### Execução via CLI
```bash
# Execução básica
python main.py --algorithm "H³-CSP" --dataset synthetic

# Com parâmetros customizados
python main.py --algorithm "H³-CSP" --dataset file \
  --beam_width 64 --k_candidates 8 --local_iters 5
```

### Suporte a Paralelização
- **Paralelismo Interno**: ❌ Não suportado (algoritmo sequencial)
- **Paralelismo de Runs**: ✅ Múltiplas execuções independentes
- **Compatibilidade**: ✅ Funciona perfeitamente com batch processing

### Metadados Detalhados
```python
metadata = {
    "iteracoes": 1,
    "algoritmo": "H³-CSP",
    "parametros_usados": {
        "beam_width": 32,
        "k_candidates": 5,
        "local_iters": 3,
        "block_small": 2,
        "block_medium": 4,
        "auto_blocks": True,
        "max_time": 300
    },
    "centro_encontrado": "ACGTACGT"
}
```

### Advanced Usage: Custom Block Strategy
```python
# Implementação customizada (futuro)
class CustomH3CSP(H3CSPAlgorithm):
    def __init__(self, strings, alphabet, **params):
        # Override block division strategy
        params["auto_blocks"] = False
        params["custom_blocks"] = self.define_custom_blocks(strings)
        super().__init__(strings, alphabet, **params)
    
    def define_custom_blocks(self, strings):
        # Custom logic for domain-specific block division
        pass
```

### Integration with Optimization
```python
# H³-CSP funciona perfeitamente com otimização Optuna
study = optuna.create_study()
study.optimize(
    lambda trial: optimize_h3_csp_params(trial, dataset),
    n_trials=100
)
```

---

**Desenvolvido para CSPBench** - Framework de Experimentação para o Closest String Problem  
📚 Para mais informações, consulte a [documentação principal](../../README.md) do framework.
