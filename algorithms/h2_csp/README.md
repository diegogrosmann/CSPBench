# H²-CSP (Hybrid Hierarchical Search)
> Anteriormente denominado H³-CSP (Hybrid Hierarchical Hamming Search). Renomeado para refletir foco em hierarquia + hibridização (H²).  
> Refactored: implementação unificada permanece em `algorithm.py` retornando `AlgorithmResult` conforme interface `CSPAlgorithm`.

The **H²-CSP (Hybrid Hierarchical Search)** algorithm is a **three‑layer adaptive approach** for the Closest String Problem that combines **hierarchical block decomposition**, **difficulty‑aware per‑block search techniques**, and **global hill‑climbing refinement**. It targets medium‑sized instances where structural locality can be exploited for quality vs. efficiency balance.

## 📋 Index

- [Algorithm Strategy](#algorithm-strategy)
- [Three-Layer Architecture](#three-layer-architecture)
- [Detailed Operation](#detailed-operation)
- [Parameters](#parameters)
- [Use Cases](#use-cases)
- [Complexity Analysis](#complexity-analysis)
- [Examples](#examples)
- [Limitations](#limitations)
- [Integration](#integration)

## 🎯 Algorithm Strategy

### Adaptive Divide‑and‑Conquer Philosophy
Complex CSP instances are decomposed into smaller contiguous blocks whose local difficulty guides selection of the most cost‑effective search technique per block.

### Key Advantages
- **Adaptive**: Technique automatically chosen per difficulty level
- **Hierarchical**: √L heuristic for near‑balanced contiguous block sizes
- **Hybrid**: Exhaustive search, beam search, and hill‑climbing
- **Scalable**: Controls combinatorial growth via decomposition
- **Deterministic**: Seeded execution (no stochastic divergence)

### Core Innovation: Smart-Core
Central adaptive layer measuring block difficulty (consensus max distance) and selecting:

- **Easy Blocks** → Exhaustive search (optimal locally)
- **Medium Blocks** → Reduced beam (beam_width/2)
- **Hard Blocks** → Full beam search

## 🏗️ Three-Layer Architecture

### Layer 1: Block Splitter
```
Entrada: String de comprimento L
├── Divisão em B ≈ ⌈√L⌉ blocos contíguos
├── Tamanho de bloco ≈ L/B
└── Saída: Lista de blocos [(0,b₁), (b₁,b₂), ..., (bₙ₋₁,L)]
```

**Example**: L=16 → B=4 blocks → [(0,4), (4,8), (8,12), (12,16)]

### Layer 2: Smart-Core (Adaptive Block Processing)
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

### Layer 3: Global Refine (Fusion + Local Search)
```
├── Fusão: Concatena melhores candidatos de cada bloco
├── Hill-Climbing: Refinamento posição-a-posição
└── Iteração até convergência ou limite de tempo
```

## ⚙️ Detailed Operation

### Phase 1: Hierarchical Decomposition (B-Splitter)
```python
def split_in_blocks(L: int) -> list[Block]:
    B = ceil(√L)  # Número de blocos
    block_size = ceil(L / B)  # Tamanho base
    # Distribui posições uniformemente
```

**Rationale**: √L balances granularity (too many small blocks) vs. overhead (too few large blocks).

### Phase 2: Adaptive Processing (Smart-Core)

#### Difficulty Measurement
```python
consensus = majority_vote(strings_block)
d_b = max_hamming_distance(consensus, strings_block)
```

#### Technique Selection
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

#### Exhaustive Search (Easy Blocks)
```python
# Para bloco de tamanho m, testa todos |Σ|^m candidatos
for candidate in product(alphabet, repeat=m):
    distance = max_hamming(candidate, block_strings)
    # Mantém k melhores
```

#### Beam Search (Medium / Hard Blocks)
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

### Phase 3: Fusion & Global Refinement

#### Block Fusion
```python
# Concatena o melhor candidato de cada bloco
center = "".join(best_candidates_per_block)
```

#### Hill-Climbing Refinement
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

## 🔧 Parameters

### Block Division Parameters

| Parameter | Type | Default | Description |
|-----------|------|--------|-----------|
| `auto_blocks` | bool | True | Usa divisão automática por √L |
| `min_block_size` | int | 2 | Tamanho mínimo de bloco |
| `max_blocks` | int | None | Máximo de blocos (None = automático) |

### Difficulty Thresholds

| Parameter | Type | Default | Description |
|-----------|------|--------|-----------|
| `block_small` | int | 2 | Limite para busca exaustiva |
| `block_medium` | int | 4 | Limite para beam search reduzido |
| `block_large` | int | 8 | Limite para beam search completo |

### Search Parameters

| Parameter | Type | Default | Description |
|-----------|------|--------|-----------|
| `beam_width` | int | 32 | Largura do beam search |
| `k_candidates` | int | 5 | Número de candidatos por bloco |
| `exhaustive_limit` | int | 10000 | Limite para busca exaustiva |

### Refinement Parameters

| Parameter | Type | Default | Description |
|-----------|------|--------|-----------|
| `local_iters` | int | 3 | Iterações de hill-climbing |
| `max_time` | int | 300 | Timeout em segundos |
| `seed` | int | None | Semente para reprodutibilidade |

## 📊 Use Cases

### 🟢 Ideal For
- **Instâncias Médias**: L=50-500, n=10-100 strings
- **Dados com Estrutura Local**: Padrões regionais preservados
- **Busca de Alta Qualidade**: Quando precisão é importante
- **Análise Comparativa**: Benchmark contra outros algoritmos

### 🟡 Suitable For
- **Datasets Balanceados**: Nenhuma região extremamente difícil
- **Problemas Estruturados**: Sequências biológicas, texto com padrões
- **Experimentação**: Ajuste fino de parâmetros
- **Desenvolvimento**: Prototipagem de novas técnicas

### 🔴 Less Suitable For
- **Instâncias Muito Pequenas**: Overhead de divisão é desnecessário
- **Instâncias Muito Grandes**: L>1000 pode ser lento
- **Tempo Real**: Múltiplas fases podem introduzir latência
- **Dados Completamente Aleatórios**: Estrutura local é inexistente

## 📈 Complexity Analysis

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
O brilhantismo do H²-CSP está na regra √L:
- **L=16**: 4 blocos de tamanho 4
- **L=64**: 8 blocos de tamanho 8  
- **L=256**: 16 blocos de tamanho 16
- **L=1024**: 32 blocos de tamanho 32

## 💡 Exemplos de Uso

### Exemplo 1: Configuração Básica
```python
from algorithms.h2_csp import H2CSPAlgorithm

strings = ["ACGTACGT", "AGCTACGT", "ATGTACGT", "ACATACGT"]
algorithm = H2CSPAlgorithm(strings, alphabet="ACGT")
center, distance, metadata = algorithm.run()

print(f"Centro H²-CSP: {center}")
print(f"Distância: {distance}")
print(f"Blocos processados: {len(algorithm.h2_csp_instance.blocks)}")
```

### Exemplo 2: Ajuste Fino de Parâmetros
```python
algorithm = H2CSPAlgorithm(
    strings,
    alphabet="ACGT",
    beam_width=64,
    k_candidates=8,
    local_iters=5,
    block_small=1,
    block_medium=3
)
center, distance, metadata = algorithm.run()
```

### Exemplo 3: Monitoramento de Progresso
```python
def progress_handler(message):
    print(f"[H²-CSP] {message}")

algorithm = H2CSPAlgorithm(strings, "ACGT")
algorithm.set_progress_callback(progress_handler)
center, distance, metadata = algorithm.run()
```

### Exemplo 4: Análise de Performance
```python
def benchmark_h2_csp(strings_list, label):
    print(f"\n=== {label} ===")
    algorithm = H2CSPAlgorithm(strings_list, "ACGT", max_time=60)
    
    start = time.time()
    center, distance, metadata = algorithm.run()
    elapsed = time.time() - start
    
    print(f"Strings: {len(strings_list)}, Comprimento: {len(strings_list[0])}")
    print(f"Blocos: {len(algorithm.h2_csp_instance.blocks)}")
    print(f"Resultado: d={distance} em {elapsed:.2f}s")
    print(f"Parâmetros: {metadata['parametros_usados']}")
```

### Exemplo 5: Comparação Multi-Algoritmo
```python
from algorithms.baseline import BaselineAlgorithm
from algorithms.csc import CSCAlgorithm
from algorithms.h2_csp import H2CSPAlgorithm

def compare_algorithms(strings):
    algorithms = [
        ("Baseline", BaselineAlgorithm),
        ("CSC", CSCAlgorithm),
        ("H²-CSP", H2CSPAlgorithm)
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
class H2CSPAlgorithm(CSPAlgorithm):
    name = "H²-CSP"
    supports_internal_parallel = False
    is_deterministic = True
```

### Configuração via YAML
```yaml
algorithm:
  name: "H²-CSP"
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
python main.py --algorithm "H²-CSP" --dataset synthetic
python main.py --algorithm "H²-CSP" --dataset file \
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
    "algoritmo": "H²-CSP",
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
class CustomH2CSP(H2CSPAlgorithm):
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
# H²-CSP funciona perfeitamente com otimização Optuna
study = optuna.create_study()
study.optimize(
    lambda trial: optimize_h2_csp_params(trial, dataset),
    n_trials=100
)
```

---

**Desenvolvido para CSPBench** - Framework de Experimentação para o Closest String Problem  
📚 Para mais informações, consulte a [documentação principal](../../README.md) do framework.
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
