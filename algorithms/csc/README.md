# CSC (Consensus String Clustering)

The **Consensus String Clustering (CSC)** algorithm is a hybrid approach to the Closest String Problem that combines **string clustering** with **block recombination** to find a high-quality center string. The algorithm groups similar strings, calculates local consensus, and then recombines segments of these consensus to generate optimized candidates.

## 📋 Table of Contents

- [Algorithmic Strategy](#algorithmic-strategy)
- [Detailed Operation](#detailed-operation)
- [Parameters and Configuration](#parameters-and-configuration)
- [Use Cases](#use-cases)
- [Algorithmic Analysis](#algorithmic-analysis)
- [Usage Examples](#usage-examples)
- [Limitations](#limitations)
- [Integration with CSPBench](#integration-with-cspbench)

## 🎯 Algorithmic Strategy

### Main Approach
CSC uses a **"divide and conquer"** strategy combined with **local pattern learning**:

1. **Similarity Clustering**: Groups strings with close Hamming distances using DBSCAN
2. **Local Consensus**: Calculates consensus string for each cluster independently
3. **Block Recombination**: Divides consensus into segments and tests all possible combinations
4. **Local Search**: Refines the best candidate through position-by-position optimization

### Advantages
- **Explores Local Structure**: Leverages regional patterns in dataset
- **Hybrid**: Combines unsupervised clustering with deterministic search
- **Scalable**: Reasonable performance even with large datasets
- **Robust**: Parameters are calculated automatically if not specified

### Philosophy
CSC assumes that similar strings can share local patterns that, when strategically combined, can lead to a better global solution than single consensus.

## ⚙️ Detailed Operation

### Step 1: Preparation and Analysis
```
Input: [ACGT, AGCT, ATGT, CCGT]
└── Hamming distance analysis
└── Automatic parameter calculation (d, n_blocks)
```

### Step 2: Clustering (DBSCAN)
```
Parameter d = 2 (distance radius)
├── Cluster 1: [ACGT, AGCT, ATGT] (similar strings)
└── Cluster 2: [CCGT] (isolated string)
```

### Step 3: Local Consensus
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

### Step 4: Block Recombination
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

### Step 5: Evaluation & Local Search
```
Best candidate: ACGT (lowest maximum distance)
└── Local search: position-wise improvements
└── Final result: ACGT
```

## 🔧 Parameters and Configuration

### Main Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `d` | int | Auto | Distance radius for DBSCAN clustering |
| `n_blocks` | int | Auto | Number of blocks for recombination |

### Automatic Calculation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_d` | int | 2 | Minimum distance for DBSCAN |
| `d_factor` | float | 0.8 | Distance average factor to calculate d |
| `min_blocks` | int | 2 | Minimum number of blocks |
| `max_blocks` | int | 4 | Maximum number of blocks |
| `n_div` | int | 6 | String number divisor for n_blocks |
| `l_div` | int | 25 | String length divisor for n_blocks |

### Automatic Parameter Calculation
```python
# d (DBSCAN radius)
d = max(min_d, floor(hamming_distance_average * d_factor))

# n_blocks (number of blocks)
n_blocks = max(min_blocks, min(max_blocks, n_strings/n_div, L_strings/l_div))
```

## 📊 Use Cases

### 🟢 Ideal For:
- **Datasets with Local Structure**: Strings that form natural groups
- **Medium Instances**: 10-100 strings, lengths 50-500
- **Moderate Noise**: Datasets with preserved local patterns
- **Exploratory Analysis**: Understanding dataset clusters

### 🟡 Suitable For:
- **Balanced Datasets**: Multiple groups of similar sizes
- **Structured Problems**: Sequences with conserved regions
- **Comparative Analysis**: Benchmark against simpler algorithms

### 🔴 Limited For:
- **Very Large Datasets**: >1000 strings (expensive clustering)
- **Very Long Strings**: >1000 characters (too many blocks)
- **High Noise**: Completely random strings
- **Real Time**: Execution can be slow for large instances

## 📈 Algorithmic Analysis

### Time Complexity
- **Clustering**: O(n² × L) where n = number of strings, L = length
- **Consensus**: O(k × m × L) where k = clusters, m = strings per cluster
- **Recombination**: O(k<sup>n_blocks</sup> × L) - exponential in number of blocks
- **Local Search**: O(iterations × L × |alphabet|)
- **Total**: O(n² × L + k<sup>n_blocks</sup> × L)

### Space Complexity
- **Storage**: O(n × L + k × L + k<sup>n_blocks</sup> × L)
- **Memory Peak**: During candidate generation

### Expected Performance
```
n=10,  L=50:   < 1 second
n=50,  L=100:  1-5 seconds
n=100, L=200:  5-30 seconds
n=500, L=500:  1-10 minutes
```

## 💡 Usage Examples

### Example 1: Basic Configuration
```python
from algorithms.csc import CSCAlgorithm

strings = ["ACGT", "AGCT", "ATGT", "CCGT"]
algorithm = CSCAlgorithm(strings, alphabet="ACGT")
center, distance, metadata = algorithm.run()

print(f"Center: {center}")
print(f"Distance: {distance}")
print(f"Clusters found: {metadata['parameters_used']}")
```

### Example 2: Custom Parameters
```python
# Force more aggressive clustering
algorithm = CSCAlgorithm(
    strings,
    alphabet="ACGT",
    d=1,  # Smaller radius = smaller clusters
    n_blocks=3  # More blocks = more combinations
)
center, distance, metadata = algorithm.run()
```

### Example 3: Via CSPBench Interface
```python
from src.core.interfaces.algorithm_interface import AlgorithmRunner

runner = AlgorithmRunner()
result = runner.run_algorithm(
    algorithm_name="CSC",
    strings=["ACGTACGT", "AGCTACGT", "ATGTACGT"],
    params={"d": 2, "n_blocks": 2}
)
```

### Example 4: Metadata Analysis
```python
center, distance, metadata = algorithm.run()

print("=== CSC Analysis ===")
print(f"Center found: {center}")
print(f"Maximum distance: {distance}")
print(f"Parameters used: {metadata['parameters_used']}")
print(f"Success: {metadata['success']}")

if not metadata['success']:
    print("⚠️ Algorithm failed, fallback used")
```

## ⚠️ Limitations

### Technical Limitations
1. **Combinatorial Explosion**: ``|consensus|^{n_blocks}`` (mitigated by ``max_candidates``)
2. **Parameter Sensitivity**: ``d`` and ``n_blocks`` directly shape search space
3. **Cluster Quality Dependency**: Weak or noisy clustering reduces benefit
4. **Candidate Memory Footprint**: All generated candidates kept before evaluation (consider streaming optimization in future)

### Practical Limitations
1. **Highly Unbalanced Datasets**: Single giant cluster negates diversification
2. **Random Strings**: No local structure → almost pure fallback mode
3. **Runtime vs Simple Heuristics**: Overhead of clustering + recombination
4. **DBSCAN Determinism**: Deterministic here (no randomness) but distance scaling still impacts cluster outcomes



## 🔗 CSPBench Integration

### Automatic Registration
The algorithm is automatically registered in the framework via decorator:

```python
@register_algorithm
class CSCAlgorithm(CSPAlgorithm):
    name = "CSC"
    supports_internal_parallel = False
    is_deterministic = True
```

### YAML Configuration
```yaml
algorithm:
  name: "CSC"
  params:
    d: 3
    n_blocks: 2
```

### CLI Execution
```bash
python main.py --algorithm CSC --dataset synthetic --d 2 --n_blocks 3
```

### Parallelization Support
- **Internal Parallelism**: ❌ Not supported
- **Run Parallelism**: ✅ Multiple executions can run in parallel
- **Compatibility**: ✅ Works with batch processing and optimization

### Returned Metadata
Important metadata keys (subset):
* ``d`` / ``n_blocks`` / ``d_auto`` / ``n_blocks_auto``
* ``n_clusters`` / ``fallback_used`` / ``degraded_mode``
* ``candidates_generated`` / ``candidates_evaluated`` / ``candidates_truncated``
* ``local_search_iterations`` / ``execution_time``
* ``avg_distance`` / ``total_distance`` / ``best_candidate_initial_distance``

### Troubleshooting

**Problem**: No clusters found
```
Cause: Radius d too small or strings mutually far.
Action: Increase d or allow automatic parameter inference.
Fallback: Global consensus + local search (see metadata.fallback_used=True)
```

**Problem**: Excessive runtime / many candidates
```
Action: Decrease n_blocks or lower max_candidates; alternatively increase d to merge clusters.
```

**Problem**: Weak solution quality
```
Action: Manually tune d; raise n_blocks slightly (trade-off); compare with Baseline algorithm.
```

---

**Developed for CSPBench** - Experimentation Framework for the Closest String Problem  
📚 For more information, see the [main documentation](../../README.md) of the framework.
