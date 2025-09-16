# Baseline: Greedy Consensus Algorithm

The **Baseline** is a simple and efficient deterministic algorithm that implements a greedy consensus strategy to solve the Closest String Problem. It serves as a fundamental reference for comparison with more sophisticated methods.

## 📊 Overview

### **Main Strategy**
- **Position-wise Consensus**: For each position, choose the most frequent symbol
- **Greedy Decision**: Make locally optimal decisions without considering global impact
- **Deterministic**: Always produces the same result for the same input
- **Efficiency**: Linear execution in O(n × L × |Σ|)

### **Operation**
1. For each position i ∈ [0, L-1]:
   - Count the frequency of each alphabet symbol
   - Select the symbol with highest frequency
   - In case of tie, choose the first symbol alphabetically
2. Build consensus string by concatenating chosen symbols
3. Calculate maximum Hamming distance to all input strings

## 🔧 Technical Characteristics

### **Complexity**
- **Time**: O(n × L × |Σ|)
  - n: number of strings
  - L: string length  
  - |Σ|: alphabet size
- **Space**: O(|Σ|) for counters + O(L) for result

### **Properties**
- ✅ **Deterministic**: Always produces the same result
- ✅ **Fast**: Near-instantaneous execution
- ✅ **Simple**: Direct and comprehensible implementation
- ✅ **Stable**: No parameters to adjust
- ❌ **Quality**: May not find global optimum
- ❌ **Independence**: Does not consider dependencies between positions

## 🎯 Casos de Uso

### **✅ When to Use**
- **Comparison Baseline**: Establish baseline for other algorithms
- **Fast Execution**: When time is extremely limited
- **Strong Consensus Data**: Sequences with well-conserved positions
- **Preprocessing**: Initial solution for iterative algorithms
- **Validation**: Verify basic framework functionality

### **❌ Limitations**
- **Local Optima**: May get stuck in suboptimal solutions
- **Ties**: Arbitrary resolution may impact quality
- **Noisy Data**: Degraded performance with high noise
- **Dependencies**: Ignores correlations between positions

## 🧮 Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tie_break` | `str` | `"lex"` | Strategy to resolve ties between equally good symbols at a position. Options: `lex` (lexicographic min), `first` (first encountered), `random` (uniform choice using algorithm RNG/seed). |

Notes:
- `random` maintains reproducibility when a `seed` is provided to the algorithm constructor.
- Changing `tie_break` can produce different but distance-equivalent consensus strings when multiple optimal symbols exist per position.

## 💻 Usage Example

### **Basic Usage**
```python
from algorithms.baseline import BaselineAlg
from src.domain.distance import create_distance_calculator

strings = ["ACGTACGT", "AGGTACGT", "ACGTAAGT"]
alphabet = "ACGT"
calc = create_distance_calculator("hamming", strings)

alg = BaselineAlg(strings=strings, alphabet=alphabet, distance_calculator=calc, seed=42, tie_break="lex")
result = alg.run()

print("Center:", result["center_string"])            # consensus string
print("Max distance:", result["max_distance"])       # d*
print("Parameters:", result["parameters"])           # actual params (incl. tie_break)
print("Execution time:", result["metadata"]["execution_time"])  # seconds
```

### **Via Framework**
```bash
# CLI execution
python main.py --algorithms Baseline --dataset synthetic

# Silent execution
python main.py --silent --algorithms Baseline --dataset synthetic --num-execs 1
```

### **Batch (YAML)**
```yaml
algorithms: ["Baseline"]
task:
  type: "execution"
  execution:
    executions:
      - name: "Baseline Test"
        dataset: dataset_1
        runs_per_algorithm_per_base: 1  # Deterministic
        timeout: 30
```

## 🔬 Algorithmic Analysis

### **Pseudocode (Frequency-Based Greedy Consensus)**
```
function greedy_consensus(strings, alphabet, tie_break="lex"):
    L = length(strings[0])
    consensus = []
    for pos in 0..L-1:
        # Count symbol frequencies at this column
        counts = {c: 0 for c in alphabet}
        for s in strings:
            counts[s[pos]] += 1
        max_count = max(counts.values())
        candidates = [c for c in alphabet if counts[c] == max_count]
        if len(candidates) == 1:
            chosen = candidates[0]
        else:
            if tie_break == 'lex':
                chosen = min(candidates)
            elif tie_break == 'first':
                chosen = candidates[0]
            elif tie_break == 'random':
                chosen = random_choice(candidates)
            else:
                chosen = min(candidates)
        consensus.append(chosen)
    return join(consensus)
```

### **Mathematical Analysis**
For a position i, let f(s,i) be the frequency of symbol s at position i:
- Choice: argmax_s f(s,i)
- Expected distance per position: ≈ n × (1 - max_s(f(s,i)/n))
- Total expected distance: Σ_i n × (1 - max_s(f(s,i)/n))

## 🎨 Visualizations

### **Consensus Analysis**
```python
# Generate consensus heatmap by position
import matplotlib.pyplot as plt
import numpy as np

def visualize_consensus(strings, alphabet):
    L = len(strings[0])
    n = len(strings)
    
    # Frequency matrix
    freq_matrix = np.zeros((len(alphabet), L))
    
    for i, symbol in enumerate(alphabet):
        for pos in range(L):
            count = sum(1 for s in strings if s[pos] == symbol)
            freq_matrix[i, pos] = count / n
    
    # Plot heatmap
    plt.imshow(freq_matrix, aspect='auto', cmap='viridis')
    plt.xlabel('Position')
    plt.ylabel('Symbol')
    plt.colorbar(label='Frequency')
    plt.title('Consensus by Position')
    plt.show()
```

## ✅ Testing & Quality

High-level tests (pytest) ensure:
- Correct consensus on controlled examples
- Deterministic output for same seed & params
- Tie-breaking strategies produce expected alternatives
- Proper error handling for empty inputs / alphabet
- Metadata fields presence & basic invariants

Run all tests:
```bash
.venv/bin/python -m pytest -k baseline -v
```

## 🔗 Integration with CSPBench

Baseline is fully integrated with the framework through:

- **Automatic Registration**: Detected via `@register_algorithm`
- **Standardized Interface**: Implements `CSPAlgorithm`
- **Parallel Execution**: Compatible with execution system
- **Reports**: Generates structured metadata
- **Monitoring**: Supports progress callbacks

---

*Baseline: The solid foundation for CSP algorithm comparison - simple, fast and reliable.*
