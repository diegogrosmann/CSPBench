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

The Baseline algorithm **has no configurable parameters**, ensuring:
- Total reproducibility
- Ease of use
- No tuning required
- Consistent behavior

## 💻 Usage Example

### **Basic Usage**
```python
from algorithms.baseline.algorithm import BaselineAlg

# Example dataset
strings = ["ACGTACGT", "AGGTACGT", "ACGTAAGT"]
alphabet = "ACGT"

# Create and run algorithm
algorithm = BaselineAlg(strings, alphabet)
center, distance, metadata = algorithm.run()

print(f"Found center: {center}")
print(f"Maximum distance: {distance}")
print(f"Metadata: {metadata}")
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

### **Pseudocode**
```
function baseline_consensus(strings, alphabet):
    L = length(strings[0])
    consensus = ""
    
    for position in range(L):
        # Count frequencies
        counts = {}
        for symbol in alphabet:
            counts[symbol] = 0
        
        for string in strings:
            symbol = string[position]
            counts[symbol] += 1
        
        # Find most frequent symbol
        max_count = 0
        best_symbol = alphabet[0]  # Tie-breaking
        
        for symbol in alphabet:
            if counts[symbol] > max_count:
                max_count = counts[symbol]
                best_symbol = symbol
        
        consensus += best_symbol
    
    return consensus
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

## 🔗 Integration with CSPBench

Baseline is fully integrated with the framework through:

- **Automatic Registration**: Detected via `@register_algorithm`
- **Standardized Interface**: Implements `CSPAlgorithm`
- **Parallel Execution**: Compatible with execution system
- **Reports**: Generates structured metadata
- **Monitoring**: Supports progress callbacks

---

*Baseline: The solid foundation for CSP algorithm comparison - simple, fast and reliable.*
