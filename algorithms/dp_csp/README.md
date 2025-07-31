# DP-CSP (Dynamic Programming Closest String Problem)

The **DP-CSP** algorithm is an **exact** solution for the Closest String Problem based on **dynamic programming**. Unlike heuristics that seek approximate solutions, DP-CSP **guarantees finding the optimal solution** - the center string with the minimum possible radius - using a systematic search over all possible distance values.

## 📋 Table of Contents

- [Algorithmic Strategy](#algorithmic-strategy)
- [Detailed Operation](#detailed-operation)
- [Parameters and Configuration](#parameters-and-configuration)
- [Use Cases](#use-cases)
- [Algorithmic Analysis](#algorithmic-analysis)
- [Usage Examples](#usage-examples)
- [Limitations](#limitations)
- [CSPBench Integration](#cspbench-integration)

## 🎯 Algorithmic Strategy

### Main Approach
DP-CSP uses an **incremental binary search** strategy combined with **decision dynamic programming**:

1. **Incremental Search**: Tests increasing values of d (0, 1, 2, ...) until finding a solution
2. **Decision DP**: For each d, checks if there exists a center string with radius ≤ d
3. **DP States**: Maintains "remaining errors" vectors for each string in the dataset
4. **Construction**: Reconstructs the optimal center string when found

### Advantages
- **Exact Solution**: Guarantees finding the minimum possible radius (optimal solution)
- **Deterministic**: Always produces the same result for the same input
- **Mathematically Rigorous**: Based on solid dynamic programming theory
- **Verifiable**: Results can be easily validated

### Philosophy
DP-CSP sacrifices execution time in favor of **absolute precision**. It is the reference for checking the quality of heuristic algorithms.

## ⚙️ Detailed Operation

### Main Algorithm: Incremental Search
```
For d = 0, 1, 2, ..., max_d:
    If center_exists_with_radius(d):
        Return found center
    Else:
        Try next d
If no d works:
    Failure (should not happen)
```

### Subproblem: Decision DP
**Input**: String set S, alphabet Σ, radius d  
**Output**: Center string c such that max(H(c,s)) ≤ d, or NULL if doesn't exist  

### Dynamic Programming States
```
State: (position, remaining_errors_vector)
where:
- position: current index in string being constructed (0 to L-1)
- remaining_errors_vector: [r₁, r₂, ..., rₙ]
  rᵢ = maximum number of errors we can still make with string i
```

### State Transitions
```
Current state: (pos, [r₁, r₂, ..., rₙ])
For each character σ ∈ Σ:
    Calculate discount dᵢ = 1 if strings[i][pos] ≠ σ, else 0
    New state: (pos+1, [r₁-d₁, r₂-d₂, ..., rₙ-dₙ])
    If min(rᵢ-dᵢ) ≥ 0: state is viable
```

### Detailed Example
```
Strings: ["ACG", "ATG", "AAG"]
Alphabet: "ACGT"
Testing d = 1:

Initial state: (0, [1,1,1])  # position 0, 1 error allowed for each

Position 0:
├── Test 'A': strings[0][0]='A', strings[1][0]='A', strings[2][0]='A'
│   └── Discount: [0,0,0] → new state: (1, [1,1,1])
├── Test 'C': discounts [1,1,1] → state: (1, [0,0,0])
├── Test 'G': discounts [1,1,1] → state: (1, [0,0,0])
└── Test 'T': discounts [1,1,1] → state: (1, [0,0,0])

Position 1 (from best previous state):
├── Test 'A': discounts [1,1,1] → unfeasible (negative states)
├── Test 'C': discounts [0,1,1] → state: (2, [1,0,0])
├── Test 'G': discounts [1,0,0] → state: (2, [0,1,1])
└── Test 'T': discounts [1,0,1] → state: (2, [0,1,0])

Position 2:
└── Choose 'G' to reach final viable state

Center found: "ACG" or "ATG" (both have radius 1)
```

### Solution Reconstruction
```
From final state, follow backtrack:
- Final state: (3, [0,0,0])
- Position 2: chosen character = 'G'
- Position 1: chosen character = 'T'  
- Position 0: chosen character = 'A'
Center: "ATG"
```

## 🔧 Parameters and Configuration

### Main Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_d` | int | Auto | Upper limit for d search (uses baseline if None) |
| `max_time` | int | 300 | Timeout in seconds to prevent infinite execution |
| `warn_threshold` | int | 9 | Alert if (d+1)^n > 10^9 states |

### Automatic max_d Calculation
```python
# If max_d not provided, use baseline as upper bound
baseline = max_distance(strings[0], strings)  # distance from first string
max_d = baseline  # guarantees at least one solution exists
```

### Safety Limits
```python
# Memory: Estimate states as (d+1)^n
# For n=5, d=10: (10+1)^5 = 161M states (~1GB RAM)
# For n=6, d=10: (10+1)^6 = 1.77B states (~14GB RAM)

# Time: Monitor elapsed_time < max_time
# Cancellation: Allow interruption via SIGTERM
```

## 📊 Use Cases

### 🟢 Ideal For:
- **Reference Verification**: Validate quality of heuristic algorithms
- **Small Instances**: n ≤ 5-8 strings, lengths ≤ 20-50
- **Theoretical Analysis**: Study exact CSP properties
- **Benchmarking**: Establish lower bounds for comparison

### 🟡 Suitable For:
- **Prototypes**: Development and testing of new algorithms
- **Critical Cases**: When absolute precision is essential
- **Specific Datasets**: Short strings with small alphabets
- **Academic Research**: Controlled experiments

### 🔴 Limited For:
- **Large Instances**: n > 10 strings (exponential explosion)
- **Long Strings**: L > 100 characters (too many states)
- **Real-time Applications**: Execution can be very slow
- **Real Datasets**: Usually too large for DP

## 📈 Algorithmic Analysis

### Time Complexity
- **External Search**: O(optimal_d × F(n,d,L))
- **Internal DP**: O(L × |Σ| × unique_states)
- **Unique States**: O((d+1)<sup>n</sup>) in worst case
- **Total**: O(d × L × |Σ| × (d+1)<sup>n</sup>)

### Space Complexity
- **State Storage**: O((d+1)<sup>n</sup>)
- **Transition Table**: O(L × |Σ| × (d+1)^n)
- **Backtracking**: O(L)
- **Total**: O(L × |Σ| × (d+1)<sup>n</sup>)

### Exponential Explosion
```
n=3: (d+1)³ maximum states
n=4: (d+1)⁴ maximum states  
n=5: (d+1)⁵ maximum states
n=6: (d+1)⁶ maximum states → 1B+ states for d≥10
n=7: (d+1)⁷ maximum states → unfeasible for d>5
```

### Estimated Performance
```
n=3, L=10, d≤5:    < 1 second
n=4, L=20, d≤8:    1-10 seconds  
n=5, L=30, d≤6:    10 seconds - 2 minutes
n=6, L=40, d≤4:    1-10 minutes
n=7, L≥50:         Probably unfeasible
```

## 💡 Usage Examples

### Example 1: Small Instance
```python
from algorithms.dp_csp import DPCSPAlgorithm

strings = ["ACG", "ATG", "AAG"]
algorithm = DPCSPAlgorithm(strings, alphabet="ACGT")
center, distance, metadata = algorithm.run()

print(f"Optimal center: {center}")
print(f"Minimum radius: {distance}")
print(f"Exact solution: {metadata['exact_solution']}")
```

### Example 2: With max_d Limit
```python
# Limit search to avoid timeout
algorithm = DPCSPAlgorithm(
    strings, 
    alphabet="ACGT",
    max_d=5,  # don't test d > 5
    max_time=60  # timeout in 1 minute
)

try:
    center, distance, metadata = algorithm.run()
    print(f"Solution found: {center} with d={distance}")
except RuntimeError as e:
    print(f"DP-CSP failed: {e}")
```

### Example 3: Benchmark Verification
```python
# Compare DP-CSP (exact) vs Baseline (heuristic)
from algorithms.baseline import BaselineAlgorithm

# Exact solution
dp_center, dp_dist, _ = DPCSPAlgorithm(strings, "ACGT").run()

# Heuristic solution  
baseline_center, baseline_dist, _ = BaselineAlgorithm(strings, "ACGT").run()

print(f"DP-CSP (optimal): d={dp_dist}")
print(f"Baseline: d={baseline_dist}")
print(f"Gap: {baseline_dist - dp_dist} ({100*(baseline_dist-dp_dist)/dp_dist:.1f}%)")
```

### Example 4: Limits Analysis
```python
import time

def test_limits():
    for n in range(3, 8):
        strings = ["A"*10, "T"*10, "G"*10, "C"*10][:n]
        
        try:
            start = time.time()
            center, dist, meta = DPCSPAlgorithm(strings, "ACGT", max_time=30).run()
            elapsed = time.time() - start
            print(f"n={n}: SUCCESS d={dist} in {elapsed:.2f}s")
        except RuntimeError as e:
            print(f"n={n}: FAILED - {e}")

test_limits()
```

## ⚠️ Limitations

### Fundamental Limitations
1. **Exponential Explosion**: (d+1)<sup>n</sup> states grow exponentially
2. **Memory Limit**: Can consume gigabytes of RAM quickly
3. **Timeout**: Execution can take hours/days for large instances
4. **Scalability**: Impractical for n > 8-10 strings

### Practical Limitations
1. **Real Datasets**: Most are too large for exact DP
2. **Online Applications**: Unacceptable latency for interactive use
3. **Computational Resources**: Requires powerful machines for n > 6
4. **Implementation**: Code complexity vs simple algorithms

### Problematic Scenarios
```python
# Case 1: Many strings (exponential explosion)
strings = ["ACGT"] * 10  # n=10, unfeasible

# Case 2: High radius needed (many states)
strings = ["AAAA", "TTTT", "GGGG", "CCCC"]  # d=4, still ok for n=4

# Case 3: Long strings (many positions)
strings = ["A"*1000, "T"*1000, "G"*1000]  # L=1000, many iterations

# Case 4: Large alphabet (more transitions)
strings = ["ABC", "DEF", "GHI"]  # alphabet with 9 letters
```

### Workarounds
```python
# Conservative max_d limit
algorithm = DPCSPAlgorithm(strings, alphabet, max_d=min(5, len(strings[0])//4))

# Aggressive timeout
algorithm = DPCSPAlgorithm(strings, alphabet, max_time=30)

# Pre-filter by size
if len(strings) > 8 or len(strings[0]) > 50:
    print("Instance too large for DP-CSP, use heuristic")
else:
    result = algorithm.run()
```

## 🔗 CSPBench Integration

### Automatic Registration
```python
@register_algorithm
class DPCSPAlgorithm(CSPAlgorithm):
    name = "DP-CSP"
    supports_internal_parallel = False
    is_deterministic = True
```

### YAML Configuration
```yaml
algorithm:
  name: "DP-CSP"
  params:
    max_d: 5
    max_time: 120
```

### CLI Execution
```bash
python main.py --algorithm DP-CSP --dataset small_synthetic --max_d 3
```

### Parallelization Support
- **Internal Parallelism**: ❌ Not supported (sequential algorithm)
- **Run Parallelism**: ✅ Multiple independent executions
- **Compatibility**: ⚠️ Careful with memory consumption in parallel

### Returned Metadata
```python
metadata = {
    "iterations": 1,
    "max_d_used": 5,
    "exact_solution": True,
    "center_found": "ACGT"
}
```

### Error Handling
```python
try:
    center, distance, metadata = algorithm.run()
    # Success: optimal solution found
except RuntimeError as e:
    # Failure: timeout, memory limit, or insufficient max_d
    print(f"DP-CSP could not solve: {e}")
    # Fallback to heuristic algorithm
```

### Troubleshooting

**Problem**: "memory usage exceeded safe limit"
```
Solution: Reduce max_d or use machine with more RAM
```

**Problem**: "execution time exceeded Xs"
```
Solution: Increase max_time or reduce instance size
```

**Problem**: "(d+1)^n exceeds practical limit"
```
Solution: Use fewer strings or heuristic algorithm
```

---

**Developed for CSPBench** - Experimentation Framework for the Closest String Problem  
📚 For more information, see the [main documentation](../../README.md) of the framework.
