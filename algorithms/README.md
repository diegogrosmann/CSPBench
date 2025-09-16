# Guide: Adding New Algorithms to CSPBench

CSPBench provides a standardized interface and detailed documentation (Google-style docstrings) across the codebase, making it easy to add new algorithms without changing the core app (`main.py`). Registration is automatic via a decorator.

## Recommended Structure

Each algorithm should have its own folder under `algorithms/`:

```
algorithms/
├── my_algorithm/
│   ├── __init__.py          # Re-exports the algorithm class
│   ├── algorithm.py         # Wrapper that implements `CSPAlgorithm`
│   ├── config.py            # Algorithm-specific default parameters
│   └── implementation.py    # Core algorithm implementation
```

## Step-by-Step Integration

1. Create the algorithm folder
    ```bash
    mkdir algorithms/my_algorithm
    ```

2. Create `config.py`
    ```python
    # Example
    MY_ALGO_DEFAULTS = {
        "param1": "default_value",
        "param2": 42,
        "max_time": 300.0,
    }
    ```

3. Implement the core logic in `implementation.py`
    ```python
    def my_algorithm_core(strings, alphabet, report=None, **params):
        """
        Run the core of MyAlgorithm.

        Args:
            strings (list[str]): Input strings.
            alphabet (str): Alphabet used.
            report (callable, optional): Progress reporter callable taking (progress: float, message: str, **data).
            **params: Algorithm-specific parameters.

        Returns:
            str: Found center string.
        """
        if report:
            report(0.0, "Starting…")
        center = alphabet[0] * len(strings[0])  # toy example
        if report:
            report(1.0, "Finished.")
        return center
    ```

4. Create the wrapper `algorithm.py`
    ```python
    from src.domain.algorithms import CSPAlgorithm, AlgorithmResult, register_algorithm
    from .config import MY_ALGO_DEFAULTS
    from .implementation import my_algorithm_core

    @register_algorithm
    class MyAlgorithm(CSPAlgorithm):
        """Wrapper integrating MyAlgorithm into CSPBench."""

        name = "MyAlgorithm"
        default_params = MY_ALGO_DEFAULTS

        def run(self) -> AlgorithmResult:
            # Use `self._monitor` if present to report progress
            reporter = None
            if getattr(self, "_monitor", None) and hasattr(self._monitor, "on_progress"):
                reporter = lambda p, m, **d: self._monitor.on_progress(p, m, **d)

            center = my_algorithm_core(self.strings, self.alphabet, report=reporter, **self.params)

            # DistanceCalculator is injected by the framework and exposed via attribute access
            max_dist = self.max_distance(center)  # delegated to DistanceCalculator

            result: AlgorithmResult = {
                "success": True,
                "center_string": center,
                "max_distance": max_dist,
                "parameters": self.get_actual_params(),
                "error": None,
                "metadata": {
                    "internal_jobs": self.internal_jobs,
                },
            }
            return result
    ```

5. Expose the algorithm in `__init__.py`
    ```python
    from .algorithm import MyAlgorithm  # noqa: F401
    ```

6. (Optional) Add a README.md describing the heuristic, parameters, and usage.

## Notes

- All algorithms must implement the `CSPAlgorithm` interface and use the `@register_algorithm` decorator.
- The algorithm will automatically appear in the CLI/Web menus via the global registry.
- See existing examples for documentation patterns: `algorithms/baseline/`, `algorithms/blf_ga/`, `algorithms/csc/`, `algorithms/h2_csp/`, `algorithms/dp_csp/`.
- Prefer English for code, comments, and docs. Use Google-style docstrings for auto-generated documentation.
