# Contributing to CSPBench

Thank you for your interest in contributing to CSPBench! This document provides guidelines and information for contributors.

## 🤝 Ways to Contribute

- **Bug Reports**: Report issues or bugs you encounter
- **Feature Requests**: Suggest new features or improvements
- **Code Contributions**: Submit pull requests with fixes or new features
- **Documentation**: Improve or translate documentation
- **Algorithm Implementations**: Add new CSP algorithm implementations
- **Testing**: Write tests and improve test coverage

## 🚀 Getting Started

### Development Setup

1. **Fork the repository**
   ```bash
   git clone https://github.com/yourusername/CSPBench.git
   cd CSPBench
   ```

2. **Set up development environment**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   pip install -e .
   ```

3. **Install development tools**
   ```bash
   pip install pre-commit black ruff mypy pytest pytest-cov
   pre-commit install
   ```

4. **Run tests**
   ```bash
   pytest
   pytest --cov=src
   ```

## 📝 Development Guidelines

### Code Style

We use modern Python best practices:

- **Black** for code formatting
- **Ruff** for linting
- **MyPy** for type checking
- **Pytest** for testing

Run quality checks:
```bash
black .
ruff check .
mypy src/
pytest
```

### Code Structure

CSPBench follows hexagonal architecture:

```
src/
├── domain/           # Core business logic (no external dependencies)
├── application/      # Use cases and application services
├── infrastructure/   # External adapters (persistence, I/O, etc.)
└── presentation/     # User interfaces (CLI, etc.)

algorithms/           # Algorithm implementations
├── your_algorithm/
│   ├── __init__.py
│   ├── algorithm.py     # Main algorithm class
│   ├── config.py        # Configuration
│   ├── implementation.py # Implementation details
│   └── README.md        # Algorithm documentation
```

### Naming Conventions

- **Files**: `snake_case.py`
- **Classes**: `PascalCase`
- **Functions/Variables**: `snake_case`
- **Constants**: `UPPER_SNAKE_CASE`
- **Private members**: `_leading_underscore`

### Documentation

- All public functions must have docstrings
- Use type hints for function parameters and return values
- Update README.md if adding new features
- Include examples in docstrings when appropriate

## 🧩 Adding New Algorithms

### Algorithm Structure

1. **Create algorithm directory**:
   ```
   algorithms/your_algorithm/
   ├── __init__.py
   ├── algorithm.py
   ├── config.py
   ├── implementation.py
   └── README.md
   ```

2. **Implement the algorithm class**:
   ```python
   from src.domain.algorithms import CSPAlgorithm, register_algorithm

   @register_algorithm
   class YourAlgorithm(CSPAlgorithm):
       name = "YourAlgorithm"
       default_params = {
           "param1": 100,
           "param2": 0.5
       }
       is_deterministic = False
       supports_internal_parallel = True

       def __init__(self, strings, alphabet, **params):
           super().__init__(strings, alphabet, **params)

       def run(self):
           # Implementation here
           return center_string, max_distance, metadata
   ```

3. **Add configuration**:
   ```python
   # config.py
   ALGORITHM_CONFIG = {
       "name": "Your Algorithm",
       "description": "Brief description",
       "parameters": {
           "param1": {"type": "int", "range": [50, 200]},
           "param2": {"type": "float", "range": [0.1, 1.0]}
       }
   }
   ```

4. **Write documentation**:
   ```markdown
   # Your Algorithm

   ## Description
   Brief description of the algorithm...

   ## Parameters
   - `param1`: Description
   - `param2`: Description

   ## References
   - Paper citations...
   ```

### Algorithm Requirements

- Must inherit from `CSPAlgorithm`
- Must implement the `run()` method
- Must be decorated with `@register_algorithm`
- Must include proper type hints
- Must handle progress reporting
- Must include comprehensive tests

## 🧪 Testing

### Test Structure

```
tests/
├── unit/           # Unit tests
├── integration/    # Integration tests
└── conftest.py     # Test configuration
```

### Writing Tests

1. **Unit tests** for individual components:
   ```python
   def test_algorithm_initialization():
       strings = ["ACGT", "AGGT", "ATGT"]
       alg = YourAlgorithm(strings, "ACGT")
       assert alg.strings == strings
   ```

2. **Integration tests** for complete workflows:
   ```python
   def test_algorithm_execution():
       result = run_algorithm("YourAlgorithm", test_dataset)
       assert result.center_string is not None
   ```

### Test Coverage

- Aim for >80% test coverage
- Test both normal and edge cases
- Include error handling tests
- Mock external dependencies

## 📦 Pull Request Process

### Before Submitting

1. **Run all quality checks**:
   ```bash
   black .
   ruff check .
   mypy src/
   pytest --cov=src
   ```

2. **Update documentation** if needed

3. **Add/update tests** for new functionality

4. **Test with different Python versions** (3.8+)

### PR Guidelines

1. **Create descriptive title**:
   - ✅ "Add support for custom distance functions"
   - ❌ "Fix bug"

2. **Provide detailed description**:
   - What changes were made
   - Why the changes were needed
   - How to test the changes

3. **Reference related issues**:
   - "Fixes #123"
   - "Related to #456"

4. **Keep PRs focused**:
   - One feature/fix per PR
   - Avoid mixing unrelated changes

### Review Process

1. **Automated checks** must pass
2. **Code review** by maintainers
3. **Manual testing** if needed
4. **Documentation review**
5. **Final approval** and merge

## 🐛 Bug Reports

### Before Reporting

1. **Search existing issues** for duplicates
2. **Try the latest version**
3. **Isolate the problem** with minimal example

### Bug Report Template

```markdown
**Describe the bug**
A clear description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run command '...'
2. With input '...'
3. See error

**Expected behavior**
What you expected to happen.

**Environment**
- OS: [e.g., Ubuntu 20.04]
- Python version: [e.g., 3.9.1]
- CSPBench version: [e.g., 0.1.0]

**Additional context**
Add any other context about the problem.
```

## 💡 Feature Requests

### Feature Request Template

```markdown
**Is your feature request related to a problem?**
A clear description of what the problem is.

**Describe the solution you'd like**
A clear description of what you want to happen.

**Describe alternatives you've considered**
Alternative solutions or features you've considered.

**Additional context**
Add any other context about the feature request.
```

## 📜 License

By contributing to CSPBench, you agree that your contributions will be licensed under the MIT License.

## 🙏 Recognition

Contributors will be:
- Listed in the CONTRIBUTORS.md file
- Mentioned in release notes
- Credited in academic publications (if applicable)

## 📞 Getting Help

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Email**: diego.grosmann@example.com

Thank you for contributing to CSPBench! 🎉
