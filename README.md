# CSPBench: A Framework for Closest String Problem Algorithms

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-CC%20BY--NC--SA%204.0-orange.svg)](LICENSE)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

CSPBench is a comprehensive framework for implementing, testing, and benchmarking algorithms for the Closest String Problem (CSP), a fundamental challenge in computational biology and bioinformatics.

## 🎯 Features

- **Hexagonal Architecture**: Clean, maintainable codebase following Domain-Driven Design principles
- **Multiple Algorithm Support**: Implementations of various CSP algorithms (Baseline, BLF-GA, CSC, H³-CSP, DP-CSP)
- **Flexible Dataset Management**: Support for synthetic, file-based, and NCBI Entrez datasets
- **Comprehensive Benchmarking**: Execution, optimization, and sensitivity analysis capabilities
- **Modern CLI Interface**: User-friendly command-line interface with interactive menus
- **Rich Reporting**: Detailed results with visualizations and statistical analysis
- **Docker Support**: Containerized execution for reproducible research
- **Extensible Design**: Easy to add new algorithms and datasets

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/diegogrosmann/CSPBench.git
cd CSPBench

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```bash
# Interactive menu
python main.py

# List available algorithms
python main.py --algorithms

# Generate synthetic datasets
python main.py --datasetsave

# Run a batch configuration
python main.py batches/processamento_padrao.yaml

# Quick algorithm test
python main.py test
```

### Docker Usage

```bash
# Build and run with Docker
make docker-build
make docker-run

# Or with docker-compose
make compose-up
```

## 📖 Documentation

### Closest String Problem

The Closest String Problem is a well-known NP-hard problem in computational biology. Given a set of strings, the goal is to find a center string that minimizes the maximum Hamming distance to all input strings.

**Formal Definition:**
- Input: Set of strings S = {s₁, s₂, ..., sₙ} of equal length L over alphabet Σ
- Output: String c that minimizes max{d(c, sᵢ) : sᵢ ∈ S}

### Supported Algorithms

| Algorithm | Type | Deterministic | Description |
|-----------|------|---------------|-------------|
| **Baseline** | Greedy | ✅ | Simple consensus algorithm |
| **BLF-GA** | Genetic Algorithm | ❌ | Block-based genetic algorithm |
| **CSC** | Heuristic | ❌ | Closest String with Constraints |
| **H³-CSP** | Heuristic | ❌ | Hierarchical heuristic approach |
| **DP-CSP** | Dynamic Programming | ✅ | Exact algorithm for small instances |

### Dataset Types

#### Synthetic Datasets
Generate artificial data for controlled experiments:
```yaml
- id: synthetic_example
  tipo: "synthetic"
  parametros:
    n: 20              # Number of sequences
    L: 50              # Sequence length
    alphabet: "ACGT"   # DNA alphabet
    noise: 0.1         # Noise level (0.0-1.0)
    seed: 42           # For reproducibility
```

#### File-based Datasets
Load existing FASTA files:
```yaml
- id: file_example
  tipo: "file"
  parametros:
    filename: "example.fasta"
```

#### NCBI Entrez Datasets
Download data directly from NCBI:
```yaml
- id: ncbi_example
  tipo: "entrez"
  parametros:
    query: "COIGene AND 600:650[SLEN]"
    db: "nucleotide"
    retmax: 20
```

## 🔬 Batch Configurations

CSPBench supports three types of batch operations:

### 1. Execution Batches
Run algorithms with fixed parameters:
```bash
python main.py batches/processamento_padrao.yaml
```

### 2. Optimization Batches
Optimize hyperparameters using Optuna:
```bash
python main.py batches/otimizacao_padrao.yaml
```

### 3. Sensitivity Analysis
Analyze parameter sensitivity using SALib:
```bash
python main.py batches/sensibilidade_padrao.yaml
```

## 📊 Results and Reports

CSPBench generates comprehensive reports including:

- **Execution Results**: Algorithm performance metrics
- **Convergence Plots**: For iterative algorithms
- **Statistical Analysis**: Performance comparisons
- **Parameter Sensitivity**: Impact of parameter variations
- **Export Formats**: JSON, CSV, and TXT outputs

## 🏗️ Architecture

CSPBench follows hexagonal architecture principles:

```
src/
├── domain/           # Core business logic
├── application/      # Use cases and services
├── infrastructure/   # External concerns
└── presentation/     # CLI and user interfaces

algorithms/           # Algorithm implementations
├── baseline/
├── blf_ga/
├── csc/
├── h3_csp/
└── dp_csp/

deploy/               # Deployment configurations
├── docker-compose.yml
├── cloudbuild.yaml
├── deploy-cloud.sh
└── k8s/

docs/                 # Documentation
config/               # Configuration files
datasets/             # Sample datasets
batches/              # Batch configurations
```

## 🧪 Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=src

# Run specific test categories
pytest tests/unit/
pytest tests/integration/
```

## 🐳 Docker

CSPBench includes production-ready Docker configuration:

### Features
- Multi-stage build for optimization
- Non-root user for security
- Health checks
- Resource limits
- Development and production configurations

### Commands
```bash
# Build and run locally
make docker
make docker-run

# Development with hot reload
make run-web

# Production deployment
make deploy

# Docker Compose
make compose-up
make compose-down
```

For detailed deployment instructions, see [`deploy/README.md`](deploy/README.md).

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Install development dependencies
pip install -r requirements.txt
pip install -e .

# Run pre-commit hooks
pre-commit install

# Run linting and formatting
black .
ruff check .
mypy src/
```

## 📝 Citation

If you use CSPBench in your research, please cite:

```bibtex
@software{cspbench2025,
  author = {Grosmann, Diego},
  title = {CSPBench: A Framework for Closest String Problem Algorithms},
  year = {2025},
  url = {https://github.com/diegogrosmann/CSPBench},
  version = {0.1.0},
  license = {CC-BY-NC-SA-4.0}
}
```

## 📄 License

This project is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0) - see the [LICENSE](LICENSE) file for details.

**Key license terms:**
- ✅ **Share and adapt** freely for non-commercial purposes
- ✅ **Attribution required** - credit the original authors
- ✅ **ShareAlike** - derivative works must use the same license
- ❌ **No commercial use** without explicit permission

This license promotes open science while protecting against unauthorized commercial exploitation.

## 🙏 Acknowledgments

- Thanks to the bioinformatics community for valuable feedback
- Inspired by various CSP algorithm implementations
- Built with modern Python practices and tools

## 📧 Contact

- **Author**: Diego Grosmann
- **Email**: diego.grosmann@example.com
- **GitHub**: [@diegogrosmann](https://github.com/diegogrosmann)

---

For more detailed documentation, please visit our [Wiki](https://github.com/diegogrosmann/CSPBench/wiki).
