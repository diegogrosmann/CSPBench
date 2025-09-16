# CSPBench - Closest String Problem Benchmark

[![CI](https://github.com/diegogrosmann/CSPBench/workflows/CI/badge.svg)](https://github.com/diegogrosmann/CSPBench/actions)
[![Coverage](https://codecov.io/gh/diegogrosmann/CSPBench/branch/main/graph/badge.svg)](https://codecov.io/gh/diegogrosmann/CSPBench)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXX)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A comprehensive framework for benchmarking and analyzing Closest String Problem (CSP) algorithms with clean architecture design and research-grade reproducibility.

## 🔬 Research Context

The Closest String Problem is a fundamental NP-hard optimization problem in computational biology and theoretical computer science. Given a set of strings, the goal is to find a center string that minimizes the maximum Hamming distance to all input strings. This problem has applications in:

- **Motif Discovery**: Finding conserved regions in DNA/protein sequences
- **Consensus Sequence Identification**: Determining representative sequences
- **Error Correction**: Correcting sequencing errors in biological data
- **Pattern Recognition**: Identifying common patterns in string data

## ✨ Key Features

- **🏗️ Clean Architecture**: Hexagonal architecture for maintainability and testability
- **📊 Comprehensive Benchmarking**: Compare multiple CSP algorithms systematically
- **🧬 Bioinformatics Integration**: Native support for FASTA files and NCBI data
- **🌐 Web Interface**: User-friendly dashboard for experiment management
- **📈 Performance Monitoring**: Real-time progress tracking and resource monitoring
- **🔄 Reproducible Research**: Automated experiment configuration and result archival
- **📦 Containerized Deployment**: Docker and cloud-ready configurations
- **🔧 Extensible Design**: Easy algorithm integration and custom metrics

## 🚀 Quick Start

### Installation

#### Option 1: pip (Recommended)
```bash
pip install cspbench
```

#### Option 2: From Source
```bash
git clone https://github.com/diegogrosmann/CSPBench.git
cd CSPBench
pip install -e .
```

#### Option 3: Docker
```bash
docker run -p 8000:8000 diegogrosmann/cspbench:latest
```

### Basic Usage

#### Command Line Interface
```bash
# Run a simple benchmark
cspbench run --algorithm baseline --dataset examples/data.fasta

# Web interface
cspbench web --port 8000

# Generate synthetic datasets
cspbench dataset generate --sequences 100 --length 50 --alphabet DNA
```

#### Python API
```python
from cspbench import CSPBench
from cspbench.algorithms import BaselineCSP
from cspbench.datasets import load_fasta

# Load data
sequences = load_fasta("path/to/sequences.fasta")

# Initialize algorithm
algorithm = BaselineCSP(max_distance=3)

# Run benchmark
benchmark = CSPBench()
results = benchmark.run(algorithm, sequences)

print(f"Best solution: {results.center_string}")
print(f"Max distance: {results.max_distance}")
print(f"Runtime: {results.execution_time:.2f}s")
```

## 📚 Documentation

- **[User Guide](https://diegogrosmann.github.io/CSPBench/user-guide/)**: Comprehensive usage instructions
- **[API Reference](https://diegogrosmann.github.io/CSPBench/api/)**: Complete API documentation
- **[Algorithm Guide](https://diegogrosmann.github.io/CSPBench/algorithms/)**: Implemented algorithms overview
- **[Developer Guide](https://diegogrosmann.github.io/CSPBench/contributing/)**: Contributing and extending CSPBench

## 🧪 Included Algorithms

| Algorithm | Type | Time Complexity | Space Complexity | Reference |
|-----------|------|-----------------|------------------|-----------|
| **Baseline** | Brute Force | O(k^l \* n \* l) | O(l) | - |
| **BLF-GA** | Genetic Algorithm | O(g \* p \* n \* l) | O(p \* l) | [Blum & Lozano, 2005] |
| **CSC** | Core String Clustering | O(n² \* l) | O(n \* l) | [Custom] |
| **DP-CSP** | Dynamic Programming | O(n \* l \* d) | O(l \* d) | [Custom] |
| **H2-CSP** | Hybrid Heuristic | O(n \* l \* log(l)) | O(n \* l) | [Custom] |

## 📊 Example Results

```python
# Benchmark comparison on synthetic data
results = benchmark.compare_algorithms(
    algorithms=[BaselineCSP(), BLFGA(), CSC(), DPCSP()],
    dataset="synthetic_100x20",
    metrics=["solution_quality", "runtime", "memory_usage"]
)

# Results automatically saved to outputs/benchmark_TIMESTAMP/
```

## 🏗️ Architecture

CSPBench follows hexagonal architecture principles:

```
src/
├── domain/              # Core business logic
│   ├── entities/        # Domain entities (Dataset, Algorithm, Result)
│   ├── repositories/    # Repository interfaces
│   └── services/        # Domain services
├── application/         # Use cases and application services
│   ├── services/        # Application services
│   └── use_cases/       # Business use cases
├── infrastructure/      # External concerns
│   ├── persistence/     # Database implementations
│   ├── algorithms/      # Algorithm implementations
│   └── external/        # External service integrations
└── presentation/        # User interfaces
    ├── cli/             # Command-line interface
    ├── web/             # Web interface
    └── api/             # REST API
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone repository
git clone https://github.com/diegogrosmann/CSPBench.git
cd CSPBench

# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run linting
ruff check .
ruff format .

# Install pre-commit hooks
pre-commit install
```

## � Citation

If you use CSPBench in your research, please cite:

```bibtex
@software{grosmann2025cspbench,
  title={CSPBench: A Comprehensive Framework for Closest String Problem Benchmarking},
  author={Grosmann, Diego},
  year={2025},
  url={https://github.com/diegogrosmann/CSPBench},
  doi={10.5281/zenodo.XXXXX}
}
```

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- National Center for Biotechnology Information (NCBI) for sequence data access
- The computational biology community for algorithm implementations and datasets
- Contributors and users who have helped improve this framework

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/diegogrosmann/CSPBench/issues)
- **Discussions**: [GitHub Discussions](https://github.com/diegogrosmann/CSPBench/discussions)
- **Email**: diego.grosmann@example.com

---

**CSPBench** - Making bioinformatics research more reproducible, one benchmark at a time. 🧬
│   ├── domain/            # Entidades e regras de negócio
│   ├── infrastructure/    # Adaptadores externos
│   └── presentation/      # CLI e Web Interface
├── algorithms/            # Implementações de algoritmos CSP
├── config/               # Configurações da aplicação
├── examples/             # Exemplos de datasets e batches
├── deployment/           # 🆕 Configurações de deployment
│   ├── cloud-run/       # Google Cloud Run
│   └── docker/          # Desenvolvimento local
└── tests/               # Testes automatizados
```

## 🛠️ Comandos Disponíveis

### Desenvolvimento
```bash
make install        # Instalar dependências
make test          # Executar testes
make lint          # Verificar qualidade do código
make run-local     # Executar CLI
make run-web       # Executar interface web
```

### Docker & Deployment
```bash
make docker-dev             # Build imagem de desenvolvimento
make docker-cloud           # Build imagem de produção
make deploy-local           # Docker Compose local
make deploy-cloud-run-secure PROJECT_ID=meu-projeto  # Deploy seguro no Cloud Run
make help-cloud             # Ver todas as opções de cloud deployment
```

## 📚 Documentação

- **[Deployment Guide](deployment/README.md)** - Configurações de deployment
- **[Cloud Run](deployment/cloud-run/README.md)** - Deploy específico para Cloud Run
- **[Web Interface](src/presentation/web/README.md)** - Interface web
- **[Algorithms](algorithms/README.md)** - Algoritmos disponíveis

## 🧬 Exemplos

### CLI
```bash
# Executar algoritmo específico
python main.py run baseline examples/datasets/example.fasta

# Executar batch de experimentos
python main.py batch examples/batches/experiment_example.yaml

# Iniciar menu interativo
python main.py
```

### Web Interface
Acesse `http://localhost:8000` após executar:
```bash
make run-web
# ou
python main.py web
```

## 🔧 Configuração

### Variáveis de Ambiente
Copie `.env.example` para `.env` e configure:
```bash
# APIs externas
NCBI_EMAIL=seu-email@exemplo.com
NCBI_API_KEY=sua-chave-ncbi

# Diretórios
DATASET_DIRECTORY=./datasets
OUTPUT_BASE_DIRECTORY=./outputs
LOG_DIRECTORY=./logs

# Web Interface
WEB_HOST=0.0.0.0
PORT=8000
```

## 🏗️ Arquitetura

O CSPBench usa **Arquitetura Hexagonal** (Ports & Adapters):

- **Domain**: Entidades core e regras de negócio
- **Application**: Casos de uso e orquestração
- **Infrastructure**: Adaptadores para I/O, persistência
- **Presentation**: CLI, Web API, interfaces externas

## 📦 Deploy

### Para desenvolvimento:
```bash
make deploy-local
```

### Para produção (Cloud Run):
```bash
# Opção 1: Deploy seguro (recomendado)
make deploy-cloud-run-secure PROJECT_ID=meu-projeto

# Opção 2: Com variáveis de ambiente
export NCBI_EMAIL="seu-email@exemplo.com"
export NCBI_API_KEY="sua-chave"
make deploy-cloud-run-env PROJECT_ID=meu-projeto

# Ver todas as opções
make help-cloud
```

Veja mais detalhes em [deployment/README.md](deployment/README.md).

## 🤝 Contribuindo

Veja [CONTRIBUTING.md](CONTRIBUTING.md) para diretrizes de contribuição.

## 📄 Licença

Este projeto está sob a licença especificada em [LICENSE](LICENSE).
