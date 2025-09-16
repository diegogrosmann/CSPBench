# CSPBench# CSPBench



[![CI](https://github.com/diegogrosmann/CSPBench/workflows/CI/badge.svg)](https://github.com/diegogrosmann/CSPBench/actions)[![CI](https://github.com/diegogrosmann/CSPBench/workflows/CI/badge.svg)](https://github.com/diegogrosmann/CSPBench/actions)

[![Coverage](https://codecov.io/gh/diegogrosmann/CSPBench/branch/main/graph/badge.svg)](https://codecov.io/gh/diegogrosmann/CSPBench)[![Coverage](https://codecov.io/gh/diegogrosmann/CSPBench/branch/main/graph/badge.svg)](https://codecov.io/gh/diegogrosmann/CSPBench)

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)



A framework for benchmarking Closest String Problem (CSP) algorithms.A framework for benchmarking Closest String Problem (CSP) algorithms.



## What is CSP?## What is CSP?



The Closest String Problem finds a center string that minimizes the maximum distance to all input strings. Applications include:The Closest String Problem finds a center string that minimizes the maximum distance to all input strings. Applications include:

- DNA motif discovery- DNA motif discovery

- Consensus sequence identification  - Consensus sequence identification  

- Error correction in biological data- Error correction in biological data



## Installation## Installation



```bash```bash

git clone https://github.com/diegogrosmann/CSPBench.gitgit clone https://github.com/diegogrosmann/CSPBench.git

cd CSPBenchcd CSPBench

pip install -e .pip install -e .

``````



## Usage## Usage



### Command Line Interface### Command Line Interface



```bash```bash

# Run single algorithm# Run single algorithm

python main.py run baseline examples/datasets/example.fastapython main.py run baseline examples/datasets/example.fasta



# Run batch experiments  # Run batch experiments  

python main.py batch examples/batches/experiment_example.yamlpython main.py batch examples/batches/experiment_example.yaml



# Interactive menu# Interactive menu

python main.pypython main.py

``````



### Web Interface### Web Interface



```bash```bash

# Start web server# Start web server

python main.py webpython main.py web



# Open browser to http://localhost:8000# Open browser to http://localhost:8000

``````



## Batch Templates## Batch Templates



Use predefined templates for experiments:Use predefined templates for experiments:



```bash```bash

# Available templates in examples/batches/:# Available templates in examples/batches/:

# - experiment_example.yaml    # Standard experiment  # - experiment_example.yaml    # Standard experiment  

# - optimizationt_example.yaml # Parameter optimization# - optimizationt_example.yaml # Parameter optimization

# - sensitivityt_example.yaml  # Sensitivity analysis# - sensitivityt_example.yaml  # Sensitivity analysis

# - TEMPLATE.yaml              # Base template# - TEMPLATE.yaml              # Base template

``````



Copy a template, modify it, and run:Copy a template, modify it, and run:

```bash```bash

cp examples/batches/experiment_example.yaml my_experiment.yamlcp examples/batches/experiment_example.yaml my_experiment.yaml

# Edit my_experiment.yaml# Edit my_experiment.yaml

python main.py batch my_experiment.yamlpython main.py batch my_experiment.yaml

``````



## Algorithms## Algorithms



| Algorithm | Type | Reference || Algorithm | Type | Reference |

|-----------|------|-----------||-----------|------|-----------|

| **Baseline** | Brute Force | - || **Baseline** | Brute Force | - |

| **BLF-GA** | Genetic Algorithm | [Blum & Lozano, 2005] || **BLF-GA** | Genetic Algorithm | [Blum & Lozano, 2005] |

| **CSC** | Core String Clustering | [Custom] || **CSC** | Core String Clustering | [Custom] |

| **DP-CSP** | Dynamic Programming | [Custom] || **DP-CSP** | Dynamic Programming | [Custom] |

| **H2-CSP** | Hybrid Heuristic | [Custom] || **H2-CSP** | Hybrid Heuristic | [Custom] |

## Development

## Development

```bash

```bash# Install dependencies

# Install dependenciespip install -e ".[dev]"

pip install -e ".[dev]"

# Run tests

# Run testspytest

pytest

# Format code

# Format coderuff format .

ruff format .```

```

## Citation

## Citation

```bibtex

```bibtex@software{grosmann2025cspbench,

@software{grosmann2025cspbench,  title={CSPBench: A Framework for Closest String Problem Benchmarking},

  title={CSPBench: A Framework for Closest String Problem Benchmarking},  author={Grosmann, Diego},

  author={Grosmann, Diego},  year={2025},

  year={2025},  url={https://github.com/diegogrosmann/CSPBench}

  url={https://github.com/diegogrosmann/CSPBench}}

}```

```

## License

## License

MIT License - see [LICENSE](LICENSE) file.

MIT License - see [LICENSE](LICENSE) file.```

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
