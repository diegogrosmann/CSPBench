# CSPBench

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17136438.svg)](https://doi.org/10.5281/zenodo.17136438)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

A framework for benchmarking Closest String Problem (CSP) algorithms, with CLI, web interface, and experiment templates.

---

## 🔎 What is CSP?

The Closest String Problem seeks a "center" string that minimizes the maximum distance to all input strings. Applications:
- Motif discovery in DNA
- Consensus sequence identification
- Error correction in biological data

---

## 🚀 Installation

```bash
git clone https://github.com/diegogrosmann/CSPBench.git
cd CSPBench
python3 -m venv venv
source venv/bin/activate
pip install -e .
cp .env.example .env
```

---

## 🧭 Usage

### CLI
```bash
# Run batch experiments
python main.py batch examples/batches/experiment_example.yaml

# Interactive menu
python main.py
```

### Web Interface (Beta)
> **Note:** The web service is in beta and may be unstable.

```bash
# Start the web server
python main.py web
# Access: http://localhost:8000
```

---

## 📦 Batch Templates

Templates available in `examples/batches/`:
- experiment_example.yaml    (Standard experiment)
- optimizationt_example.yaml (Parameter optimization)
- sensitivityt_example.yaml  (Sensitivity analysis)
- TEMPLATE.yaml              (Base template)

Example:
```bash
cp examples/batches/TEMPLATE.yaml my_experiment.yaml
# Edit my_experiment.yaml
python main.py batch my_experiment.yaml
```

---

## 🧠 Algorithms

- Baseline (Greedy consensus)
- BLF-GA (Genetic Algorithm)
- CSC (Core String Clustering)
- DP-CSP (Dynamic Programming)
- H2-CSP (Hybrid Heuristic)

- Guide to add your own: see `algorithms/README.md` (step‑by‑step on creating a new algorithm and integrating via `@register_algorithm`).

---

## 🧑‍💻 Development

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Lint and format
ruff check .
ruff format .

# (Optional) pre-commit
pre-commit install
```

### Make (if available)
> **Note:** This section is in beta and may contain errors.

If you have `make` installed, you can use the following commands for common development tasks:
```bash
make install        # Install dependencies
make test           # Run tests
make lint           # Check code quality
make run-local      # Run CLI
make run-web        # Run web interface
make docker-dev     # Build development image
make docker-cloud   # Build production image
make deploy-local   # Local Docker Compose
make help-cloud     # Cloud deployment options
```

---

## 🗂️ Structure (summary)

```
.
├── main.py
├── src/
│   ├── application/      # Use cases and orchestration
│   ├── domain/           # Core entities and business rules
│   ├── infrastructure/   # External adapters
│   └── presentation/     # CLI and Web Interface
├── algorithms/           # CSP algorithm implementations
├── config/               # Application configurations
├── examples/             # Dataset and batch examples
├── deployment/           # Deployment configurations
│   ├── cloud-run/        # Google Cloud Run
│   └── docker/           # Local development
└── tests/                # Automated tests
```

Documentation:
- `algorithms/README.md` — How to add a new algorithm (plugin guide)

---

## ⚙️ Configuration

Create `.env` from `.env.example`:
```bash
# External APIs
NCBI_EMAIL=your-email@example.com
NCBI_API_KEY=your-ncbi-key

# Directories
DATASET_DIRECTORY=./datasets
OUTPUT_BASE_DIRECTORY=./outputs
LOG_DIRECTORY=./logs

# Web Interface
WEB_HOST=0.0.0.0
PORT=8000
```

---

## ☁️ Deployment

### Development (local)
```bash
make deploy-local
```
---

## 🤝 Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidelines.

---

## 📝 Citation

If you use CSPBench in your research, please cite:
```bibtex
@software{grosmann2025cspbench,
  title={CSPBench: A Comprehensive Framework for Closest String Problem Benchmarking},
  author={Grosmann, Diego},
  year={2025},
  url={https://github.com/diegogrosmann/CSPBench}
}
```

---

## 📜 License

MIT License — see [LICENSE](LICENSE).

---

## 🙏 Acknowledgments

- NCBI (sequence data)
- Computational biology community
- Project contributors and users

---

## 📞 Support

- Issues: https://github.com/diegogrosmann/CSPBench/issues
- Discussions: https://github.com/diegogrosmann/CSPBench/discussions
- Email: diego.grosmann@ifma.edu.br

---

CSPBench — Making bioinformatics research more reproducible. 🧬
