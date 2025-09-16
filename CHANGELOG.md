# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2025-01-16

### Added
- Initial stable release of CSPBench framework
- Complete hexagonal architecture implementation
- Web interface with real-time monitoring
- Command-line interface with comprehensive commands
- Five CSP algorithms (Baseline, BLF-GA, CSC, DP-CSP, H2-CSP)
- FASTA file support and NCBI data integration
- Automated testing suite with >90% coverage
- Docker containerization support
- Cloud deployment configurations (Google Cloud Run)
- Comprehensive documentation and API reference
- Pre-commit hooks and code quality tools
- MIT license for open-source distribution

### Features
- **Algorithm Implementations**:
  - Baseline brute-force algorithm
  - BLF-GA (Genetic Algorithm approach)
  - CSC (Core String Clustering)
  - DP-CSP (Dynamic Programming solution)
  - H2-CSP (Hybrid Heuristic approach)

- **User Interfaces**:
  - Interactive web dashboard
  - Command-line interface
  - REST API endpoints
  - Real-time WebSocket monitoring

- **Data Management**:
  - FASTA file import/export
  - NCBI sequence database integration
  - Synthetic dataset generation
  - Result archival and comparison

- **Development Tools**:
  - Comprehensive test suite
  - Code quality enforcement (ruff, mypy, black)
  - Pre-commit hooks
  - Automated CI/CD pipeline
  - Docker development environment

### Infrastructure
- Clean hexagonal architecture design
- SQLAlchemy-based persistence layer
- FastAPI web framework
- Typer CLI framework
- WebSocket real-time communication
- Comprehensive logging system
- Configuration management
- Error handling and recovery

### Documentation
- Complete API documentation
- User guide and tutorials
- Algorithm implementation details
- Development and contribution guidelines
- Deployment instructions
- Research context and citations

## [0.1.0] - 2024-12-01

### Added
- Initial project structure
- Basic algorithm implementations
- Core domain models
- Simple CLI interface

---

## Release Notes

### v1.0.0 - "Stable Foundation"

This is the first stable release of CSPBench, marking a significant milestone in the development of a comprehensive framework for Closest String Problem research and benchmarking.

**Key Highlights:**
- **Research-Ready**: Implements five different CSP algorithms with proper benchmarking capabilities
- **User-Friendly**: Multiple interfaces (web, CLI, API) for different use cases
- **Developer-Friendly**: Clean architecture, comprehensive tests, and excellent documentation
- **Production-Ready**: Docker support, cloud deployment, and robust error handling
- **Open Source**: MIT license enabling wide adoption and contribution

**For Researchers:**
- Reproducible experiment configurations
- Automated result archival and comparison
- Support for biological sequence data (FASTA, NCBI)
- Comprehensive performance metrics

**For Developers:**
- Extensible architecture for adding new algorithms
- Well-documented API for integration
- Modern Python practices and tooling
- Comprehensive test coverage

**Next Steps:**
- Community feedback and contributions
- Additional algorithm implementations
- Enhanced visualization capabilities
- Performance optimizations

---

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/diegogrosmann/CSPBench/tags).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.