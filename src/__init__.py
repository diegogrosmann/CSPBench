"""
CSPBench - Framework for Closest String Problem Algorithm Benchmarking.

This package provides a comprehensive, modular framework for testing and comparing
algorithms that solve the Closest String Problem (CSP). The implementation follows
hexagonal architecture principles with clear separation of concerns across
domain, application, infrastructure, and presentation layers.

Features:
    - Multiple CSP algorithm implementations with standardized interfaces
    - Comprehensive benchmarking and comparison tools with statistical analysis
    - Modern web-based interface with real-time execution monitoring
    - Real-time progress tracking via WebSocket connections
    - Advanced dataset management and generation utilities
    - Interactive results visualization and export capabilities
    - Batch processing with YAML configuration support
    - Docker containerization for consistent deployment

Architecture:
    The framework follows hexagonal (ports and adapters) architecture:

    - Domain Layer: Core business logic, algorithm interfaces, and entities
    - Application Layer: Use cases, service orchestration, and workflow management
    - Infrastructure Layer: Data persistence, external adapters, and technical concerns
    - Presentation Layer: User interfaces (CLI, TUI, Web) and API endpoints

Example:
    Basic usage for programmatic access::

        from src.application.services.work_service import get_work_service

        # Create and execute a benchmarking work
        work_service = get_work_service()
        work_id = work_service.create_batch_work(
            datasets=['sample_dataset'],
            algorithms=['brute_force', 'heuristic'],
            tasks=['closest_string']
        )

        # Monitor progress via web interface or programmatically
        status = work_service.get_status(work_id)

        # Access results when completed
        results = work_service.get_results(work_id)

Web Interface:
    Start the web interface with::

        python -m src.presentation.web.run_web

    Or access via Docker::

        docker run -p 8000:8000 cspbench:latest

Version:
    0.1.0

Authors:
    CSPBench Development Team

License:
    This project is part of academic research. See LICENSE file for details.
"""

__version__ = "0.1.0"
__author__ = "CSPBench Development Team"

# Lazy import strategy: avoid importing subpackages automatically to prevent side-effects
__all__ = ["__version__", "__author__"]
