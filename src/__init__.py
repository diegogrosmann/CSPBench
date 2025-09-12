"""
CSPBench - Framework for Closest String Problem Algorithm Benchmarking.

This package provides a modular framework for testing and comparing algorithms
that solve the Closest String Problem (CSP). The implementation follows
hexagonal architecture principles with clear separation of concerns across
domain, application, infrastructure, and presentation layers.

Features:
    - Multiple CSP algorithm implementations
    - Comprehensive benchmarking and comparison tools
    - Web-based interface for execution monitoring
    - Real-time progress tracking via WebSocket
    - Dataset management and generation utilities
    - Results visualization and export capabilities

Architecture:
    - Domain Layer: Core business logic and algorithm interfaces
    - Application Layer: Use cases and service orchestration
    - Infrastructure Layer: Data persistence and external adapters
    - Presentation Layer: User interfaces (CLI, TUI, Web)

Example:
    Basic usage example::

        from src.application.services.work_service import get_work_service
        
        work_service = get_work_service()
        work_id = work_service.create_batch_work(
            datasets=['sample_dataset'],
            algorithms=['brute_force', 'heuristic'],
            tasks=['closest_string']
        )
        
        # Monitor progress via web interface or programmatically
        status = work_service.get_status(work_id)

Version:
    0.1.0

Authors:
    CSPBench Development Team
"""

__version__ = "0.1.0"
__author__ = "CSPBench Development Team"

# Lazy import strategy: avoid importing subpackages automatically to prevent side-effects
__all__ = ["__version__", "__author__"]
