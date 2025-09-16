"""Command Line Interface (CLI) Presentation Layer

This package provides the command-line interface components for CSPBench,
offering users a terminal-based way to interact with the framework's functionality.

The CLI layer includes:
- Interactive progress monitoring with real-time updates
- Batch file execution and management
- Dataset generation wizards
- Work item management commands
- Algorithm registry browsing

Key Components:
    commands.py: CLI command registration and implementation
    progress_monitor.py: Real-time htop-style progress monitoring
    interactive_menu.py: Interactive batch file selection interface
    dataset_wizard.py: Interactive dataset generation wizard

Features:
    - Comprehensive batch execution with progress tracking
    - Real-time monitoring interface using curses
    - Interactive dataset generation and management
    - Work item lifecycle management (pause, resume, restart)
    - Algorithm discovery and listing
    - Web interface launcher

Usage:
    The CLI is accessed through the main entry point and provides
    various commands for different aspects of the CSPBench workflow::

        # Execute batch configuration
        python main.py batch config.yaml

        # Monitor work progress
        python main.py monitor work_id

        # List available algorithms
        python main.py algorithms

        # Start web interface
        python main.py web

        # Generate datasets interactively
        python main.py datasetsave

Note:
    This CLI implementation is designed to work seamlessly with the
    application services layer and provides consistent behavior
    across different execution contexts.
"""

# TODO: Implement additional CLI utilities after core migration
__all__ = []
