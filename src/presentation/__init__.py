"""
CSPBench Presentation Layer.

This module contains all user interface implementations and presentation logic
for the CSPBench framework. It provides multiple interface options for
different use cases and deployment scenarios.

Components:
    - CLI (Command Line Interface): For automated scripts and batch processing
    - TUI (Terminal User Interface): For interactive terminal-based usage
    - Web Interface: Modern web-based UI with real-time monitoring

The presentation layer is designed to be independent of business logic,
communicating with the application layer through well-defined service interfaces.
This separation allows for easy addition of new interface types without
affecting core functionality.

Architecture:
    - Controllers: Handle user input and coordinate with services
    - Views: Render data and manage user interface state
    - Routes: Define web endpoints and request handling
    - WebSocket: Real-time communication for progress monitoring

Note:
    This module is currently under active development as part of the
    migration to the new hexagonal architecture.
"""

# TODO: Implement complete presentation layer after migration
# Current implementation focuses on web interface with WebSocket support
__all__ = []
