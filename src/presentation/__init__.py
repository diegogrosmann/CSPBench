"""
CSPBench Presentation Layer - User Interface Components.

This module contains all user interface implementations and presentation logic
for the CSPBench framework. It provides multiple interface options designed
for different use cases and deployment scenarios while maintaining consistent
functionality across all interfaces.

Interface Components:
    - CLI (Command Line Interface): Optimized for automated scripts, CI/CD pipelines,
      and batch processing scenarios with comprehensive argument parsing
    - TUI (Terminal User Interface): Interactive terminal-based interface for 
      developers and researchers who prefer console-based workflows
    - Web Interface: Modern, responsive web-based UI with real-time monitoring,
      WebSocket support, and comprehensive visualization capabilities

Architecture Principles:
    The presentation layer follows hexagonal architecture principles:
    
    - Controllers: Handle user input validation and coordinate with application services
    - Views: Render data presentation and manage user interface state transitions
    - Routes: Define web endpoints, request/response handling, and API documentation
    - WebSocket: Enable real-time communication for live progress monitoring
    - Models: Pydantic models for request/response validation and serialization

Key Features:
    - Unified service interfaces across all presentation modes
    - Real-time progress monitoring with WebSocket connections
    - Comprehensive error handling and user feedback
    - Security validation and input sanitization
    - Responsive design for various screen sizes and devices
    - Export capabilities for results in multiple formats

Security Considerations:
    - Input validation through Pydantic models and custom validators
    - Path sanitization for file operations to prevent traversal attacks
    - Rate limiting and request validation for web endpoints
    - CORS configuration for secure cross-origin requests
    - File upload restrictions and content validation

Development Status:
    The presentation layer is actively developed as part of the migration to
    hexagonal architecture. Current focus is on the web interface with full
    WebSocket support for real-time monitoring capabilities.

Example Usage:
    Web Interface::
    
        from src.presentation.web.app import app
        import uvicorn
        
        # Start development server
        uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)
        
    CLI Interface (planned)::
    
        from src.presentation.cli import CLIApplication
        
        cli = CLIApplication()
        cli.run(['benchmark', '--config', 'batch.yaml'])

Future Enhancements:
    - Complete CLI implementation with argparse integration
    - Advanced TUI with rich text interface and interactive widgets
    - Mobile-responsive web interface improvements
    - API documentation generation and interactive testing
    - Plugin system for custom visualization components
"""

# TODO: Implement complete presentation layer components after core migration
# Current implementation focuses on web interface with comprehensive WebSocket support
# CLI and TUI implementations are planned for future releases

__all__ = [
    # Will be populated as components are implemented
    # "CLIApplication",
    # "TUIApplication", 
    # "WebApplication"
]
